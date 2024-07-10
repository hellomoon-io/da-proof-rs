use anyhow::{anyhow, ensure};
use bisetmap::BisetMap;
use grid::Grid;
use itertools::Itertools;
use num::integer::sqrt;
use reed_solomon_erasure::galois_16::ReedSolomon;
use rs_merkle::algorithms::Sha256;
use rs_merkle::{Hasher, MerkleTree};

use crate::util::Axis;

/// Type aliases

pub type Share<const SHARE_SIZE: usize> = [u8; SHARE_SIZE];
pub type MerkleRoot = [u8; 32];
type Square<const SHARE_SIZE: usize> = Grid<Share<SHARE_SIZE>>;

/// A "square" of data aligned into `Share`s.
///
/// *Invariant*: the data is square, i.e. the grid is equally as wide as it is high.
#[derive(Debug, Clone)]
pub struct DataSquare<const SHARE_SIZE: usize>(Square<SHARE_SIZE>);

/// An extended "square" of data.
///
/// This includes the erasure encoded data derived from the original `DataSquare`.
#[derive(Debug, Clone)]
pub struct ExtendedDataSquare<const SHARE_SIZE: usize> {
    original_width: usize,
    square: DataSquare<SHARE_SIZE>,
    row_roots: Vec<MerkleRoot>,
    col_roots: Vec<MerkleRoot>,
}

/// An error that is returned when Byzantine (invalid, malicious, etc.) data is detected.
///
/// This can be a result of the shares in an axis not matching the known Merkle roots,
/// or recovery data that doesn't match the original.
#[derive(Debug)]
pub struct ByzantineDataError<const SHARE_SIZE: usize> {
    /// The axis on which Byzantine data was detected
    axis: Axis,
    /// The index of the axis on which Byzantine data was detected
    idx: usize,
    /// The shares that could be used to construct a fraud proof.
    shares: Option<Vec<Share<SHARE_SIZE>>>,
}

impl<const SHARE_SIZE: usize> ByzantineDataError<SHARE_SIZE> {
    /// Constructs a new `ByzantineDataError` with `shares` populated.
    fn new_with_shares(eds: &ExtendedDataSquare<SHARE_SIZE>, axis: Axis, idx: usize) -> Self {
        Self {
            axis,
            idx,
            shares: Some(
                (match axis {
                    Axis::Row => eds.square.0.iter_row(idx),
                    Axis::Col => eds.square.0.iter_col(idx),
                })
                .map(ToOwned::to_owned)
                .collect_vec(),
            ),
        }
    }
}

impl<const SHARE_SIZE: usize> std::fmt::Display for ByzantineDataError<SHARE_SIZE> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "byzantine data detected @ {} {}", self.axis, self.idx)
    }
}

impl<const SHARE_SIZE: usize> std::error::Error for ByzantineDataError<SHARE_SIZE> {}

/// Computes the root of the Merkle tree with the leaves as shares in the square.
fn compute_root<'a, const SHARE_SIZE: usize>(
    square: &Square<SHARE_SIZE>,
    axis: Axis,
    idx: usize,
) -> anyhow::Result<MerkleRoot> {
    let iter = match axis {
        Axis::Row => square.iter_row(idx),
        Axis::Col => square.iter_col(idx),
    };

    let leaves: Vec<MerkleRoot> = iter.map(|share| Sha256::hash(share)).collect();
    let tree = MerkleTree::<Sha256>::from_leaves(&leaves);

    tree.root()
        .ok_or(anyhow!("unable to calculate root for {axis} {idx}"))
}

impl<const SHARE_SIZE: usize> DataSquare<SHARE_SIZE> {
    /// Constructs a `DataSquare` from a flat slice containing `Share`s.
    ///
    /// Returns `Err(_)` if the length of the input is not a perfect square.
    pub fn new(data: &[Share<SHARE_SIZE>]) -> anyhow::Result<Self> {
        let width = sqrt(data.len());

        ensure!(
            width * width == data.len(),
            "data length must be a perfect square, received {}",
            data.len()
        );

        let mut inner: Square<SHARE_SIZE> = Grid::init(width, width, [0u8; SHARE_SIZE]);

        (0..width).for_each(|row_idx| {
            inner
                .iter_row_mut(row_idx)
                .enumerate()
                .for_each(|(i, elt)| *elt = data[row_idx * width + i]);
        });

        Ok(Self(inner))
    }

    /// Erasure encodes `self` using `width` rows & columns.
    ///
    /// Uses a Reed-Solomon encoding w/ a 16-bit galois field.
    fn erasure_encode(&mut self, width: usize) -> anyhow::Result<()> {
        let rs = ReedSolomon::new(width, width)?;

        // Encode each *original* row and column.
        //
        // They don't depend on each other, so we can do it in the same loop.
        (0..width).try_for_each(|i| {
            {
                // Annoyingly, the Reed-Solomon library we're using here only supports encoding
                // byte arrays of size 2 (i.e. a 16 bit integer). This means we need to "compress"
                // and "decompress" the data during encoding and decoding.
                //
                // TODO: figure out a (potentially unsafe) way to do this without so much copying.
                let mut row_vec = self
                    .0
                    .iter_row(i)
                    .map(|share| {
                        share
                            .iter()
                            .tuples()
                            .map(|(b0, b1)| [*b0, *b1])
                            .collect_vec()
                    })
                    .collect_vec();

                rs.encode(&mut row_vec)?;

                self.0
                    .iter_row_mut(i)
                    .zip(row_vec.iter())
                    .for_each(|(dest, src)| {
                        dest.iter_mut().enumerate().for_each(|(i, b)| {
                            *b = src[i / 2][i % 2];
                        })
                    });
            };
            // TODO: maybe move this to a separate function, the reptition is kinda ugly.
            {
                let mut col_vec = self
                    .0
                    .iter_col(i)
                    .map(|share| {
                        share
                            .iter()
                            .tuples()
                            .map(|(b0, b1)| [*b0, *b1])
                            .collect_vec()
                    })
                    .collect_vec();

                rs.encode(&mut col_vec)?;

                self.0
                    .iter_col_mut(i)
                    .zip(col_vec.iter())
                    .for_each(|(dest, src)| {
                        dest.iter_mut().enumerate().for_each(|(i, b)| {
                            *b = src[i / 2][i % 2];
                        })
                    });
            };
            Ok::<(), anyhow::Error>(())
        })?;

        // Encode the bottom right quadrant of the array using the rows from the column
        // encoding of the original data.
        //
        // NOTE: this should be equivalent to encoding the columns of the row encoding
        // of the original data (this is checked when not in release mode).
        (width..width * 2).try_for_each(|i| {
            let mut row_vec = self
                .0
                .iter_row(i)
                .map(|share| {
                    share
                        .iter()
                        .tuples()
                        .map(|(b0, b1)| [*b0, *b1])
                        .collect_vec()
                })
                .collect_vec();
            rs.encode(&mut row_vec)?;
            self.0
                .iter_row_mut(i)
                .zip(row_vec.iter())
                .for_each(|(dest, src)| {
                    dest.iter_mut().enumerate().for_each(|(i, b)| {
                        *b = src[i / 2][i % 2];
                    })
                });
            Ok::<(), anyhow::Error>(())
        })?;

        // Invariant: encoding Q3 from the rows of Q2 should be equivalent to encoding it from the
        // columns of Q1.
        #[cfg(debug_assertions)]
        (width..width * 2).try_for_each(|i| {
            let mut col_vec = self
                .0
                .iter_col(i)
                .map(|share| {
                    share
                        .iter()
                        .tuples()
                        .map(|(b0, b1)| [*b0, *b1])
                        .collect_vec()
                })
                .collect_vec();
            rs.encode(&mut col_vec)?;
            self.0
                .iter_col_mut(i)
                .zip(col_vec.iter())
                .for_each(|(dest, src)| {
                    dest.iter_mut().enumerate().for_each(|(i, b)| {
                        assert_eq!(*b, src[i / 2][i % 2]);
                    })
                });
            Ok::<(), anyhow::Error>(())
        })?;

        Ok(())
    }

    /// Computes the roots for a particular axis and index in `self`.
    fn compute_roots(&self, axis: Axis) -> anyhow::Result<Vec<MerkleRoot>> {
        let dim = match axis {
            Axis::Row => self.0.rows(),
            Axis::Col => self.0.cols(),
        };
        (0..dim)
            .map(|i| compute_root(&self.0, axis, i))
            .collect::<anyhow::Result<Vec<MerkleRoot>>>()
    }

    /// Extends self using a 2D Reed-Solomon encoding.
    ///
    /// The resulting backing grid has dimension `2 * width` x `2 * width`.
    pub fn extend(mut self) -> anyhow::Result<ExtendedDataSquare<SHARE_SIZE>> {
        let original_width = self.0.cols();
        (0..self.0.cols()).for_each(|_| self.0.push_col(vec![[0u8; SHARE_SIZE]; self.0.rows()]));
        (0..self.0.rows()).for_each(|_| self.0.push_row(vec![[0u8; SHARE_SIZE]; self.0.cols()]));
        self.erasure_encode(original_width)?;
        let row_roots = self.compute_roots(Axis::Row)?;
        let col_roots = self.compute_roots(Axis::Col)?;
        Ok(ExtendedDataSquare {
            square: self,
            original_width,
            row_roots,
            col_roots,
        })
    }
}

/// Describes a step in the repair process.
#[derive(Debug)]
struct RepairStep {
    solved: bool,
    progress_made: bool,
}

impl Default for RepairStep {
    fn default() -> Self {
        Self {
            solved: true,
            progress_made: false,
        }
    }
}

impl RepairStep {
    fn combine(&self, other: &Self) -> Self {
        Self {
            solved: self.solved && other.solved,
            progress_made: self.progress_made || other.progress_made,
        }
    }
}

impl<const SHARE_SIZE: usize> ExtendedDataSquare<SHARE_SIZE> {
    /// TODO: verify encoding

    /// Before executing a repair, check if the slice is valid.
    ///
    /// No-op if the slice is incomplete (according to the user).
    fn repair_preflight_slice<'a>(
        &self,
        axis: Axis,
        idx: usize,
        missing: &BisetMap<usize, usize>,
    ) -> anyhow::Result<()> {
        let has_missing = match axis {
            Axis::Row => missing.key_exists(&idx),
            Axis::Col => missing.value_exists(&idx),
        };
        if !has_missing {
            let roots = match axis {
                Axis::Row => &self.row_roots,
                Axis::Col => &self.col_roots,
            };
            let expected_root = roots
                .get(idx)
                .ok_or(ByzantineDataError::new_with_shares(self, axis, idx))?;
            let calculated_root = compute_root(&self.square.0, axis, idx)?;
            if expected_root != &calculated_root {
                Err(ByzantineDataError::new_with_shares(self, axis, idx).into())
            } else {
                Ok(())
            }
        } else {
            Ok(())
        }
    }

    /// Before executing a repair, check if `self` contains any known Byzantine data.
    fn repair_preflight(&self, missing: &BisetMap<usize, usize>) -> anyhow::Result<()> {
        (0..self.original_width * 2)
            .try_for_each(|i| self.repair_preflight_slice(Axis::Row, i, &missing))?;
        (0..self.original_width * 2)
            .try_for_each(|i| self.repair_preflight_slice(Axis::Col, i, &missing))?;
        Ok(())
    }

    /// Attempt to solve the current slice.
    fn solve_slice<'a>(
        &mut self,
        axis: Axis,
        idx: usize,
        missing: &BisetMap<usize, usize>,
    ) -> anyhow::Result<RepairStep> {
        let rs = ReedSolomon::new(self.original_width, self.original_width)?;

        let (missing_indices, iter, expected_root) = match axis {
            Axis::Row => (
                missing.get(&idx),
                self.square.0.iter_row(idx),
                self.row_roots[idx],
            ),
            Axis::Col => (
                missing.rev_get(&idx),
                self.square.0.iter_col(idx),
                self.col_roots[idx],
            ),
        };

        if !missing_indices.is_empty() {
            let mut shares = iter
                .enumerate()
                .map(|(jdx, share)| {
                    if missing_indices.contains(&jdx) {
                        None
                    } else {
                        let temp = share
                            .iter()
                            .tuples()
                            .map(|(b0, b1)| [*b0, *b1])
                            .collect_vec();
                        Some(temp)
                    }
                })
                .collect_vec();

            if let Err(_) = rs.reconstruct(&mut shares) {
                return Ok(RepairStep {
                    solved: false,
                    progress_made: false,
                });
            };

            // TODO: so much copying...
            //
            // Ideally we use a different library for the RS encoding
            let decompressed = shares
                .iter()
                .map(|share| {
                    let reconstructed = share
                        .as_ref()
                        .ok_or(anyhow!("reconstruction succeeded but got None."))?;

                    let mut share = [0u8; SHARE_SIZE];
                    share.iter_mut().enumerate().for_each(|(i, dest)| {
                        *dest = reconstructed[i / 2][i % 2];
                    });

                    Ok::<Share<SHARE_SIZE>, anyhow::Error>(share)
                })
                .collect::<anyhow::Result<Vec<Share<SHARE_SIZE>>>>()?;

            let calculated_root = {
                let leaves = decompressed
                    .iter()
                    .map(|share| Sha256::hash(share))
                    .collect_vec();
                let tree = MerkleTree::<Sha256>::from_leaves(&leaves);
                tree.root()
                    .ok_or(ByzantineDataError::new_with_shares(self, axis, idx))
            }?;

            if calculated_root != expected_root {
                return Err(ByzantineDataError::new_with_shares(self, axis, idx).into());
            }

            missing_indices.iter().try_for_each(|jdx| {
                let (row_idx, col_idx) = match axis {
                    Axis::Row => (idx, *jdx),
                    Axis::Col => (*jdx, idx),
                };

                let share = unsafe { self.square.0.get_unchecked_mut(row_idx, col_idx) };
                *share = decompressed[*jdx];

                Ok::<(), anyhow::Error>(())
            })?;

            missing_indices.iter().for_each(|jdx| {
                // TODO: check orthogonal axis root if slice is now completed
                let (row_idx, col_idx) = match axis {
                    Axis::Row => (idx, *jdx),
                    Axis::Col => (*jdx, idx),
                };

                missing.remove(&row_idx, &col_idx)
            });

            Ok(RepairStep {
                solved: true,
                progress_made: true,
            })
        } else {
            Ok(RepairStep {
                solved: true,
                progress_made: false,
            })
        }
    }

    /// Attempt to solve the EDS by reaching a fixed point.
    ///
    /// Returns `Err(_)` if Byzantine data is detected or no progress can be made.
    fn solve(&mut self, missing: &BisetMap<usize, usize>) -> anyhow::Result<()> {
        loop {
            let repair_step: RepairStep =
                (0..self.original_width * 2).try_fold(RepairStep::default(), |acc, i| {
                    let row_repair_step = self.solve_slice(Axis::Row, i, missing)?;
                    let col_repair_step = self.solve_slice(Axis::Col, i, missing)?;

                    Ok::<_, anyhow::Error>(acc.combine(&row_repair_step).combine(&col_repair_step))
                })?;

            if repair_step.solved {
                break;
            } else if !repair_step.progress_made {
                return Err(anyhow!("no progress made in solution"));
            }
        }

        Ok(())
    }

    /// Repair the EDS, assuming the indices in `missing` are erased/corrupted.
    pub fn repair(&mut self, missing: &[(usize, usize)]) -> anyhow::Result<()> {
        // Assemble missing maps
        let missing_map = BisetMap::new();
        missing.iter().for_each(|(row, col)| {
            missing_map.insert(*row, *col);
        });
        self.repair_preflight(&missing_map)?;
        self.solve(&missing_map)?;
        Ok(())
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use grid::grid;

    #[test]
    fn two_by_two() {
        let data = [[0u8, 1], [2, 3], [4, 5], [6, 7]];
        let ds = DataSquare::new(&data).unwrap();
        assert_eq!(
            ds.0,
            grid![
                [[0u8, 1], [2, 3]]
                [[4, 5], [6, 7]]
            ]
        );
    }

    #[test]
    fn extension() {
        let arr: [[u8; 64]; 4] = (0..4)
            .map(|i| [i + 1; 64])
            .collect::<Vec<_>>()
            .try_into()
            .unwrap();
        let ds = DataSquare::new(&arr).unwrap();
        let _eds = ds.extend(); // The debug assert should fail if this test doesn't pass
    }

    #[test]
    fn extension_repair_max() {
        let arr: [[u8; 64]; 4] = (0..4)
            .map(|i| [i + 1; 64])
            .collect::<Vec<_>>()
            .try_into()
            .unwrap();
        let ds = DataSquare::new(&arr).unwrap();
        let eds0 = ds.extend().unwrap();
        let mut eds1 = eds0.clone();
        let missing = [
            (0, 0),
            (0, 2),
            (0, 3),
            (1, 0),
            (1, 1),
            (1, 2),
            (1, 3),
            (2, 0),
            (2, 1),
            (2, 2),
            (3, 0),
            (3, 3),
        ];

        // Zero out the shares that we are "missing"
        for (row_idx, col_idx) in missing {
            let share = unsafe { eds1.square.0.get_unchecked_mut(row_idx, col_idx) };
            *share = [0u8; 64]
        }

        eds1.repair(&missing).unwrap();

        assert_eq!(eds0.square.0, eds1.square.0)
    }

    #[test]
    fn extension_repair_preflight_incorrect_root() {
        let arr: [[u8; 64]; 4] = (0..4)
            .map(|i| [i + 1; 64])
            .collect::<Vec<_>>()
            .try_into()
            .unwrap();
        let ds = DataSquare::new(&arr).unwrap();
        let mut eds = ds.extend().unwrap();
        eds.row_roots[0][0] = !eds.row_roots[0][0];
        let err = eds.repair(&[]).unwrap_err();
        assert!(err.is::<ByzantineDataError<64>>())
    }

    #[test]
    fn extension_repair_incorrect_root() {
        let arr: [[u8; 64]; 4] = (0..4)
            .map(|i| [i + 1; 64])
            .collect::<Vec<_>>()
            .try_into()
            .unwrap();
        let ds = DataSquare::new(&arr).unwrap();
        let mut eds = ds.extend().unwrap();
        eds.row_roots[0][0] = !eds.row_roots[0][0];
        let err = eds.repair(&[(0, 0)]).unwrap_err();
        assert!(err.is::<ByzantineDataError<64>>())
    }

    #[test]
    fn extension_repair_incorrect_ortho_root() {
        let arr: [[u8; 64]; 4] = (0..4)
            .map(|i| [i + 1; 64])
            .collect::<Vec<_>>()
            .try_into()
            .unwrap();
        let ds = DataSquare::new(&arr).unwrap();
        let mut eds = ds.extend().unwrap();
        eds.col_roots[0][0] = !eds.col_roots[0][0];
        let err = eds.repair(&[(0, 0)]).unwrap_err();
        assert!(err.is::<ByzantineDataError<64>>())
    }
}
