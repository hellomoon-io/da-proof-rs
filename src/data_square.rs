use crate::bisetmap::BisetMap;
use crate::common::MerkleProof;
use crate::util::{as_chunks, as_chunks_mut, as_flat, invariant, sqrt_ceil};
use borsh::{BorshDeserialize, BorshSerialize};
use generic_array::{ArrayLength, GenericArray};
use grid::Grid;
use itertools::Itertools;
use num::integer::sqrt;
use reed_solomon_erasure::galois_16::ReedSolomon;
use rs_merkle::algorithms::Sha256;
use rs_merkle::{Hasher, MerkleTree};

use crate::common::{new_share, Axis, MerkleRoot, NonZeroMultipleOf64, Share, Square};
use crate::error::{
    ByzantineData, DataNotSquare, ShareSizeDoesNotMatch, TooManyErasures,
    UnableToConstructMerkleRoot,
};

/// A "square" of data aligned into `Share`s.
///
/// *Invariant*: the data is square, i.e. the grid is equally as wide as it is high.
#[derive(Debug, Clone, PartialEq)]
pub struct DataSquare<Sz: NonZeroMultipleOf64>(Square<Sz>);

impl<Sz: NonZeroMultipleOf64> BorshSerialize for DataSquare<Sz> {
    fn serialize<W: std::io::Write>(&self, writer: &mut W) -> std::io::Result<()> {
        // Serialize the size of the shares
        Sz::U64.serialize(writer)?;
        // Serialize the width of the square
        (self.0.rows() as u64).serialize(writer)?;
        self.0.iter_rows().try_for_each(|mut row| {
            row.try_for_each(|elt| elt.iter().try_for_each(|byte| byte.serialize(writer)))
        })
    }
}

impl<Sz: NonZeroMultipleOf64> BorshDeserialize for DataSquare<Sz> {
    fn deserialize_reader<R: std::io::Read>(reader: &mut R) -> std::io::Result<Self> {
        let sz = u64::deserialize_reader(reader)?;
        invariant!(
            sz == Sz::U64,
            std::io::Error::new(
                std::io::ErrorKind::Other,
                ShareSizeDoesNotMatch {
                    actual: sz,
                    expected: Sz::U64
                }
            )
        );
        let width = u64::deserialize_reader(reader)?;
        let mut grid = Grid::new(width as usize, width as usize);
        (0..width).try_for_each(|row| {
            (0..width).try_for_each(|col| {
                let mut share = GenericArray::<u8, Sz>::default();
                (0..Sz::USIZE).try_for_each(|i| {
                    share[i] = u8::deserialize_reader(reader)?;
                    Ok::<(), std::io::Error>(())
                })?;

                grid[(row as usize, col as usize)] = share;

                Ok::<(), std::io::Error>(())
            })
        })?;
        Ok(Self(grid))
    }
}

/// An extended "square" of data.
///
/// This includes the erasure encoded data derived from the original `DataSquare`.
#[derive(Debug, Clone, PartialEq)]
pub struct ExtendedDataSquare<Sz: NonZeroMultipleOf64> {
    original_width: usize,
    square: DataSquare<Sz>,
    row_roots: Vec<MerkleRoot>,
    col_roots: Vec<MerkleRoot>,
}

impl<Sz: NonZeroMultipleOf64> BorshSerialize for ExtendedDataSquare<Sz> {
    fn serialize<W: std::io::Write>(&self, writer: &mut W) -> std::io::Result<()> {
        (self.original_width as u64).serialize(writer)?;
        self.square.serialize(writer)?;
        self.row_roots.serialize(writer)?;
        self.col_roots.serialize(writer)
    }
}

impl<Sz: NonZeroMultipleOf64> BorshDeserialize for ExtendedDataSquare<Sz> {
    fn deserialize_reader<R: std::io::Read>(reader: &mut R) -> std::io::Result<Self> {
        let original_width = u64::deserialize_reader(reader)?;
        let square = DataSquare::<Sz>::deserialize_reader(reader)?;
        let row_roots = Vec::<MerkleRoot>::deserialize_reader(reader)?;
        let col_roots = Vec::<MerkleRoot>::deserialize_reader(reader)?;
        Ok(Self {
            original_width: original_width as usize,
            square,
            row_roots,
            col_roots,
        })
    }
}

/// A two way map of missing shares in an EDS
type MissingMap = BisetMap<usize, usize>;

impl<Sz: NonZeroMultipleOf64> ByzantineData<Sz> {
    /// Constructs a new `ByzantineDataError` with `shares` populated.
    fn new_from_eds(eds: &ExtendedDataSquare<Sz>, axis: Axis, idx: usize) -> Self {
        Self {
            axis,
            idx,
            shares: (match axis {
                Axis::Row => eds.square.0.iter_row(idx),
                Axis::Col => eds.square.0.iter_col(idx),
            })
            .map(ToOwned::to_owned)
            .collect_vec(),
        }
    }
}

impl TooManyErasures {
    fn new_from_missing_map(missing_map: &MissingMap) -> Self {
        Self {
            missing: missing_map.flat_iter().collect_vec(),
        }
    }
}

/// Computes the root of the Merkle tree with the leaves as shares in the square.
fn compute_root<Sz: ArrayLength>(
    square: &Square<Sz>,
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
        .ok_or(UnableToConstructMerkleRoot::Axis(axis, idx).into())
}

fn compute_proof<Sz: ArrayLength>(
    square: &Square<Sz>,
    axis: Axis,
    row: usize,
    col: usize,
) -> MerkleProof {
    let (iter, ortho_idx) = match axis {
        Axis::Row => (square.iter_row(row), col),
        Axis::Col => (square.iter_col(col), row),
    };
    let leaves: Vec<MerkleRoot> = iter.map(|share| Sha256::hash(share)).collect();
    let tree = MerkleTree::<Sha256>::from_leaves(&leaves);
    tree.proof(&[ortho_idx]).into()
}

impl<Sz: NonZeroMultipleOf64> DataSquare<Sz> {
    /// Constructs a new `DataSquare` from a flat slice containing `Share`s.
    ///
    /// Pads the square if the number of shares is not a perfect square.
    pub fn new_pad(data: &[Share<Sz>]) -> Self {
        let width = sqrt_ceil(data.len());

        let ga = GenericArray::<u8, Sz>::default();

        let mut inner: Square<Sz> = Grid::init(width, width, ga);

        (0..width).for_each(|row_idx| {
            inner
                .iter_row_mut(row_idx)
                .enumerate()
                .for_each(|(i, elt)| {
                    *elt = match data.get(row_idx * width + i) {
                        Some(share) => share.clone(),
                        None => new_share(),
                    }
                });
        });

        Self(inner)
    }

    /// Constructs a `DataSquare` from a flat slice containing `Share`s.
    ///
    /// Returns `Err(_)` if the length of the input is not a perfect square.
    pub fn new(data: &[Share<Sz>]) -> anyhow::Result<Self> {
        let width = sqrt(data.len());

        invariant!(
            width * width == data.len(),
            DataNotSquare { len: data.len() }.into()
        );

        Ok(Self::new_pad(data))
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
                let mut row_view = self
                    .0
                    .iter_row_mut(i)
                    .map(|share| as_chunks_mut::<u8, 2>(share))
                    .collect::<anyhow::Result<Vec<_>>>()?;
                rs.encode(&mut row_view)?;
            };
            {
                let mut col_view = self
                    .0
                    .iter_col_mut(i)
                    .map(|share| as_chunks_mut::<u8, 2>(share))
                    .collect::<anyhow::Result<Vec<_>>>()?;
                rs.encode(&mut col_view)?;
            };
            Ok::<(), anyhow::Error>(())
        })?;

        // Encode the bottom right quadrant of the array using the rows from the column
        // encoding of the original data.
        //
        // NOTE: this should be equivalent to encoding the columns of the row encoding
        // of the original data (this is checked when not in release mode).
        (width..width * 2).try_for_each(|i| {
            let mut row_view = self
                .0
                .iter_row_mut(i)
                .map(|share| as_chunks_mut::<u8, 2>(share))
                .collect::<anyhow::Result<Vec<_>>>()?;
            rs.encode(&mut row_view)?;
            Ok::<(), anyhow::Error>(())
        })?;

        // Invariant: encoding Q3 from the rows of Q2 should be equivalent to encoding it from the
        // columns of Q1.
        #[cfg(debug_assertions)]
        (width..width * 2).try_for_each(|i| {
            let col_view = self
                .0
                .iter_col(i)
                .map(|share| as_chunks::<u8, 2>(share))
                .collect::<anyhow::Result<Vec<_>>>()?;

            let mut parity = vec![new_share::<Sz>(); width];
            let mut parity_view = parity
                .iter_mut()
                .map(|x| as_chunks_mut::<u8, 2>(x))
                .collect::<anyhow::Result<Vec<_>>>()?;

            rs.encode_sep(&col_view[..width], &mut parity_view[..])?;

            assert_eq!(&col_view[width..], &parity_view[..]);
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
    pub fn extend(mut self) -> anyhow::Result<ExtendedDataSquare<Sz>> {
        let original_width = self.0.cols();
        (0..self.0.cols()).for_each(|_| self.0.push_col(vec![new_share(); self.0.rows()]));
        (0..self.0.rows()).for_each(|_| self.0.push_row(vec![new_share(); self.0.cols()]));
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
    fn join(&self, other: &Self) -> Self {
        Self {
            solved: self.solved && other.solved,
            progress_made: self.progress_made || other.progress_made,
        }
    }
}

impl<Sz: NonZeroMultipleOf64> ExtendedDataSquare<Sz> {
    /// TODO: verify encoding

    /// Before executing a repair, check if the slice is valid.
    ///
    /// No-op if the slice is incomplete (according to the user).
    fn repair_preflight_slice<'a>(
        &self,
        axis: Axis,
        idx: usize,
        missing: &MissingMap,
    ) -> anyhow::Result<()> {
        let has_missing = match axis {
            Axis::Row => missing.exists_left(&idx),
            Axis::Col => missing.exists_right(&idx),
        };
        if !has_missing {
            let roots = match axis {
                Axis::Row => &self.row_roots,
                Axis::Col => &self.col_roots,
            };
            let expected_root = roots
                .get(idx)
                .ok_or_else(|| ByzantineData::new_from_eds(self, axis, idx))?;
            let calculated_root = compute_root(&self.square.0, axis, idx)?;
            if expected_root != &calculated_root {
                Err(ByzantineData::new_from_eds(self, axis, idx).into())
            } else {
                Ok(())
            }
        } else {
            Ok(())
        }
    }

    /// Before executing a repair, check if `self` contains any known Byzantine data.
    fn repair_preflight(&self, missing: &MissingMap) -> anyhow::Result<()> {
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
        missing: &mut MissingMap,
    ) -> anyhow::Result<RepairStep> {
        let rs = ReedSolomon::new(self.original_width, self.original_width)?;

        let (missing_indices_opt, iter, expected_root) = match axis {
            Axis::Row => (
                missing.get_left(&idx).cloned(),
                self.square.0.iter_row(idx),
                self.row_roots[idx],
            ),
            Axis::Col => (
                missing.get_right(&idx).cloned(),
                self.square.0.iter_col(idx),
                self.col_roots[idx],
            ),
        };

        if let Some(missing_indices) = missing_indices_opt {
            // Copy the shares into a new array where `shares[i] = Some(_)` if the element is
            // present and `None` if it is not.
            let mut shares = iter
                .enumerate()
                .map(|(jdx, share)| {
                    if missing_indices.contains(&jdx) {
                        Ok(None)
                    } else {
                        Ok(Some(
                            as_chunks::<u8, 2>(share)?.iter().cloned().collect_vec(),
                        ))
                    }
                })
                .collect::<anyhow::Result<Vec<_>>>()?;

            // We "ignore" failures here because we want to keep trying to solve.
            if let Err(_) = rs.reconstruct(&mut shares) {
                return Ok(RepairStep {
                    solved: false,
                    progress_made: false,
                });
            };

            let decompressed = shares
                .iter()
                .map(|reconstructed_opt| {
                    let reconstructed = reconstructed_opt.as_ref().unwrap();
                    let mut share: Share<Sz> = new_share();
                    share.clone_from_slice(as_flat(&reconstructed)?);
                    Ok(share)
                })
                .collect::<anyhow::Result<Vec<Share<Sz>>>>()?;

            let calculated_root = {
                let leaves = decompressed
                    .iter()
                    .map(|share| Sha256::hash(share))
                    .collect_vec();
                let tree = MerkleTree::<Sha256>::from_leaves(&leaves);
                tree.root()
                    .ok_or_else(|| ByzantineData::new_from_eds(self, axis, idx))
            }?;

            if calculated_root != expected_root {
                return Err(ByzantineData::new_from_eds(self, axis, idx).into());
            }

            missing_indices.iter().try_for_each(|jdx| {
                let (row_idx, col_idx) = match axis {
                    Axis::Row => (idx, *jdx),
                    Axis::Col => (*jdx, idx),
                };

                let share = unsafe { self.square.0.get_unchecked_mut(row_idx, col_idx) };
                *share = decompressed[*jdx].clone();

                Ok::<(), anyhow::Error>(())
            })?;

            missing_indices.iter().try_for_each(|jdx| {
                // Get the missing status for the orthogonal slice of the index we just updated.
                //
                // If there are no more missing shares in the slice, check the slice against the
                // expected root and the RS encoding.
                let (row_idx, col_idx, ortho_missing_opt, ortho_iter, ortho_expected_root) =
                    match axis {
                        Axis::Row => (
                            idx,
                            *jdx,
                            missing.get_right(jdx),
                            self.square.0.iter_col(*jdx),
                            self.col_roots[*jdx],
                        ),
                        Axis::Col => (
                            *jdx,
                            idx,
                            missing.get_left(jdx),
                            self.square.0.iter_row(*jdx),
                            self.row_roots[*jdx],
                        ),
                    };

                // The orthogonal slice is complete if the only missing share is the current one
                //
                // We already upated it in the square, but we haven't removed it from the map of
                // missing indices.
                let ortho_complete = ortho_missing_opt.map_or(false, |ortho_missing| {
                    ortho_missing.len() == 1 && ortho_missing.contains(jdx)
                });

                if ortho_complete {
                    let ortho_calculated_root = {
                        // This clone just clones the references, not the underlying data.
                        let leaves = ortho_iter
                            .clone()
                            .map(|share| Sha256::hash(share))
                            .collect_vec();
                        let tree = MerkleTree::<Sha256>::from_leaves(&leaves);
                        tree.root()
                            .ok_or_else(|| ByzantineData::new_from_eds(self, axis.ortho(), *jdx))
                    }?;

                    // Invariant: the computed Merkle root of the slice should match the "known" one.
                    invariant!(
                        ortho_expected_root == ortho_calculated_root,
                        ByzantineData::new_from_eds(self, axis.ortho(), *jdx).into()
                    );

                    let shares_view = ortho_iter
                        .clone()
                        .map(|x| as_chunks(x))
                        .collect::<anyhow::Result<Vec<&[[u8; 2]]>>>()?;

                    let mut parity = vec![new_share::<Sz>(); self.original_width];

                    let mut parity_view = parity
                        .iter_mut()
                        .map(|x| as_chunks_mut(x))
                        .collect::<anyhow::Result<Vec<&mut [[u8; 2]]>>>()?;
                    rs.encode_sep(&shares_view[..self.original_width], &mut parity_view)?;

                    // Invariant: the parity data of the shares should match the recomputed parity.
                    invariant!(
                        &shares_view[..self.original_width] == &parity_view[..],
                        ByzantineData::new_from_eds(self, axis.ortho(), *jdx).into()
                    )
                }

                missing.remove(&row_idx, &col_idx);

                Ok::<(), anyhow::Error>(())
            })?;

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
    fn solve(&mut self, missing: &mut MissingMap) -> anyhow::Result<()> {
        loop {
            let repair_step: RepairStep =
                (0..self.original_width * 2).try_fold(RepairStep::default(), |acc, i| {
                    let row_repair_step = self.solve_slice(Axis::Row, i, missing)?;
                    let col_repair_step = self.solve_slice(Axis::Col, i, missing)?;

                    Ok::<_, anyhow::Error>(acc.join(&row_repair_step).join(&col_repair_step))
                })?;

            if repair_step.solved {
                return Ok(());
            } else if !repair_step.progress_made {
                return Err(TooManyErasures::new_from_missing_map(missing).into());
            }
        }
    }

    /// Repair the EDS, assuming the indices in `missing` are erased/corrupted.
    pub fn repair(&mut self, missing: &[(usize, usize)]) -> anyhow::Result<()> {
        // Assemble missing maps
        let mut missing_map = BisetMap::new();
        missing.iter().for_each(|(row, col)| {
            missing_map.insert(*row, *col);
        });
        self.repair_preflight(&missing_map)?;
        self.solve(&mut missing_map)?;
        Ok(())
    }

    pub fn data_root(&self) -> anyhow::Result<MerkleRoot> {
        let leaves = self
            .row_roots
            .iter()
            .chain(self.col_roots.iter())
            .map(ToOwned::to_owned)
            .collect_vec();

        let tree = MerkleTree::<Sha256>::from_leaves(&leaves[..]);
        tree.root().ok_or(UnableToConstructMerkleRoot::Data.into())
    }

    pub fn merkle_proof_share(&self, axis: Axis, row: usize, col: usize) -> MerkleProof {
        compute_proof(&self.square.0, axis, row, col)
    }

    pub fn merkle_proof_axis(&self, axis: Axis, idx: usize) -> MerkleProof {
        let leaves = self
            .row_roots
            .iter()
            .chain(self.col_roots.iter())
            .map(ToOwned::to_owned)
            .collect_vec();
        let tree = MerkleTree::<Sha256>::from_leaves(&leaves[..]);
        let idx = match axis {
            Axis::Row => idx,
            Axis::Col => self.original_width * 2 + 1,
        };
        tree.proof(&[idx]).into()
    }

    pub fn verify_merkle_proof_share(
        &self,
        axis: Axis,
        row: usize,
        col: usize,
        merkle_proof: &MerkleProof,
    ) -> bool {
        let (root, idx) = match axis {
            Axis::Row => (self.row_roots[row], row),
            Axis::Col => (self.col_roots[col], col),
        };

        merkle_proof.0.verify(
            root,
            &[idx],
            &[Sha256::hash(&self.square.0[(row, col)])],
            self.original_width * 2,
        )
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use borsh::{BorshDeserialize, BorshSerialize};
    use generic_array::arr;
    use typenum::{U128, U256};

    fn make_ds<const L: usize>() -> DataSquare<U256> {
        let mut i = 0;
        let arr: [Share<_>; L] = [(); L].map(|()| {
            i = i + 1;
            arr![i; 256]
        });
        DataSquare::new_pad(&arr)
    }

    #[test]
    fn extension() {
        let ds = make_ds::<4>();
        let _eds = ds.extend(); // The debug assert should fail if this test doesn't pass
    }

    #[test]
    fn new_pad_pads() {
        let ds = make_ds::<3>();
        let empty_share = new_share::<U256>();
        assert_eq!(ds.0[(1, 1)], empty_share)
    }

    #[test]
    fn extension_repair_max() {
        let ds = make_ds::<4>();
        let eds0 = ds.extend().unwrap();
        let mut eds1 = eds0.clone();
        // -----------------
        // | x |   | x | x |
        // |---|---|---|---|
        // | x | x | x | x |
        // |---|---|---|---|
        // | x | x | x |   |
        // |---|---|---|---|
        // | x |   |   | x |
        // -----------------
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
            *share = arr![0u8; 256]
        }

        eds1.repair(&missing).unwrap();

        assert_eq!(eds0.square.0, eds1.square.0)
    }

    #[test]
    fn extension_repair_too_many_erasures() {
        let ds = make_ds::<4>();
        let eds0 = ds.extend().unwrap();
        let mut eds1 = eds0.clone();
        // -----------------
        // | x |   | x | x |
        // |---|---|---|---|
        // | x | x | x | x |
        // |---|---|---|---|
        // | x | x | x |   |
        // |---|---|---|---|
        // | x | x |   | x |
        // -----------------
        // (no row has >=2 known shares)
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
            (3, 1),
            (3, 3),
        ];

        // Zero out the shares that we are "missing"
        for (row_idx, col_idx) in missing {
            let share = unsafe { eds1.square.0.get_unchecked_mut(row_idx, col_idx) };
            *share = new_share()
        }

        let err = eds1.repair(&missing).unwrap_err();
        assert!(err.is::<TooManyErasures>())
    }

    #[test]
    fn extension_repair_preflight_incorrect_root() {
        let ds = make_ds::<4>();
        let mut eds = ds.extend().unwrap();
        eds.row_roots[0][0] = !eds.row_roots[0][0];
        let err = eds.repair(&[]).unwrap_err();
        assert!(err.is::<ByzantineData<U256>>())
    }

    #[test]
    fn extension_repair_incorrect_root() {
        let ds = make_ds::<4>();
        let mut eds = ds.extend().unwrap();
        eds.row_roots[0][0] = !eds.row_roots[0][0];
        let err = eds.repair(&[(0, 0)]).unwrap_err();
        assert!(err.is::<ByzantineData<U256>>())
    }

    #[test]
    fn extension_repair_incorrect_ortho_root() {
        let ds = make_ds::<4>();
        let mut eds = ds.extend().unwrap();
        eds.col_roots[0][0] = !eds.col_roots[0][0];
        let err = eds.repair(&[(0, 0)]).unwrap_err();
        assert!(err.is::<ByzantineData<U256>>())
    }

    #[test]
    fn extension_repair_incorrect_ortho_encoding() {
        let ds = make_ds::<4>();
        let mut eds = ds.extend().unwrap();
        eds.square.0[(1, 0)] = new_share();

        eds.col_roots[0] = compute_root(&eds.square.0, Axis::Col, 0).unwrap();

        let err = eds.repair(&[(0, 0)]).unwrap_err();
        assert!(err.is::<ByzantineData<U256>>())
    }

    #[test]
    fn ds_serialization_deserialization() {
        let ds = make_ds::<4>();
        let mut buf = Vec::<u8>::new();
        ds.serialize(&mut buf).unwrap();
        let desered = DataSquare::<U256>::deserialize(&mut &buf[..]).unwrap();
        assert_eq!(ds, desered);
    }

    #[test]
    fn ds_deserialize_incorrect_size() {
        let ds = make_ds::<4>();
        let mut buf = Vec::<u8>::new();
        ds.serialize(&mut buf).unwrap();
        let err = DataSquare::<U128>::deserialize(&mut &buf[..]).unwrap_err();
        let err = err
            .get_ref()
            .unwrap()
            .downcast_ref::<ShareSizeDoesNotMatch>()
            .unwrap();

        assert_eq!(err.expected, 128);
        assert_eq!(err.actual, 256)
    }

    #[test]
    fn eds_serialization_deserialization() {
        let ds = make_ds::<4>();
        let eds = ds.extend().unwrap();
        let mut buf = Vec::<u8>::new();
        eds.serialize(&mut buf).unwrap();
        let desered = ExtendedDataSquare::<U256>::deserialize(&mut &buf[..]).unwrap();
        assert_eq!(eds, desered);
    }

    #[test]
    fn merkle_proof_verification_ok() {
        let ds = make_ds::<4>();
        let eds = ds.extend().unwrap();
        let proof = eds.merkle_proof_share(Axis::Row, 1, 1);
        assert!(eds.verify_merkle_proof_share(Axis::Row, 1, 1, &proof));
    }

    #[test]
    fn merkle_proof_verification_err() {
        let ds = make_ds::<4>();
        let eds = ds.extend().unwrap();
        let proof = eds.merkle_proof_share(Axis::Row, 1, 1);
        assert!(!eds.verify_merkle_proof_share(Axis::Row, 1, 2, &proof));
    }
}
