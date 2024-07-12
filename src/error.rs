use generic_array::ArrayLength;

use crate::common::{Axis, Share};

/// An error that is returned when Byzantine (invalid, malicious, etc.) data is detected.
///
/// This can be a result of the shares in an axis not matching the known Merkle roots,
/// or recovery data that doesn't match the original.
#[derive(Debug)]
pub struct ByzantineData<Sz: ArrayLength> {
    /// The axis on which Byzantine data was detected
    pub axis: Axis,
    /// The index of the axis on which Byzantine data was detected
    pub idx: usize,
    /// The shares that could be used to construct a fraud proof.
    pub shares: Vec<Share<Sz>>,
}

impl<Sz: ArrayLength> std::fmt::Display for ByzantineData<Sz> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "byzantine data detected @ {} {}", self.axis, self.idx)
    }
}

impl<Sz: ArrayLength + std::fmt::Debug> std::error::Error for ByzantineData<Sz> {}

/// An error that is returned when there are too many erasures to recover the data.
#[derive(Debug)]
pub struct TooManyErasures {
    /// The indices of the data that is still missing in the array.
    pub missing: Vec<(usize, usize)>,
}

impl std::fmt::Display for TooManyErasures {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "too many erasures, missing {} shares",
            self.missing.len()
        )
    }
}

impl std::error::Error for TooManyErasures {}

/// An error that is returned when the user attempts to construct a square with a non-square slice.
#[derive(Debug)]
pub struct DataNotSquare {
    pub len: usize,
}

impl std::fmt::Display for DataNotSquare {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "data length must be perfect square, not {}", self.len)
    }
}

impl std::error::Error for DataNotSquare {}

/// An error that is returned when the user attempts to flatten/chunk a slice on incongruent bounds.
#[derive(Debug)]
pub struct InvalidChunkSize {
    pub len: usize,
    pub chunk_size: usize,
}

impl std::fmt::Display for InvalidChunkSize {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "unable to evenly chunk slice of length {} into {}",
            self.len, self.chunk_size
        )
    }
}

impl std::error::Error for InvalidChunkSize {}
