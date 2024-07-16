use generic_array::ArrayLength;

use crate::common::{Axis, NonZeroMultipleOf64, Share};

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
    /// The shares that could be used to construct a fraud proof
    pub shares: Vec<Option<Share<Sz>>>,
}

impl<Sz: NonZeroMultipleOf64> ByzantineData<Sz> {}

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

/// An error that is returned when the user deserializes a share with an incorrect size.
#[derive(Debug)]
pub struct ShareSizeDoesNotMatch {
    pub actual: u64,
    pub expected: u64,
}

impl std::fmt::Display for ShareSizeDoesNotMatch {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "expected a share size of {}, got {}",
            self.expected, self.actual
        )
    }
}

impl std::error::Error for ShareSizeDoesNotMatch {}

/// An error that is returned when unable to compute a Merkle root
#[derive(Debug)]
pub enum UnableToConstructMerkleRoot {
    Axis(Axis, usize),
    Data,
}

impl std::fmt::Display for UnableToConstructMerkleRoot {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "unable to construct merkle root ")?;
        match self {
            UnableToConstructMerkleRoot::Axis(axis, col) => write!(f, "for {axis}, {col}"),
            UnableToConstructMerkleRoot::Data => write!(f, "for data"),
        }
    }
}

impl std::error::Error for UnableToConstructMerkleRoot {}
