use std::ops::Rem;

use generic_array::{ArrayLength, GenericArray};
use grid::Grid;
use typenum::{IsNotEqual, U0, U64};

pub const TWO_TO_SIXTEEEN: usize = 1 << 16;

pub trait NonZeroMultipleOf64: ArrayLength + std::fmt::Debug {}
impl<T: ArrayLength + Rem<U64, Output = U0> + IsNotEqual<U0> + std::fmt::Debug> NonZeroMultipleOf64
    for T
{
}

/// Type aliases

pub type Share<Sz> = GenericArray<u8, Sz>;
pub type MerkleRoot = [u8; 32];
pub(crate) type Square<Sz> = Grid<Share<Sz>>;

#[derive(Debug, Clone, Copy)]
pub enum Axis {
    Row,
    Col,
}

/// Creates a new empty share of size `Sz`.
pub fn new_share<Sz: ArrayLength + std::fmt::Debug>() -> Share<Sz> {
    GenericArray::default()
}

impl std::fmt::Display for Axis {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Axis::Row => write!(f, "row"),
            Axis::Col => write!(f, "col"),
        }
    }
}

impl Axis {
    pub fn ortho(&self) -> Self {
        match self {
            Axis::Row => Axis::Col,
            Axis::Col => Axis::Row,
        }
    }
}
