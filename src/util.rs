use std::mem::size_of;

use num::PrimInt;

#[derive(Debug, Clone, Copy)]
pub enum Axis {
    Row,
    Col,
}

impl std::fmt::Display for Axis {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Axis::Row => write!(f, "Row"),
            Axis::Col => write!(f, "Col"),
        }
    }
}

pub trait BitLength: Sized {
    const BITS: usize = size_of::<Self>() * 8;

    fn bit_length(&self) -> usize;
}

impl<T: PrimInt> BitLength for T {
    fn bit_length(&self) -> usize {
        Self::BITS - self.leading_zeros() as usize
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[quickcheck_macros::quickcheck]
    fn u32_bit_length(v: u32) -> bool {
        let mut v_ = v;
        let mut sz = 0;

        while v_ > 0 {
            sz += 1;
            v_ >>= 1;
        }

        sz == v.bit_length()
    }
}
