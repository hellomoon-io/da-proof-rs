use core::slice;
use std::mem::size_of;

use generic_array::GenericArray;
use num::integer::{sqrt, Roots};
use num::PrimInt;

pub trait BitLength: Sized {
    const BITS: usize = size_of::<Self>() * 8;

    fn bit_length(&self) -> usize;
}

impl<T: PrimInt> BitLength for T {
    fn bit_length(&self) -> usize {
        Self::BITS - self.leading_zeros() as usize
    }
}

/// Early returns for invariant conditions
macro_rules! invariant {
    ($cond:expr, $err:expr) => {
        if !$cond {
            return Err($err);
        }
    };
}

pub(crate) use invariant;

use crate::common::{NonZeroMultipleOf64, Share};
use crate::error::InvalidChunkSize;

/// Returns a chunked view of `slice`.
///
/// Returns `Err(InvalidChunkSize(_))` if `N == 0` or `N % slice.len() != 0`.
pub fn as_chunks<T, const N: usize>(slice: &[T]) -> anyhow::Result<&[[T; N]]> {
    invariant!(
        N != 0 && slice.len() % N == 0,
        InvalidChunkSize {
            len: slice.len(),
            chunk_size: N
        }
        .into()
    );

    // This is safe because the two conditions we care about are satisfied:
    // - `N != 0` - this would result in a division-by-zero error
    // - `len % N == 0` - this would result in left over chunks
    Ok(unsafe { slice::from_raw_parts(slice.as_ptr().cast(), slice.len() / N) })
}

/// Returns a mutable chunked view of `slice`.
///
/// Returns `Err(InvalidChunkSize(_))` if `N == 0` or `N % slice.len() != 0`.
pub fn as_chunks_mut<T, const N: usize>(slice: &mut [T]) -> anyhow::Result<&mut [[T; N]]> {
    invariant!(
        N != 0 && slice.len() % N == 0,
        InvalidChunkSize {
            len: slice.len(),
            chunk_size: N
        }
        .into()
    );

    // See above safety statement
    Ok(unsafe { slice::from_raw_parts_mut(slice.as_mut_ptr().cast(), slice.len() / N) })
}

/// Returns a flattened view of a chunked slice.
///
/// Returns `Err(InvalidChunkSize(_))` if `N == 0`.
pub fn as_flat<T, const N: usize>(slice: &[[T; N]]) -> anyhow::Result<&[T]> {
    invariant!(
        N != 0,
        InvalidChunkSize {
            len: slice.len(),
            chunk_size: N
        }
        .into()
    );

    // This is safe because `N != 0`
    Ok(unsafe { slice::from_raw_parts(slice.as_ptr().cast(), slice.len() * N) })
}

/// Returns a mutable flattened view of a chunked slice.
///
/// Returns `Err(InvalidChunkSize(_))` if `N == 0`.
pub fn as_flat_mut<T, const N: usize>(slice: &mut [[T; N]]) -> anyhow::Result<&mut [T]> {
    invariant!(
        N != 0,
        InvalidChunkSize {
            len: slice.len(),
            chunk_size: N
        }
        .into()
    );

    Ok(unsafe { slice::from_raw_parts_mut(slice.as_mut_ptr().cast(), slice.len() * N) })
}

/// Compute the square root of `v`, rounding up if `v` is not a perfect square.
pub fn sqrt_ceil<T: Roots + std::marker::Copy>(v: T) -> T {
    let res = sqrt(v);
    if res * res == v {
        res
    } else {
        res + T::one()
    }
}

/// Returns a view of `slice` as shares.
///
/// Returns `Err(InvalidChunkSize(_))` if `slice.len() % Sz == 0`.
pub fn as_shares<Sz: NonZeroMultipleOf64>(slice: &[u8]) -> anyhow::Result<&[Share<Sz>]> {
    invariant!(
        slice.len() % Sz::USIZE == 0,
        InvalidChunkSize {
            len: slice.len(),
            chunk_size: Sz::USIZE
        }
        .into()
    );

    Ok(GenericArray::<u8, Sz>::chunks_from_slice(slice).0)
}

/// Returns a view of `slice` as shares.
///
/// The second argument of the tuple is the "remainder".
pub fn as_shares_rem<Sz: NonZeroMultipleOf64>(slice: &[u8]) -> (&[Share<Sz>], &[u8]) {
    GenericArray::<u8, Sz>::chunks_from_slice(slice)
}

#[cfg(test)]
mod test {
    use typenum::U64;

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

    #[quickcheck_macros::quickcheck]
    fn sqrt_ceil_u32(v: u32) -> bool {
        let r = sqrt_ceil(v);
        v <= 1 || (r as u64) * (r as u64) >= (v as u64) && (r - 1) * (r - 1) < v
    }

    #[test]
    fn as_chunks_ok() {
        let arr = [0u8, 1u8, 2u8, 3u8];
        let chunks = as_chunks::<u8, 2>(&arr).unwrap();
        assert_eq!(chunks, &[[0u8, 1u8], [2u8, 3u8]])
    }

    #[test]
    fn as_chunks_mutation_persists() {
        let mut arr = [0u8, 1u8, 2u8, 3u8];
        let chunks = as_chunks_mut::<u8, 2>(&mut arr).unwrap();
        chunks[0][0] = 0xFF;
        assert_eq!(arr[0], 0xFF);
    }

    #[test]
    fn as_chunks_err() {
        let arr = [0u8, 1u8, 2u8, 3u8];
        let err = as_chunks::<u8, 3>(&arr).unwrap_err();
        assert!(err.is::<InvalidChunkSize>())
    }

    #[test]
    fn as_flat_ok() {
        let arr = [[0u8, 1u8], [2u8, 3u8]];
        let flat = as_flat(&arr).unwrap();
        assert_eq!(flat, &[0u8, 1u8, 2u8, 3u8])
    }

    #[test]
    fn as_flat_mutation_persists() {
        let mut arr = [[0u8, 1u8], [2u8, 3u8]];
        let flat = as_flat_mut(&mut arr).unwrap();
        flat[0] = 0xFF;
        assert_eq!(arr[0][0], 0xFF)
    }

    #[test]
    fn as_shares_ok() {
        let mut i = 0u8;
        let arr = [(); 128].map(|()| {
            let v = i / 64;
            i += 1;
            v
        });
        let shares = as_shares::<U64>(&arr).unwrap();
        shares[0].iter().for_each(|v| assert_eq!(v, &0u8));
        shares[1].iter().for_each(|v| assert_eq!(v, &1u8));
    }
}
