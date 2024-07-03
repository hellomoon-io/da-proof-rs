//! This module contains types for doing finite field arithmetic in GF(2^16).
//!
//! Much of the implementation is based off of the paper "Fast Software Implementations of Finite Field Operations" by C. Huang and L. Xu.
//! (https://lihaoxu.eng.wayne.edu/NISL/Papers/Tech/GF.pdf)
//!
//! A 16-bit Galois field is used in the data availability proof.

use crate::gf16_lut::LUT;

use crate::util::BitLength;

const TWO_TO_SIXTEEEN: usize = 1 << 16;

/// An irreducible polynomial of degree 16.
///
/// Equivalent to x^16 + x^5 + x^3 + x^2 + 1
///
/// Source: https://www.wolframalpha.com/input?i=primitive+polynomials+of+GF%2865536%29
pub const IRREDUCIBLE_POLYNOMIAL_DEGREE_16: u32 = 0b10000000000101101;

#[derive(Clone, Copy, Eq, PartialEq, Hash)]
pub struct GF16(pub(crate) u16);

/// Formats a GF16 number `x` as `GF16(x)`.
impl std::fmt::Debug for GF16 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "GF16(")?;
        <u16 as std::fmt::Debug>::fmt(&self.0, f)?;
        write!(f, ")")?;
        Ok(())
    }
}

/// Formats a GF16 number `x` as `x`.
impl std::fmt::Display for GF16 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        <u16 as std::fmt::Display>::fmt(&self.0, f)
    }
}

/// Formats a GF16 number `x` as hex.
impl std::fmt::LowerHex for GF16 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        <u16 as std::fmt::LowerHex>::fmt(&self.0, f)
    }
}

impl From<u16> for GF16 {
    fn from(value: u16) -> Self {
        Self(value)
    }
}

impl From<GF16> for u16 {
    fn from(value: GF16) -> Self {
        value.0
    }
}

/// Computes long multiplication where addition is equivalent to XOR.
///
/// This function is slow and should not be generally used -- it's mostly meant for building the
/// lookup table for fast arithmetic.
fn mul_xor(lhs: u16, rhs: u16) -> u32 {
    // The result needs to be wider than the operands in the case of overflow.
    let mut res: u32 = 0;
    let mut i = 0;

    // Rust panics if the RHS of a shift will result in 0.
    while rhs.checked_shr(i).unwrap_or(0) > 0 {
        if (rhs & 1 << i) != 0 {
            res ^= (lhs as u32) << i;
        }
        i += 1;
    }
    res
}

/// Performs long division where subtraction is equivalent to XOR.
///
/// This function is slow and should not be generally used -- it's mostly meant for building the
/// lookup table for fast arithmetic.
fn div_xor(mut lhs: u32, rhs: u32) -> u32 {
    let lhs_bl = lhs.bit_length();
    let rhs_bl = rhs.bit_length();

    if lhs_bl >= rhs_bl {
        for i in (0..(lhs_bl - rhs_bl + 1) as isize).rev() {
            if lhs & (1 << i + (rhs_bl as isize) - 1) != 0 {
                lhs ^= (rhs as u32) << i;
            }
        }
    }

    lhs
}

impl GF16 {
    /// Performs slow multiplication within GF(2^16).
    ///
    /// This function is slow and should not be generally used -- it's mostly meant for building the
    /// lookup table for fast arithmetic.
    fn mul_slow(self, rhs: Self) -> Self {
        Self(
            div_xor(mul_xor(self.0, rhs.0), IRREDUCIBLE_POLYNOMIAL_DEGREE_16)
                .try_into() // This is safe because we've reduced the result to be within GF(2^16)
                .unwrap(),
        )
    }
}

impl num::traits::Zero for GF16 {
    fn zero() -> Self {
        0.into()
    }

    fn is_zero(&self) -> bool {
        self.0 == 0
    }
}

#[derive(Debug)]
pub struct GF16LUT {
    pub(crate) exp: [GF16; TWO_TO_SIXTEEEN],
    pub(crate) log: [GF16; TWO_TO_SIXTEEEN],
}

impl GF16LUT {
    pub fn generate() -> Self {
        let mut exp = [GF16(0); TWO_TO_SIXTEEEN];
        let mut log = [GF16(0); TWO_TO_SIXTEEEN];

        let _: GF16 = (0..TWO_TO_SIXTEEEN).fold(GF16::from(1), |acc, i| {
            exp[i] = acc;
            log[acc.0 as usize] = GF16(i.try_into().unwrap());
            acc.mul_slow(2.into())
        });

        Self { exp, log }
    }
}

/// Addition in GF16 is equivalent to XORing, since there is no carrying each bit.
///
/// This is essentially equivalent to polynomial addition mod 2 for each bit.
impl std::ops::Add for GF16 {
    type Output = GF16;

    fn add(self, rhs: Self) -> Self::Output {
        GF16(self.0 ^ rhs.0)
    }
}

/// See `Add` implementation.
impl std::ops::AddAssign for GF16 {
    fn add_assign(&mut self, rhs: Self) {
        self.0 ^= rhs.0
    }
}

/// Subtraction in GF16 is equivalent to XORing, since there is no carrying each bit.
///
/// This is essentially equivalent to polynomial subtraction mod 2 for each bit.
impl std::ops::Sub for GF16 {
    type Output = GF16;

    fn sub(self, rhs: Self) -> Self::Output {
        GF16(self.0 ^ rhs.0)
    }
}

/// See `Sub` implementation.
impl std::ops::SubAssign for GF16 {
    fn sub_assign(&mut self, rhs: Self) {
        self.0 ^= rhs.0
    }
}

/// Constant time multiplication using pre-generated lookup tables.
impl std::ops::Mul for GF16 {
    type Output = GF16;

    fn mul(self, rhs: Self) -> Self::Output {
        if self.0 == 0 || rhs.0 == 0 {
            0.into()
        } else {
            LUT.exp[(LUT.log[self.0 as usize].0 as usize + LUT.log[rhs.0 as usize].0 as usize)
                % (TWO_TO_SIXTEEEN - 1)]
        }
    }
}

impl std::ops::MulAssign for GF16 {
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}

/// Constant time division using pre-generated lookup tables.
impl std::ops::Div for GF16 {
    type Output = GF16;

    fn div(self, rhs: Self) -> Self::Output {
        if rhs.0 == 0 {
            panic!("attempt to divide GF16 by zero")
        }

        if self.0 == 0 {
            0.into()
        } else {
            LUT.exp[(LUT.log[self.0 as usize].0 as usize + (TWO_TO_SIXTEEEN - 1)
                - LUT.log[rhs.0 as usize].0 as usize)
                % (TWO_TO_SIXTEEEN - 1)]
        }
    }
}

/// Constant time exponentiation using pre-generated lookup tables.
impl num::traits::Pow<GF16> for GF16 {
    type Output = GF16;

    fn pow(self, rhs: GF16) -> Self::Output {
        LUT.exp[(LUT.log[self.0 as usize].0 as usize * rhs.0 as usize) % (TWO_TO_SIXTEEEN - 1)]
    }
}

/// Constant time multiplicative inverse using pre-generated lookup tables.
impl num::traits::Inv for GF16 {
    type Output = GF16;

    fn inv(self) -> Self::Output {
        LUT.exp[(TWO_TO_SIXTEEEN - 1) - LUT.log[self.0 as usize].0 as usize]
    }
}

mod poly {
    use num::Zero;
    use std::cmp::max;

    use super::GF16;

    #[derive(Debug, Clone, PartialEq)]
    pub struct Poly(Vec<GF16>);

    #[derive(Debug, Clone)]
    pub struct Term {
        coef: GF16,
        degree: usize,
    }

    impl Term {
        pub fn new(coef: GF16, degree: usize) -> Self {
            Self { coef, degree }
        }
    }

    impl std::ops::Add for Term {
        type Output = Poly;

        fn add(self, rhs: Term) -> Self::Output {
            let mut vec_self = vec![GF16::zero(); self.degree + 1];
            let mut vec_rhs = vec![GF16::zero(); rhs.degree + 1];
            vec_self[0] = self.coef;
            vec_rhs[0] = rhs.coef;
            Poly(vec_self) + Poly(vec_rhs)
        }
    }

    impl Poly {
        pub fn degree(&self) -> usize {
            self.0.len() - 1
        }

        /// The coefficients in the polynomial in descending degree order.
        ///
        /// i.e. 2x^2 + 1 would be represented as `[2, 0, 1]`.
        pub fn coeffs(&self) -> &Vec<GF16> {
            &self.0
        }

        /// Scales the polynomial by `scalar`.
        pub fn scale(&self, scalar: GF16) -> Self {
            Poly(self.0.iter().map(|v| *v * scalar).collect())
        }

        /// Scales the polynomial in place by `scalar`.
        pub fn scale_mut(&mut self, scalar: GF16) {
            self.0.iter_mut().for_each(|v| *v *= scalar)
        }

        /// Evaluates the polynomial at `x` using Horner's method.
        pub fn eval(&self, x: GF16) -> GF16 {
            let mut iter = self.0.iter();
            match iter.next() {
                Some(v) => iter.fold(*v, |acc, term| (acc * x) + *term),
                None => 0.into(),
            }
        }
    }

    /// Polynomial addition using the standard "whiteboard" method.
    ///
    /// The result stays within GF(2^16).
    impl std::ops::Add for Poly {
        type Output = Self;

        fn add(self, rhs: Self) -> Self::Output {
            let sz = max(self.0.len(), rhs.0.len());
            let mut res = vec![GF16::zero(); sz];
            for i in 0..self.0.len() {
                res[i + sz - self.0.len()] = self.0[i];
            }
            for i in 0..rhs.0.len() {
                res[i + sz - rhs.0.len()] += rhs.0[i];
            }
            Poly(res)
        }
    }

    /// Adds a term to the polynomial.
    impl std::ops::Add<Term> for Poly {
        type Output = Self;

        fn add(self, rhs: Term) -> Self::Output {
            // First we construct a temporary polynomial including just `rhs`, and add that to `self`.
            let mut rhs_vec = vec![GF16::zero(); rhs.degree + 1];
            rhs_vec[0] = rhs.coef;
            self + Poly(rhs_vec)
        }
    }

    /// Polynomial multiplication using the standard "whiteboard" method.
    ///
    /// The result stays within GF(2^16).
    impl std::ops::Mul for Poly {
        type Output = Self;

        fn mul(self, rhs: Self) -> Self::Output {
            let sz = self.0.len() + rhs.0.len() - 1;
            let mut res = vec![GF16::zero(); sz];
            for j in 0..rhs.0.len() {
                for i in 0..self.0.len() {
                    res[i + j] += self.0[i] * rhs.0[j];
                }
            }
            Poly(res)
        }
    }

    /// Synthetic division for polynomials, optimized for GF(2^16).
    impl std::ops::Div for Poly {
        type Output = (Self, Self);

        fn div(mut self, rhs: Self) -> Self::Output {
            let divisor_degree = rhs.0.len() - 1;
            if self.0.len() < divisor_degree {
                (Self(vec![]), self)
            } else {
                for i in 0..(self.0.len() - divisor_degree) {
                    let coef = self.0[i];
                    if coef != GF16::zero() {
                        for j in 1..rhs.0.len() {
                            if rhs.0[j] != GF16::zero() {
                                self.0[i + j] += rhs.0[j] * coef;
                            }
                        }
                    }
                }

                let sep = self.0.len() - (rhs.0.len() - 1);
                (Self(self.0[..sep].into()), Self(self.0[sep..].into()))
            }
        }
    }

    /// Negation is the identity function in GF(2^16), since all operations are carried out modulo 2.
    impl std::ops::Neg for Poly {
        type Output = Self;

        fn neg(self) -> Self::Output {
            self
        }
    }
}

#[cfg(test)]
mod test {
    use std::collections::HashSet;
    use std::hash::RandomState;

    use num::pow::Pow;
    use num::traits::Inv;
    use num::Zero;

    use super::poly::*;
    use super::*;

    #[test]
    fn lut_entries_unique() {
        let lut = GF16LUT::generate();
        assert_eq!(
            HashSet::<GF16, RandomState>::from_iter(lut.exp.clone().into_iter()).len(),
            super::TWO_TO_SIXTEEEN - 1
        );
        assert_eq!(lut.exp[0], lut.exp[super::TWO_TO_SIXTEEEN - 1]);
        assert_eq!(
            HashSet::<GF16, RandomState>::from_iter(lut.log.clone().into_iter()).len(),
            super::TWO_TO_SIXTEEEN
        );
    }

    #[quickcheck_macros::quickcheck]
    fn mul_equivalent_to_mul_slow(a: u16, b: u16) -> bool {
        let a = GF16(a);
        let b = GF16(b);
        a * b == a.mul_slow(b)
    }

    #[quickcheck_macros::quickcheck]
    fn mul_div_isomorphism(a: u16, b: u16) -> bool {
        let a = GF16(a);
        let b = GF16(b);
        b == 0.into() || (a * b) / b == a
    }

    #[quickcheck_macros::quickcheck]
    fn mul_inv_one(a: u16) -> bool {
        let a = GF16(a);
        a == 0.into() || a * a.inv() == 1.into()
    }

    #[quickcheck_macros::quickcheck]
    #[ignore] // This is pretty slow since it does a lot of slow multiplication
    fn pow_isomorphism(a: u16, b_: u16) {
        let a = GF16(a);
        let b = GF16(b_);
        let powed = a.pow(b);
        if a == 0.into() || b == 0.into() {
            assert_eq!(powed, 1.into(), "{} ** {} == {} != 1", a, b, powed);
        } else {
            let muled = (1..(b_ as usize) + 1).fold(GF16(1), |acc, _| acc * a);
            assert_eq!(a.pow(b), muled, "{} ** {} == {} != {}", a, b, powed, muled);
        }
    }

    #[test]
    fn term_add() {
        let poly = Term::new(GF16(0x01), 4) + Term::new(GF16(0x0f), 3);
        assert_eq!(
            poly.coeffs(),
            &vec![
                0x01.into(),
                0x0f.into(),
                GF16::zero(),
                GF16::zero(),
                GF16::zero()
            ]
        )
    }

    #[test]
    fn poly_add_term() {
        let poly = Term::new(GF16(0x01), 4)
            + Term::new(GF16(0x0f), 3)
            + Term::new(GF16(0x36), 2)
            + Term::new(GF16(0x78), 1)
            + Term::new(GF16(0x40), 0);
        assert_eq!(
            poly.coeffs(),
            &vec![
                0x01.into(),
                0x0f.into(),
                0x36.into(),
                0x78.into(),
                0x40.into()
            ]
        )
    }

    #[test]
    fn poly_div() {
        let dividend =
            Term::new(0x12.into(), 6) + Term::new(0x34.into(), 5) + Term::new(0x56.into(), 4);
        let divisor = Term::new(GF16(0x01), 4)
            + Term::new(0x0f.into(), 3)
            + Term::new(0x36.into(), 2)
            + Term::new(0x78.into(), 1);

        let (q, r) = dividend.clone() / divisor.clone();

        // The following results are from the galois Python library.
        //
        // See the Jupyter notebook in the root of the project.
        assert_eq!(
            q,
            Term::new(0x12.into(), 2) + Term::new(0xDA.into(), 1) + Term::new(0x78C.into(), 0)
        );
        assert_eq!(
            r,
            Term::new(0x3988.into(), 3) + Term::new(0xBED8.into(), 2) + Term::new(0x560D.into(), 1)
        );

        assert_eq!(q * divisor + r, dividend);
    }
}
