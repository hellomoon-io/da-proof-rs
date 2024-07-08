use num::pow::Pow;

use crate::gf16::poly::{Poly, Term};
use crate::gf16::GF16;

pub fn new_poly(symbol_count: u16) -> Poly {
    (0..symbol_count).fold(Poly::from(Term::new(1.into(), 0)), |acc, i| {
        let multiplier = Term::new(1.into(), 1) + Term::new(GF16::from(2).pow(i.into()), 0);
        acc * multiplier
    })
}

#[cfg(test)]
mod test {
    #[test]
    fn new_poly_4() {
        let poly = super::new_poly(4);
        assert_eq!(
            poly.coeffs(),
            &vec![
                0x01.into(),
                0x0F.into(),
                0x36.into(),
                0x78.into(),
                0x40.into()
            ]
        )
    }
}
