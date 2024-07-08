use std::env::VarError;

use anyhow::ensure;
use num::pow::Pow;
use num::traits::Inv;
use num::{One, ToPrimitive, Zero};

use crate::common::TWO_TO_SIXTEEEN;
use crate::gf16::poly::{Poly, Term};
use crate::gf16::GF16;

mod generator;

pub struct ReedSolomon {
    symbol_count: u16,
    polynomial: Poly,
}

#[derive(Debug)]
pub struct Message(Vec<u16>);

trait AsPoly {
    fn as_poly<'a>(&'a self) -> &'a Poly {
        unsafe { &*(self as *const Self as *const Poly) }
    }
}

impl AsPoly for Message {}

#[derive(Debug, Clone)]
pub struct Syndromes(Vec<u16>);

impl Syndromes {
    pub fn check(&self) -> bool {
        match self.0.iter().max() {
            Some(0) | None => true,
            _ => false,
        }
    }
}

impl AsPoly for Syndromes {}

impl ReedSolomon {
    pub fn new(symbol_count: u16) -> Self {
        Self {
            symbol_count,
            polynomial: generator::new_poly(symbol_count),
        }
    }

    pub fn symbol_count(&self) -> u16 {
        self.symbol_count
    }

    pub fn encode(&self, msg: &[u16]) -> anyhow::Result<Message> {
        ensure!(
            msg.len() + (self.symbol_count() as usize) < TWO_TO_SIXTEEEN,
            "Message length ({}) + symbol count ({}) is too large.",
            msg.len(),
            self.symbol_count()
        );

        let msg = GF16::slice_from_u16_slice(msg);

        let mut temp = vec![GF16::zero(); msg.len() + self.polynomial.degree()];
        temp[..msg.len()].clone_from_slice(msg);

        // TODO: get rid of spurious clone.
        //
        // Not a huge deal since the vector is relatively small but still performance-sensitive.
        let (_, rem) = Poly::from(temp) / self.polynomial.clone();

        let vec: Vec<GF16> = msg
            .iter()
            .chain(Into::<Vec<GF16>>::into(rem).iter())
            .map(ToOwned::to_owned)
            .collect();

        Ok(Message(GF16::u16_vec_from_vec(vec)))
    }

    fn calculate_syndromes(&self, msg: &Message) -> Syndromes {
        let mut res = vec![GF16::zero(); self.symbol_count() as usize + 1];
        res.iter_mut().skip(1).enumerate().for_each(|(i, v)| {
            let i = i.to_i64().unwrap();
            let y = msg.as_poly().eval(GF16(2).pow(i));
            *v = y
        });

        Syndromes(GF16::u16_vec_from_vec(res))
    }

    fn compute_errata_locator(&self, coeff_positions: &[usize]) -> Poly {
        let poly_one = Poly::from(Term::new(1.into(), 0));
        coeff_positions.iter().fold(poly_one.clone(), |poly, pos| {
            let i = pos.to_i64().unwrap();
            let summand = Term::new(GF16(2).pow(i), 1) + Term::new(0.into(), 0);
            let multiplicand = poly_one.clone() + summand;

            poly * multiplicand
        })
    }

    fn compute_error_evaluator(&self, syndromes: &Syndromes, errata_locator: &Poly) -> Poly {
        let mut syndromes = syndromes.clone();
        syndromes.0.reverse();

        let dividend = syndromes.as_poly() * errata_locator;
        let divisor = Poly::from(Term::new(1.into(), (self.symbol_count() + 1) as usize));
        // TODO: make this quicker
        (dividend / divisor).1
    }

    fn correct_errata(&self, msg: &Message, errata_positions: &[usize]) -> anyhow::Result<Message> {
        let coeff_positions: Vec<usize> = errata_positions
            .iter()
            .map(|pos| msg.0.len() - pos - 1)
            .collect();
        let errata_locator = self.compute_errata_locator(&coeff_positions[..]);

        let syndromes = self.calculate_syndromes(msg);
        let mut error_evaluator = self.compute_error_evaluator(&syndromes, &errata_locator);
        error_evaluator.0.reverse();

        let x = (0..coeff_positions.len()).map(|pos| {
            let l = TWO_TO_SIXTEEEN as i64 - 1 - pos as i64;
            GF16(2).pow(-l)
        });

        let mut e = vec![GF16::zero(); msg.0.len()];

        x.clone().enumerate().try_for_each(|(i, xi)| {
            let xi_inv = xi.inv();

            let err_loc_prime = x
                .clone()
                .enumerate()
                .filter_map(|(j, xj)| {
                    if i != j {
                        Some(GF16(1) - (xi_inv * xj))
                    } else {
                        None
                    }
                })
                .fold(GF16::one(), |acc, coeff| acc * coeff);
            let y = xi.pow(1) * error_evaluator.eval(xi_inv);

            ensure!(
                err_loc_prime != GF16::zero(),
                "Unable to determine error location."
            );

            let magnitude = y / err_loc_prime;

            e[errata_positions[i]] = magnitude;
            Ok(())
        })?;

        Ok(Message(GF16::u16_vec_from_vec(
            (msg.as_poly() + &Poly(e)).0,
        )))
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn uncorrupted_syndromes_are_zero() {
        let rs = ReedSolomon::new(8);
        let encoded = rs.encode(&[0xdead, 0xbeef]).unwrap();
        let syndromes = rs.calculate_syndromes(&encoded);
        assert!(syndromes.check())
    }

    #[test]
    fn corrupted_syndromes_are_nonzero() {
        let rs = ReedSolomon::new(8);
        let mut encoded = rs.encode(&[1234, 4321]).unwrap();
        encoded.0[0] = 0x1234;
        let syndromes = rs.calculate_syndromes(&encoded);
        assert!(!syndromes.check())
    }

    #[test]
    fn correct_errata() {
        let rs = ReedSolomon::new(8);
        let mut encoded = rs.encode(&[1234, 4321]).unwrap();
        encoded.0[0] = 0;
        println!("CORRUPTED: {:?}", encoded);
        let corrected = rs.correct_errata(&encoded, &[0]);
        todo!("CORRECTED: {:?}", corrected);
    }
}
