//! The degree-four extension used at the Hachi / Greyhound boundary.
//!
//! This basis is adapted to the subfield of the cyclotomic ring. It is not
//! the byte representation of the extension in the original Hachi prototype.

use ark_ff::{Fp2, Fp2Config, Fp4, Fp4Config};

pub const Q: u64 = (1u64 << 32) - 99;
pub(crate) const fn fq(value: u64) -> Fq {
    lattirust_arithmetic::ring::f_p::FqConfig::<Q>::from_u64(value)
}

pub type Fq = lattirust_arithmetic::ring::f_p::Fq<Q>;

pub struct Fq2Config;
impl Fp2Config for Fq2Config {
    type Fp = Fq;
    const NONRESIDUE: Fq = fq(2);
    const FROBENIUS_COEFF_FP2_C1: &'static [Fq] = &[fq(1), fq(4294967196)];
}

pub type Fq2 = Fp2<Fq2Config>;

pub struct Fq4Config;
impl Fp4Config for Fq4Config {
    type Fp2Config = Fq2Config;
    // Fq4 = Fq2[T] / (T² - √2).
    const NONRESIDUE: Fq2 = Fq2::new(fq(0), fq(1));
    const FROBENIUS_COEFF_FP4_C1: &'static [Fq] =
        &[fq(1), fq(983270775), fq(4294967196), fq(3311696422)];
}

pub type Fq4 = Fp4<Fq4Config>;

#[cfg(test)]
mod tests {
    use super::*;
    use ark_ff::Field;

    #[test]
    fn quartic_frobenius_matches_exponentiation() {
        let a = Fq4::new(
            Fq2::new(7u64.into(), 3u64.into()),
            Fq2::new(2u64.into(), 11u64.into()),
        );
        let mut by_powers = a;
        for i in 0..=4 {
            let mut by_frobenius = a;
            by_frobenius.frobenius_map_in_place(i);
            assert_eq!(by_frobenius, by_powers);
            by_powers = by_powers.pow([Q]);
        }
    }
}
