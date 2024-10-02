use ark_ff::{Field, MontBackend};
use ark_std::{mem::swap, ops::Mul};

use crate::{
    cyclotomic_ring::{CyclotomicConfig, CyclotomicPolyRingGeneral, CyclotomicPolyRingNTTGeneral},
    poly_ring::PolyRing,
    Cyclotomic, OverField,
};

mod decomposition;
mod ntt;

mod fq_def {
    #![allow(non_local_definitions)]
    use ark_ff::{Fp256, MontBackend};

    #[derive(ark_ff::MontConfig)]
    #[modulus = "3618502788666131213697322783095070105623107215331596699973092056135872020481"]
    #[generator = "3"]
    pub struct FqConfig;
    pub type Fq = Fp256<MontBackend<FqConfig, 4>>;
}

pub use fq_def::*;

pub type RqNTT = CyclotomicPolyRingNTTGeneral<StarkRingConfig, 4, { ntt::N }>;
pub type RqPoly = CyclotomicPolyRingGeneral<StarkRingConfig, 4, { ntt::D }>;

pub struct StarkRingConfig;

impl CyclotomicConfig<4> for StarkRingConfig {
    type BaseFieldConfig = MontBackend<FqConfig, 4>;

    type BaseCRTField = Fq;

    const CRT_FIELD_EXTENSION_DEGREE: usize = 1;

    fn reduce_in_place(coefficients: &mut Vec<ark_ff::Fp<Self::BaseFieldConfig, 4>>) {
        let (left, right) = coefficients.split_at_mut(ntt::D);
        for (coeff_left, coeff_right) in left.iter_mut().zip(right) {
            *coeff_left -= coeff_right;
        }

        coefficients.resize(ntt::D, <Fq as Field>::ZERO);
    }

    #[inline(always)]
    fn crt_in_place(coefficients: &mut Vec<Fq>) {
        ntt::stark_prime_crt_in_place(coefficients);
    }

    #[inline(always)]
    fn crt(coefficients: Vec<Fq>) -> Vec<Fq> {
        ntt::stark_prime_crt(coefficients)
    }

    #[inline(always)]
    fn icrt(evaluations: Vec<Fq>) -> Vec<Fq> {
        ntt::stark_prime_icrt(evaluations)
    }

    #[inline(always)]
    fn icrt_in_place(evaluations: &mut Vec<Fq>) {
        ntt::stark_prime_icrt_in_place(evaluations);
    }
}

impl OverField for RqPoly {}
impl OverField for RqNTT {}

impl Mul<Fq> for RqNTT {
    type Output = Self;

    fn mul(self, rhs: Fq) -> Self {
        Self(self.0.map(|x| x * rhs))
    }
}

impl From<Fq> for RqNTT {
    fn from(value: Fq) -> Self {
        Self::from_scalar(value)
    }
}

impl Cyclotomic for RqPoly {
    fn rot(&mut self) {
        let mut buf = -self.0[Self::dimension() - 1];

        for i in 0..Self::dimension() {
            swap(&mut buf, &mut self.0[i]);
        }
    }
}

#[cfg(test)]
mod test {
    use ark_poly::{
        univariate::{DenseOrSparsePolynomial, DensePolynomial, SparsePolynomial},
        DenseUVPolynomial,
    };
    use ark_std::UniformRand;
    use rand::thread_rng;

    use super::*;
    use crate::{balanced_decomposition::Decompose, PolyRing, Ring};

    #[test]
    fn test_implements_decompose() {
        fn takes_decompose<T: Decompose>(_x: T) {}

        let x = RqPoly::ONE;

        takes_decompose(x);
    }

    #[test]
    fn test_crt_one() {
        let one = RqPoly::ONE;

        assert_eq!(RqNTT::from(one), RqNTT::ONE)
    }

    #[test]
    fn test_icrt_one() {
        let one = RqNTT::ONE;

        assert_eq!(RqPoly::from(one), RqPoly::ONE)
    }

    #[test]
    fn test_reduce() {
        let mut rng = thread_rng();
        let mut coeffs: Vec<Fq> = (0..(2 * ntt::D)).map(|_| Fq::rand(&mut rng)).collect();

        let poly = DensePolynomial::from_coefficients_slice(&coeffs);

        // X^16 + 1
        let the_cyclotomic = SparsePolynomial::from_coefficients_slice(&[
            (16, <Fq as Field>::ONE),
            (0, <Fq as Field>::ONE),
        ]);

        let (_, r) =
            DenseOrSparsePolynomial::divide_with_q_and_r(&poly.into(), &the_cyclotomic.into())
                .unwrap();

        StarkRingConfig::reduce_in_place(&mut coeffs);

        assert_eq!(r.coeffs(), &coeffs);
    }

    #[test]
    fn test_mul_crt() {
        let mut rng = thread_rng();

        let coeff_1 = RqPoly::rand(&mut rng);
        let coeff_2 = RqPoly::rand(&mut rng);

        let ntt_form_1 = RqNTT::from(coeff_1);
        let ntt_form_2 = RqNTT::from(coeff_2);

        let ntt_mul = ntt_form_1 * ntt_form_2;
        let coeffs_mul = coeff_1 * coeff_2;

        // ntt_mul.coeffs() performs INTT while coeffs_mul.coeffs() just returns the coefficients
        assert_eq!(RqPoly::from(ntt_mul).coeffs(), coeffs_mul.coeffs());
    }

    #[test]
    fn test_cyclotomic() {
        let mut rng = thread_rng();

        let mut a = RqPoly::rand(&mut rng);
        let initial_a = a;

        let x = RqPoly::x();

        for i in 1..RqPoly::dimension() {
            a.rot();
            assert_eq!(a, initial_a * x.pow([i as u64]));
        }
    }
}
