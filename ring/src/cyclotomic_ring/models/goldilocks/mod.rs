use ark_ff::{Field, Fp3, Fp3Config, MontBackend, MontFp};

use crate::{
    cyclotomic_ring::{CyclotomicConfig, CyclotomicPolyRingGeneral, CyclotomicPolyRingNTTGeneral},
    traits::FromRandomBytes,
    Ring,
};

mod ntt;

mod fq_def {
    #![allow(non_local_definitions)]
    use ark_ff::{Fp64, MontBackend, MontConfig};
    #[derive(MontConfig)]
    #[modulus = "18446744069414584321"]
    #[generator = "7"]
    pub struct FqConfig;
    pub type Fq = Fp64<MontBackend<FqConfig, 1>>;
}

pub use fq_def::*;

pub struct Goldilocks3Config;

pub type RqNTT = CyclotomicPolyRingNTTGeneral<GoldilocksRingConfig, 1, { ntt::N }>;
pub type RqPoly = CyclotomicPolyRingGeneral<GoldilocksRingConfig, 1, { ntt::D }>;

impl Fp3Config for Goldilocks3Config {
    type Fp = Fq;

    const NONRESIDUE: Self::Fp = MontFp!("1099511627776");

    const TWO_ADICITY: u32 = 32;

    const TRACE_MINUS_ONE_DIV_TWO: &'static [u64] =
        &[9223372049739677694, 9223372049739677692, 2147483646];

    // 7 ^ t.
    const QUADRATIC_NONRESIDUE_TO_T: ark_ff::Fp3<Self> = Fp3::new(
        MontFp!("3607031617444012685"),
        <Fq as Field>::ZERO,
        <Fq as Field>::ZERO,
    );

    // Do we even need these?
    // NQR ^ (MODULUS^i - 1)/3, i=0,1,2 with NQR = u = (0,1,0)
    const FROBENIUS_COEFF_FP3_C1: &'static [Fq] = &[];
    // NQR ^ (2*MODULUS^i - 2)/3, i=0,1,2 with NQR = u = (0,1,0)
    const FROBENIUS_COEFF_FP3_C2: &'static [Fq] = &[];
}

pub type Fq3 = Fp3<Goldilocks3Config>;

impl Ring for Fq3 {
    const ZERO: Self = <Fq3 as Field>::ZERO;
    const ONE: Self = <Fq3 as Field>::ONE;
}

impl FromRandomBytes<Fq3> for Fq3 {
    #[inline(always)]
    fn byte_size() -> usize {
        3 * 8
    }

    fn try_from_random_bytes(bytes: &[u8]) -> Option<Fq3> {
        Self::from_random_bytes(bytes)
    }
}

pub struct GoldilocksRingConfig;

impl CyclotomicConfig<1> for GoldilocksRingConfig {
    type BaseFieldConfig = MontBackend<FqConfig, 1>;
    type BaseCRTField = Fq3;

    const CRT_FIELD_EXTENSION_DEGREE: usize = 3;

    fn reduce_in_place(coefficients: &mut Vec<Fq>) {
        for i in 0..ntt::D / 2 {
            let a24_i = coefficients
                .get(ntt::D + i)
                .copied()
                .unwrap_or(<Fq as Field>::ZERO);
            let a24_12_i = coefficients
                .get(ntt::D + ntt::D / 2 + i)
                .copied()
                .unwrap_or(<Fq as Field>::ZERO);
            coefficients[i] -= a24_i;
            coefficients[i] -= a24_12_i;
        }

        for i in ntt::D / 2..ntt::D {
            let a_12_i = coefficients
                .get(ntt::D / 2 + i)
                .copied()
                .unwrap_or(<Fq as Field>::ZERO);
            coefficients[i] += a_12_i;
        }

        coefficients.resize(ntt::D, <Fq as Field>::ZERO);
    }

    #[inline(always)]
    fn crt_in_place(coefficients: &mut Vec<Fq>) {
        ntt::goldilocks_crt_in_place(coefficients);
    }

    #[inline(always)]
    fn crt(coefficients: Vec<Fq>) -> Vec<Fq3> {
        ntt::goldilocks_crt(coefficients)
    }

    #[inline(always)]
    fn icrt(evaluations: Vec<Fq3>) -> Vec<Fq> {
        ntt::goldilocks_icrt(evaluations)
    }

    #[inline(always)]
    fn icrt_in_place(evaluations: &mut Vec<Fq>) {
        ntt::goldilocks_icrt_in_place(evaluations);
    }
}

mod test {
    #![allow(unused_imports)]
    use crate::PolyRing;

    use super::*;
    use ark_ff::One;
    use ark_poly::{
        univariate::{DenseOrSparsePolynomial, DensePolynomial, SparsePolynomial},
        DenseMVPolynomial, DenseUVPolynomial,
    };
    use ark_std::UniformRand;
    use rand::thread_rng;

    // Note: if ord X = 24 then X can't be a cubic residue.
    #[test]
    fn test_nonresidue_is_order_24() {
        let nonresidue: Fq = Goldilocks3Config::NONRESIDUE;

        let mut pow = nonresidue;

        for _i in 0..22 {
            pow *= nonresidue;
            assert_ne!(pow, <Fq as Field>::ONE);
        }

        assert_eq!(pow * nonresidue, <Fq as Field>::ONE);
    }

    #[test]
    fn test_reduce() {
        let mut rng = thread_rng();
        let mut coeffs: Vec<Fq> = (0..(2 * ntt::D)).map(|_| Fq::rand(&mut rng)).collect();

        let poly = DensePolynomial::from_coefficients_slice(&coeffs);

        // X^24 - X^12 + 1
        let the_cyclotomic = SparsePolynomial::from_coefficients_slice(&[
            (24, <Fq as Field>::ONE),
            (12, -<Fq as Field>::ONE),
            (0, <Fq as Field>::ONE),
        ]);

        let (_, r) =
            DenseOrSparsePolynomial::divide_with_q_and_r(&poly.into(), &the_cyclotomic.into())
                .unwrap();

        GoldilocksRingConfig::reduce_in_place(&mut coeffs);

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
}
