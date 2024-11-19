use ark_ff::{Field, Fp3, Fp3Config, MontBackend, MontFp};
use ark_std::{mem::swap, ops::Mul, vec::*};

pub mod fq9;
pub mod ntt;

use crate::{
    cyclotomic_ring::{
        CyclotomicConfig, CyclotomicPolyRingGeneral, CyclotomicPolyRingNTTGeneral, CRT, ICRT,
    },
    impl_crt_icrt_for_a_ring,
    traits::FromRandomBytes,
    Cyclotomic, OverField, PolyRing, Ring,
};

use self::fq9::{Fp9, Fp9Config};

mod fq_def {
    #![allow(non_local_definitions)]
    use ark_ff::{Fp64, MontBackend, MontConfig};
    #[derive(MontConfig)]
    #[modulus = "2013265921"]
    #[generator = "31"]
    pub struct FqConfig;
    pub type Fq = Fp64<MontBackend<FqConfig, 1>>;
}

pub use fq_def::*;

pub type RqNTT = CyclotomicPolyRingNTTGeneral<BabyBearRingConfig, 1, { ntt::N }>;
pub type RqPoly = CyclotomicPolyRingGeneral<BabyBearRingConfig, 1, { ntt::D }>;

pub struct BabyBear3ExtConfig;

impl Fp3Config for BabyBear3ExtConfig {
    type Fp = Fq;

    const NONRESIDUE: Self::Fp = MontFp!("503591070");

    const TWO_ADICITY: u32 = 27;

    // Parameter below are placeholders and not needed
    const TRACE_MINUS_ONE_DIV_TWO: &'static [u64] = &[];
    const QUADRATIC_NONRESIDUE_TO_T: Fp3<Self> = Fp3::new(
        <Fq as Field>::ZERO,
        <Fq as Field>::ZERO,
        <Fq as Field>::ZERO,
    );
    const FROBENIUS_COEFF_FP3_C1: &'static [Self::Fp] = &[];
    const FROBENIUS_COEFF_FP3_C2: &'static [Self::Fp] = &[];
}

pub struct BabyBear9ExtConfig;

impl Fp9Config for BabyBear9ExtConfig {
    type Fp3Config = BabyBear3ExtConfig;

    const NONRESIDUE: Fp3<Self::Fp3Config> = Fp3::new(
        <<BabyBear3ExtConfig as Fp3Config>::Fp as Field>::ZERO,
        <<BabyBear3ExtConfig as Fp3Config>::Fp as Field>::ONE,
        <<BabyBear3ExtConfig as Fp3Config>::Fp as Field>::ZERO,
    );

    const SQRT_PRECOMP: Option<ark_ff::SqrtPrecomputation<fq9::Fp9<Self>>> = None;

    const FROBENIUS_COEFF_FP9_C1: &'static [<Self::Fp3Config as Fp3Config>::Fp] = &[];
    const FROBENIUS_COEFF_FP9_C2: &'static [<Self::Fp3Config as Fp3Config>::Fp] = &[];
}

pub type Fq9 = Fp9<BabyBear9ExtConfig>;

impl Ring for Fq9 {
    const ZERO: Self = <Fq9 as Field>::ZERO;
    const ONE: Self = <Fq9 as Field>::ONE;
}

impl FromRandomBytes<Fq9> for Fq9 {
    #[inline(always)]
    fn byte_size() -> usize {
        9 * 8
    }

    fn try_from_random_bytes(bytes: &[u8]) -> Option<Fq9> {
        Self::from_random_bytes(bytes)
    }
}

pub struct BabyBearRingConfig;

impl CyclotomicConfig<1> for BabyBearRingConfig {
    type BaseFieldConfig = MontBackend<FqConfig, 1>;
    type BaseCRTField = Fq9;

    const CRT_FIELD_EXTENSION_DEGREE: usize = 9;

    fn reduce_in_place(coefficients: &mut Vec<Fq>) {
        for i in 0..ntt::D / 2 {
            let a72_i = coefficients
                .get(ntt::D + i)
                .copied()
                .unwrap_or(<Fq as Field>::ZERO);
            let a72_36_i = coefficients
                .get(ntt::D + ntt::D / 2 + i)
                .copied()
                .unwrap_or(<Fq as Field>::ZERO);
            coefficients[i] -= a72_i;
            coefficients[i] -= a72_36_i;
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
    fn crt_in_place(coefficients: &mut [Fq]) {
        ntt::babybear_crt_in_place(coefficients);
    }

    #[inline(always)]
    fn crt(coefficients: Vec<Fq>) -> Vec<Fq9> {
        ntt::babybear_crt(coefficients)
    }

    #[inline(always)]
    fn icrt(evaluations: Vec<Fq9>) -> Vec<Fq> {
        ntt::babybear_icrt(evaluations)
    }

    #[inline(always)]
    fn icrt_in_place(evaluations: &mut [Fq]) {
        ntt::babybear_icrt_in_place(evaluations)
    }
}

impl OverField for RqPoly {}
impl OverField for RqNTT {}

impl Mul<Fq9> for RqNTT {
    type Output = Self;

    fn mul(self, rhs: Fq9) -> Self {
        Self(self.0.map(|x| x * rhs))
    }
}

impl From<Fq9> for RqNTT {
    fn from(value: Fq9) -> Self {
        Self::from_scalar(value)
    }
}

impl Cyclotomic for RqPoly {
    fn rot(&mut self) {
        let last = self.0[Self::dimension() - 1];
        let mut buf = -last;

        for i in 0..Self::dimension() {
            swap(&mut buf, &mut self.0[i]);
        }

        self.0[36] += last;
    }
}

impl_crt_icrt_for_a_ring!(RqNTT, RqPoly, BabyBearRingConfig);

#[cfg(test)]
mod tests {
    use ark_poly::{
        univariate::{DenseOrSparsePolynomial, DensePolynomial, SparsePolynomial},
        DenseUVPolynomial,
    };
    use ark_std::UniformRand;
    use rand::SeedableRng;
    use rand_chacha::ChaCha8Rng;

    use super::*;
    use crate::{
        balanced_decomposition::Decompose,
        cyclotomic_ring::crt::{CRT, ICRT},
        PolyRing,
    };

    #[test]
    fn test_implements_decompose() {
        fn takes_decompose<T: Decompose>(_x: T) {}

        let x = RqPoly::ONE;

        takes_decompose(x);
    }

    #[test]
    fn test_crt_one() {
        let one = RqPoly::ONE;

        assert_eq!(one.crt(), RqNTT::ONE)
    }

    #[test]
    fn test_icrt_one() {
        let one = RqNTT::ONE;

        assert_eq!(one.icrt(), RqPoly::ONE)
    }

    #[test]
    fn test_nonresidue_is_order_24() {
        let nonresidue: Fq = BabyBear3ExtConfig::NONRESIDUE;

        let mut pow = nonresidue;

        for _i in 0..22 {
            pow *= nonresidue;
            assert_ne!(pow, <Fq as Field>::ONE);
        }

        assert_eq!(pow * nonresidue, <Fq as Field>::ONE);
    }

    #[test]
    fn test_reduce() {
        let mut rng = ChaCha8Rng::seed_from_u64(0);
        let mut coeffs: Vec<_> = (0..(2 * ntt::D)).map(|_| Fq::rand(&mut rng)).collect();

        let poly = DensePolynomial::from_coefficients_slice(&coeffs);

        // X^72 - X^36 + 1
        let the_cyclotomic = SparsePolynomial::from_coefficients_slice(&[
            (72, <Fq as Field>::ONE),
            (36, -<Fq as Field>::ONE),
            (0, <Fq as Field>::ONE),
        ]);

        let (_, r) =
            DenseOrSparsePolynomial::divide_with_q_and_r(&poly.into(), &the_cyclotomic.into())
                .unwrap();

        BabyBearRingConfig::reduce_in_place(&mut coeffs);

        assert_eq!(r.coeffs, coeffs);
    }

    #[test]
    fn test_mul_crt() {
        let mut rng = ChaCha8Rng::seed_from_u64(0);
        let coeff_1 = RqPoly::rand(&mut rng);
        let coeff_2 = RqPoly::rand(&mut rng);

        let ntt_form_1: RqNTT = coeff_1.crt();
        let ntt_form_2: RqNTT = coeff_2.crt();

        let ntt_mul = ntt_form_1 * ntt_form_2;
        let coeffs_mul = coeff_1 * coeff_2;

        // ntt_mul.coeffs() performs INTT while coeffs_mul.coeffs() just returns the coefficients
        assert_eq!(ntt_mul.icrt(), coeffs_mul);
    }

    #[test]
    fn test_cyclotomic() {
        let mut rng = ChaCha8Rng::seed_from_u64(0);

        let mut a = RqPoly::rand(&mut rng);
        let initial_a = a;

        let x = RqPoly::x();

        for i in 1..RqPoly::dimension() {
            a.rot();
            assert_eq!(a, initial_a * x.pow([i as u64]));
        }
    }
}
