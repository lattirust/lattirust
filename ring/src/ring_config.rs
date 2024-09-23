// RpConfig

use ark_ff::{BigInt, Field, Fp, FpConfig, MontBackend};
use ark_poly::{EvaluationDomain, Radix2EvaluationDomain};
use num_traits::Zero;

use crate::{ntt, z_q::FqConfig, Zq};

/// A configuration trait for fully splitting cyclotomic rings.
/// Used for setting up CRT data compile-time.
pub trait RpConfig<const N: usize> {
    /// The FpConfig of the underlying Fp field.
    type FpConfig: FpConfig<N>;
    /// The type of evaluation domain to compute the NTT on.
    type EvaluationDomain: EvaluationDomain<Fp<Self::FpConfig, N>>;

    /// Constant instance of the evaluation domain.
    const DOMAIN: Self::EvaluationDomain;

    /// The parameter of the cyclotomic polynomial.
    /// The ring determined by the configuration is meant to be
    /// Fp[X]/(Phi_D(X)).
    const D: usize;
    /// Euler's totient of D: \varphi(D). I.e. number of components
    const PHI_D: usize;
    /// An element of Fp of order D.
    const PRIMITIVE_DTH_ROOT_OF_UNITY: Fp<Self::FpConfig, N>;

    /// Given an index i\in [0,..., D-1] tells if its coprime with D,
    /// which implies that \psi^i is again a primitive Dth root of unity
    /// if \psi is.
    fn is_primitive_index(i: usize) -> bool;

    /// Given coefficients of a polynomial of degree D
    /// reduces it mod Phi_D(X).
    fn reduce_in_place(coefficients: &mut Vec<Fp<Self::FpConfig, N>>);

    /// Computes the evaluations of the polynomial with
    /// coefficients `coefficients` on the primitive Dth roots of unity.
    fn crt_in_place(coefficients: &mut Vec<Fp<Self::FpConfig, N>>) {
        coefficients.resize(Self::PHI_D, Fp::<Self::FpConfig, N>::ZERO);

        // We resize the coefficient vector with D - PHI_D zeros
        // to have a polynomial of "degree" D.
        coefficients.extend((0..(Self::D - Self::PHI_D)).map(|_| Fp::<Self::FpConfig, N>::ZERO));
        Self::DOMAIN.fft_in_place(coefficients);

        // Once we've done the NTT we remove the evaluations at PRIMITIVE_INDICES.
        // Those evaluations correspond to the evaluations at non-primitive roots of unity.
        let evals: Vec<Fp<Self::FpConfig, N>> = coefficients
            .iter()
            .enumerate()
            .filter_map(|(i, &x)| Self::is_primitive_index(i).then_some(x))
            .collect();

        *coefficients = evals
    }

    /// Given primitive Dth root of unity evaluations of a polynomial
    /// computes its coefficient form.
    fn icrt_in_place(evaluations: &mut Vec<Fp<Self::FpConfig, N>>) {
        // Enlarge the evaluation vector to be of length D.
        // Add dummy zero evaluations at the non primitive indices.
        // let mut evals = evaluations.clone();
        let mut evals = Vec::new();
        evals.resize(Self::D, Fp::<Self::FpConfig, N>::ZERO);

        for (i, value) in evals.iter_mut().enumerate() {
            if Self::is_primitive_index(i) {
                *value = evaluations.remove(0);
            }
        }
        Self::DOMAIN.ifft_in_place(&mut evals);

        // Reduce the resulting polynomial mod Phi_D(X).
        Self::reduce_in_place(&mut evals);
        *evaluations = evals;
    }
}

/// A generic configuration for power-of-two fully splitting cyclotomic rings
/// of small modulus (< 2^64).
pub struct Pow2Rp64Config<const Q: u64, const PHI_D: usize>;

impl<const Q: u64, const PHI_D: usize> RpConfig<1> for Pow2Rp64Config<Q, PHI_D> {
    type FpConfig = MontBackend<FqConfig<Q>, 1>; // TODO: change FqConfig to a suitable FpConfig
                                                 // that admits multiple primes

    type EvaluationDomain = Radix2EvaluationDomain<Zq<Q>>;

    const DOMAIN: Self::EvaluationDomain = Radix2EvaluationDomain {
        size: { 2 * PHI_D as u64 },
        log_size_of_group: {
            assert!(
                PHI_D.is_power_of_two(),
                "PHI_D most be a power of two to use Pow2Rp64Config"
            );
            if 2 * PHI_D == 0 {
                0
            } else {
                1usize.leading_zeros() - (2 * PHI_D).leading_zeros()
            }
        },
        size_as_field_element: { Fp::new(BigInt::new([2 * PHI_D as u64])) },
        size_inv: {
            let size = 2 * PHI_D as u64;

            Fp::new(BigInt::new([ntt::const_pow_mod::<Q>(size, Q - 2)]))
        },
        group_gen: {
            let two_adic_root = Self::FpConfig::TWO_ADIC_ROOT_OF_UNITY.0 .0[0];
            let two_adicity = { Q - 1 };

            Fp::new(BigInt::new([ntt::const_pow_mod::<Q>(
                two_adic_root,
                two_adicity / (2 * PHI_D as u64),
            )]))
        },
        group_gen_inv: {
            let two_adic_root = Self::FpConfig::TWO_ADIC_ROOT_OF_UNITY.0 .0[0];
            let two_adicity = { Q - 1 };

            let the_root = ntt::const_pow_mod::<Q>(two_adic_root, two_adicity / (2 * PHI_D as u64));

            Fp::new(BigInt::new([ntt::const_pow_mod::<Q>(the_root, Q - 2)]))
        },
        offset: <Fp<Self::FpConfig, 1> as Field>::ONE,
        offset_inv: <Fp<Self::FpConfig, 1> as Field>::ONE,
        offset_pow_size: <Fp<Self::FpConfig, 1> as Field>::ONE,
    };

    const D: usize = 2 * PHI_D;

    const PHI_D: usize = PHI_D;

    const PRIMITIVE_DTH_ROOT_OF_UNITY: Fp<Self::FpConfig, 1> = {
        let two_adic_root = Self::FpConfig::TWO_ADIC_ROOT_OF_UNITY.0 .0[0];
        let two_adicity = { Q - 1 };

        Fp::new(BigInt::new([ntt::const_pow_mod::<Q>(
            two_adic_root,
            two_adicity / (2 * PHI_D as u64),
        )]))
    };

    fn is_primitive_index(i: usize) -> bool {
        i % 2 == 1
    }

    fn reduce_in_place(coefficients: &mut Vec<Fp<Self::FpConfig, 1>>) {
        coefficients.resize(2 * PHI_D, Fp::zero());
        let (left, right) = coefficients.split_at_mut(PHI_D);
        for (coeff_left, coeff_right) in left.iter_mut().zip(right) {
            *coeff_left -= coeff_right;
        }

        // Truncate the resulting vector to be of size D.
        coefficients.resize(PHI_D, Fp::zero());
    }
}
