mod coeff_form;
mod ntt;
mod ntt_form;

use ark_ff::{BigInt, Field, Fp, FpConfig, MontBackend, Zero};
use ark_poly::Radix2EvaluationDomain;
use ark_std::vec::*;

use crate::{
    cyclotomic_ring::RpConfig,
    zn::z_q::{FqConfig, Zq},
};

pub use coeff_form::*;
pub use ntt_form::*;

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
