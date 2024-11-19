use ark_ff::Fp64;
use ark_std::ops::Mul;
use num_bigint::BigUint;

use super::{Pow2CyclotomicPolyRing, Pow2Rp64Config};
use crate::{
    cyclotomic_ring::{CyclotomicConfig, CyclotomicPolyRingNTTGeneral, RpConfig, CRT, ICRT},
    traits::{WithL2Norm, WithLinfNorm},
    OverField, PolyRing,
};

pub type Fp64Pow2<const Q: u64, const PHI_D: usize> =
    Fp64<<Pow2Rp64Config<Q, PHI_D> as RpConfig<1>>::FpConfig>; // This looks ugly change this

/// A cyclotomic ring with a cyclotomic polynomial degree of a power of two
/// and a modulus less than 2^64
pub type Pow2CyclotomicPolyRingNTT<const Q: u64, const PHI_D: usize> =
    CyclotomicPolyRingNTTGeneral<Pow2Rp64Config<Q, PHI_D>, 1, PHI_D>;

impl<const Q: u64, const PHI_D: usize> WithL2Norm for Pow2CyclotomicPolyRingNTT<Q, PHI_D> {
    fn l2_norm_squared(&self) -> BigUint {
        self.coeffs().l2_norm_squared()
    }
}

impl<const Q: u64, const PHI_D: usize> WithLinfNorm for Pow2CyclotomicPolyRingNTT<Q, PHI_D> {
    fn linf_norm(&self) -> BigUint {
        self.coeffs().linf_norm()
    }
}

impl<const Q: u64, const PHI_D: usize> Mul<Fp64Pow2<Q, PHI_D>>
    for Pow2CyclotomicPolyRingNTT<Q, PHI_D>
{
    type Output = Self;

    fn mul(self, rhs: Fp64Pow2<Q, PHI_D>) -> Self::Output {
        self.mul(Self::from_scalar(rhs))
    }
}

impl<const Q: u64, const PHI_D: usize> OverField for Pow2CyclotomicPolyRingNTT<Q, PHI_D> {}

impl<const Q: u64, const PHI_D: usize> From<Fp64Pow2<Q, PHI_D>>
    for Pow2CyclotomicPolyRingNTT<Q, PHI_D>
{
    fn from(value: Fp64Pow2<Q, PHI_D>) -> Self {
        Self::from_scalar(value)
    }
}

impl<const Q: u64, const PHI_D: usize> CRT for Pow2CyclotomicPolyRing<Q, PHI_D> {
    type CRTForm = Pow2CyclotomicPolyRingNTT<Q, PHI_D>;

    fn crt(mut self) -> Pow2CyclotomicPolyRingNTT<Q, PHI_D> {
        <Pow2Rp64Config<Q, PHI_D> as CyclotomicConfig<1>>::crt_in_place(&mut self.0);

        Pow2CyclotomicPolyRingNTT::from_array(self.0)
    }
}

impl<const Q: u64, const PHI_D: usize> ICRT for Pow2CyclotomicPolyRingNTT<Q, PHI_D> {
    type ICRTForm = Pow2CyclotomicPolyRing<Q, PHI_D>;

    fn icrt(mut self) -> Pow2CyclotomicPolyRing<Q, PHI_D> {
        <Pow2Rp64Config<Q, PHI_D> as CyclotomicConfig<1>>::icrt_in_place(&mut self.0);

        Pow2CyclotomicPolyRing::from_coeffs_vec(self.0.to_vec())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cyclotomic_ring::{crt::CRT, models::pow2_debug::Pow2CyclotomicPolyRing, ICRT};
    use ark_std::UniformRand;
    use rand::SeedableRng;
    use rand_chacha::ChaCha8Rng;

    const FERMAT_Q: u64 = 65537;

    #[test]
    fn test_crt_conversion() {
        let mut rng = ChaCha8Rng::seed_from_u64(0);

        // Test for a few different sizes
        test_crt_size::<8>(&mut rng);
        test_crt_size::<16>(&mut rng);
        test_crt_size::<32>(&mut rng);
    }

    fn test_crt_size<const SIZE: usize>(rng: &mut impl rand::Rng) {
        // Create random polynomial in coefficient form
        let coeff_1 = Pow2CyclotomicPolyRing::<FERMAT_Q, SIZE>::rand(rng);

        // Convert to NTT form using CRT trait
        let ntt_form: Pow2CyclotomicPolyRingNTT<FERMAT_Q, SIZE> = coeff_1.crt();

        // Convert back using ICRT trait
        let recovered: Pow2CyclotomicPolyRing<FERMAT_Q, SIZE> = ntt_form.icrt();

        // Verify round-trip conversion preserves the polynomial
        assert_eq!(coeff_1, recovered);
    }
}
