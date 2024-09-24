use std::ops::Mul;

use ark_ff::{Field, Fp64};

use crate::{
    cyclotomic_ring::{CyclotomicPolyRingNTTGeneral, RpConfig},
    traits::{WithL2Norm, WithLinfNorm},
    PolyRing, WithRot,
};

use super::Pow2Rp64Config;

pub type Fp64Pow2<const Q: u64, const PHI_D: usize> =
    Fp64<<Pow2Rp64Config<Q, PHI_D> as RpConfig<1>>::FpConfig>; // This looks ugly change this

/// A cyclotomic ring with a cyclotomic polynomial degree of a power of two
/// and a modulus less than 2^64
pub type Pow2CyclotomicPolyRingNTT<const Q: u64, const PHI_D: usize> =
    CyclotomicPolyRingNTTGeneral<Pow2Rp64Config<Q, PHI_D>, 1, PHI_D>;

impl<const Q: u64, const PHI_D: usize> WithL2Norm for Pow2CyclotomicPolyRingNTT<Q, PHI_D> {
    fn l2_norm_squared(&self) -> u128 {
        self.coeffs().l2_norm_squared()
    }
}

impl<const Q: u64, const PHI_D: usize> WithLinfNorm for Pow2CyclotomicPolyRingNTT<Q, PHI_D> {
    fn linf_norm(&self) -> u128 {
        self.coeffs().linf_norm()
    }
}

impl<const Q: u64, const PHI_D: usize> WithRot for Pow2CyclotomicPolyRingNTT<Q, PHI_D> {
    fn multiply_by_xi(&self, i: usize) -> Vec<Self::BaseRing> {
        let mut xi = if i < PHI_D {
            vec![<Fp64Pow2<Q, PHI_D> as Field>::ZERO; PHI_D]
        } else {
            vec![<Fp64Pow2<Q, PHI_D> as Field>::ZERO; i + 1] // Shouldn't get to this point for the
                                                             // purpose of RotSum
        };
        xi[i] = <Fp64Pow2<Q, PHI_D> as Field>::ONE;

        let xi_poly = Pow2CyclotomicPolyRingNTT::from(xi);
        let result = (*self * xi_poly).0;
        result.iter().copied().collect::<Vec<_>>()
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
