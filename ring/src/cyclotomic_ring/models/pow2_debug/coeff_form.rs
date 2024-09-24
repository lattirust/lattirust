use ark_ff::Field;

use crate::{
    balanced_decomposition::{decompose_balanced_polyring, Decompose},
    cyclotomic_ring::{models::pow2_debug::Fp64Pow2, CyclotomicPolyRingGeneral},
    traits::{WithL2Norm, WithLinfNorm},
    PolyRing, WithRot,
};

use super::Pow2Rp64Config;

pub type Pow2CyclotomicPolyRing<const Q: u64, const PHI_D: usize> =
    CyclotomicPolyRingGeneral<Pow2Rp64Config<Q, PHI_D>, 1, PHI_D>;

impl<const Q: u64, const PHI_D: usize> Decompose for Pow2CyclotomicPolyRing<Q, PHI_D> {
    fn decompose(&self, b: u128, padding_size: Option<usize>) -> Vec<Self> {
        decompose_balanced_polyring(self, b, padding_size)
    }
}
impl<const Q: u64, const PHI_D: usize> WithL2Norm for Pow2CyclotomicPolyRing<Q, PHI_D> {
    fn l2_norm_squared(&self) -> u128 {
        self.coeffs().l2_norm_squared()
    }
}

impl<const Q: u64, const PHI_D: usize> WithLinfNorm for Pow2CyclotomicPolyRing<Q, PHI_D> {
    fn linf_norm(&self) -> u128 {
        self.coeffs().linf_norm()
    }
}

impl<const Q: u64, const PHI_D: usize> WithRot for Pow2CyclotomicPolyRing<Q, PHI_D> {
    fn multiply_by_xi(&self, i: usize) -> Vec<Self::BaseRing> {
        let bs = self.0;
        let len = bs.ncols();
        assert_eq!(len, PHI_D);
        let mut result = vec![<Fp64Pow2<Q, PHI_D> as Field>::ZERO; len];
        for (j, &coeff) in bs.iter().enumerate() {
            if j + i < len {
                result[(j + i) % len] += coeff;
            } else {
                result[(j + i) % len] -= coeff;
            }
        }
        result
    }
}
