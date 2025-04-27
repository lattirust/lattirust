use std::ops::Mul;

use crate::linear_algebra::Vector;
use crate::ring::Ring;
use crate::traits::{FromRandomBytes, WithConjugationAutomorphism, WithL2Norm, WithLinfNorm};

pub trait PolyRing:
    Ring
    + Mul<Self::BaseRing, Output = Self>
    + From<Vec<Self::BaseRing>>
    + WithConjugationAutomorphism
    + WithL2Norm
    + WithLinfNorm
    + FromRandomBytes<Self>
    + From<Self::BaseRing>
{
    type BaseRing: Ring;

    fn coefficients(&self) -> Vec<Self::BaseRing>;

    fn try_from_coefficients(coeffs: &[Self::BaseRing]) -> Option<Self>;

    fn flattened(vec: &Vector<Self>) -> Vector<Self::BaseRing> {
        Self::flattened_coeffs(vec).into()
    }

    fn flattened_coeffs(vec: &Vector<Self>) -> Vec<Self::BaseRing> {
        vec.into_iter()
            .flat_map(|x| x.coefficients())
            .collect::<Vec<Self::BaseRing>>()
    }

    fn dimension() -> usize;

    fn from_scalar(scalar: Self::BaseRing) -> Self;

    #[inline]
    fn x() -> Self {
        Self::from(vec![Self::BaseRing::ZERO, Self::BaseRing::ONE])
    }
}

#[macro_export]
macro_rules! test_polyring {
    ($T:ty, $N:expr) => {
        #[test]
        fn test_coefficients() {
            let rng = &mut ark_std::test_rng();
            for _ in 0..$N {
                let poly = <$T as UniformRand>::rand(rng);
                let coeffs = poly.coefficients();
                assert!(coeffs.len() <= <$T as PolyRing>::dimension());
                let poly_ = <$T as PolyRing>::try_from_coefficients(&coeffs).unwrap();
                assert_eq!(poly, poly_);
            }
        }
    };
}
