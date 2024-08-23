use std::ops::Mul;

use crate::linear_algebra::Vector;
use crate::ring::representatives::{SignedRepresentative, UnsignedRepresentative};
use crate::ring::Ring;
use crate::traits::{FromRandomBytes, WithConjugationAutomorphism, WithL2Norm, WithLinfNorm};

pub trait ConvertibleRing:
    Ring + Into<UnsignedRepresentative> + Into<SignedRepresentative> + From<SignedRepresentative>
{
}

pub trait PolyRing:
    Ring
    + Mul<Self::BaseRing, Output = Self>
    + From<Vec<Self::BaseRing>>
    + WithConjugationAutomorphism
    + WithL2Norm
    + WithLinfNorm
    + FromRandomBytes<Self>
    + From<u128>
    + From<Self::BaseRing>
{
    type BaseRing: ConvertibleRing;

    fn coeffs(&self) -> Vec<Self::BaseRing>;
    fn flattened(vec: &Vector<Self>) -> Vector<Self::BaseRing> {
        Self::flattened_coeffs(vec).into()
    }
    fn flattened_coeffs(vec: &Vector<Self>) -> Vec<Self::BaseRing> {
        vec.into_iter()
            .flat_map(|x| x.coeffs())
            .collect::<Vec<Self::BaseRing>>()
    }
    fn dimension() -> usize;

    fn from_scalar(scalar: Self::BaseRing) -> Self;

    #[inline]
    fn x() -> Self {
        Self::from(vec![Self::BaseRing::ZERO, Self::BaseRing::ONE])
    }
}

// Norms
impl<F: ConvertibleRing> WithL2Norm for F {
    fn l2_norm_squared(&self) -> u128 {
        Into::<SignedRepresentative>::into(*self).0.pow(2) as u128
    }
}

impl<F: ConvertibleRing> WithLinfNorm for F {
    fn linf_norm(&self) -> u128 {
        Into::<SignedRepresentative>::into(*self).0.unsigned_abs()
    }
}
