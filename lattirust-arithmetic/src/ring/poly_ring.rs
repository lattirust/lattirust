use std::ops::Mul;

use ark_ff::{Field, Fp, FpConfig};
use ark_serialize::CanonicalSerialize;
use num_traits::Zero;

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
}

impl<C: FpConfig<N>, const N: usize> FromRandomBytes<Fp<C, N>> for Fp<C, N> {
    fn byte_size() -> usize {
        Self::zero().uncompressed_size() + 9 // TODO: check if this is correct; this is inferred from Fp<C, N>::from_random_bytes()
    }

    fn try_from_random_bytes(bytes: &[u8]) -> Option<Self> {
        Self::from_random_bytes(&mut bytes.as_ref())
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
        Into::<SignedRepresentative>::into(*self).0.abs() as u128
    }
}
