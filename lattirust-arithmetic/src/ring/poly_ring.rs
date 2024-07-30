use std::ops::Mul;

use ark_ff::{Field, Fp, FpConfig};
use ark_serialize::CanonicalSerialize;
use num_traits::Zero;

use crate::linear_algebra::{Matrix, Vector};
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

pub trait WithRot: PolyRing {
    fn rot(&self) -> Matrix<Self::BaseRing>;

    fn rot_sum(&self, bs: &Vec<Self::BaseRing>) -> Vec<Self::BaseRing> {
        let degree = Self::dimension(); // if tau is 1 in latticefold paper lemma 2.1.
        let mut result = vec![Self::BaseRing::ZERO; degree];
        let rot_matrix = self.rot();

        for (i, b_i) in bs.iter().enumerate() {
            let vec_xi_a = rot_matrix.column(i);
            for (j, coeff) in vec_xi_a.iter().enumerate() {
                result[j] = result[j] + (*coeff * *b_i);
            }
        }
        result
    }

    // Multiply by x^i depends on the irreducible polynomial
    fn multiply_by_xi(bs: &Vec<Self::BaseRing>, i: usize) -> Vec<Self::BaseRing>;
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
