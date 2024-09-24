use ark_crypto_primitives::sponge::Absorb;
use ark_ff::PrimeField;

use crate::balanced_decomposition::{decompose_balanced, Decompose};
use crate::representatives::{SignedRepresentative, UnsignedRepresentative};
use crate::traits::{FromRandomBytes, WithL2Norm, WithLinfNorm};
use crate::Ring;
use lattirust_linear_algebra::{Matrix, Vector};

pub trait ConvertibleRing:
    Ring + Into<UnsignedRepresentative> + Into<SignedRepresentative> + From<SignedRepresentative>
{
}

pub trait PolyRing:
    Ring
    // + Mul<Self::BaseRing, Output = Self>
    + From<Vec<Self::BaseRing>>
    + FromRandomBytes<Self>
    + From<u128>
    // + From<Self::BaseRing>
{
    type BaseRing: Ring;

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

impl<R: ConvertibleRing> Decompose for R {
    fn decompose(&self, b: u128, padding_size: Option<usize>) -> Vec<Self> {
        decompose_balanced(self, b, padding_size)
    }
}

pub trait OverField: PolyRing<BaseRing: PrimeField + Absorb> {}

pub trait WithRot: PolyRing {
    fn rot(&self) -> Matrix<Self::BaseRing> {
        let degree = Self::dimension() - 1;
        let coeffs = self.coeffs();
        let mut columns = Vec::with_capacity(degree);

        for i in 0..degree {
            let vec_xi_a = if i == 0 {
                Vector::from_vec(coeffs.clone())
            } else {
                Vector::from_vec(self.multiply_by_xi(i))
            };
            columns.push(vec_xi_a);
        }

        Matrix::from_columns(columns.as_slice())
    }

    fn rot_sum(&self, bs: &[Self::BaseRing]) -> Vec<Self::BaseRing> {
        let degree = Self::dimension(); // if tau is 1 in latticefold paper lemma 2.1.
        let mut result = vec![Self::BaseRing::ZERO; degree];
        let rot_matrix = self.rot();

        for (i, b_i) in bs.iter().enumerate() {
            let vec_xi_a = rot_matrix.column(i);
            for (j, coeff) in vec_xi_a.iter().enumerate() {
                result[j] += *coeff * *b_i;
            }
        }
        result
    }

    // Multiply by x^i depends on the irreducible polynomial
    fn multiply_by_xi(&self, i: usize) -> Vec<Self::BaseRing>;
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
