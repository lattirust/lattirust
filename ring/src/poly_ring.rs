use ark_crypto_primitives::sponge::Absorb;
use ark_ff::Field;
use ark_std::ops::Mul;

use crate::{traits::FromRandomBytes, Ring};
use lattirust_linear_algebra::{Matrix, Vector};

pub trait PolyRing: Ring + From<Vec<Self::BaseRing>> + FromRandomBytes<Self> + From<u128> {
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

pub trait OverField:
    PolyRing<BaseRing: Field<BasePrimeField: Absorb>>
    + Mul<Self::BaseRing, Output = Self>
    + From<Self::BaseRing>
{
}

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
