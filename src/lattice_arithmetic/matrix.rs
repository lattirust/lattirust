#![allow(non_snake_case)]

use ark_ff::Field;
use ark_std::UniformRand;
use nalgebra::Scalar;
use num_traits::Zero;
use rand::Rng;
use serde::{Deserialize, Serialize};

use crate::lattice_arithmetic::poly_ring::PolyRing;
use crate::lattice_arithmetic::ring::Ring;
use crate::lattice_arithmetic::traits::WithL2Norm;

pub type Matrix<R> = nalgebra::DMatrix<R>;

pub type SparseMatrix<R> = nalgebra_sparse::CscMatrix<R>; // We typically have more rows than columns, hence CSC.

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct SymmetricMatrix<F: Clone>(Vec<Vec<F>>);

impl<F: Clone> From<Vec<Vec<F>>> for SymmetricMatrix<F> {
    fn from(value: Vec<Vec<F>>) -> Self {
        assert!(value.iter().enumerate().all(|(i, v_i)| v_i.len() == i + 1), "cannot convert value: Vec<Vec<F>> to SymmetricMatrix<F>, row has wrong number of entries");
        Self(value)
    }
}

impl<F: Clone + Scalar> From<Matrix<F>> for SymmetricMatrix<F> {
    fn from(value: Matrix<F>) -> Self {
        assert_eq!(value.transpose(), value);
        Self(
            value.row_iter().enumerate().map(
                |(i, v_i)|
                    v_i.columns(0, i + 1).into_iter().map(|v| v.clone()).collect()
            ).collect()
        )
    }
}

impl<F: Clone + Scalar> Into<Matrix<F>> for SymmetricMatrix<F> {
    fn into(self) -> Matrix<F> {
        Matrix::<F>::from_fn(self.size(), self.size(), |i, j| self.at(i, j).clone())
    }
}

impl<F: Zero + Clone> SymmetricMatrix<F> {
    pub fn zero(n: usize) -> SymmetricMatrix<F> {
        SymmetricMatrix::<F>((0..n).map(|i| vec![F::zero(); i + 1]).collect())
    }
}

impl<F: Clone + Eq> PartialEq<Matrix<F>> for SymmetricMatrix<F> {
    fn eq(&self, other: &Matrix<F>) -> bool {
        self.0.iter().enumerate().all(
            |(i, self_i)|
                self_i.into_iter().enumerate().all(
                    |(j, self_ij)| other[(i, j)] == *self_ij
                )
        )
    }
}

impl<F: Clone> SymmetricMatrix<F> {
    #[inline]
    pub fn size(&self) -> usize { self.0.len() }

    #[inline]
    pub fn at(&self, i: usize, j: usize) -> &F {
        debug_assert!(i < self.0.len() && j < self.0.len());
        if i <= j { &self.0[i][j] } else { &self.0[j][i] }
    }
    #[inline]
    pub fn at_mut(&mut self, i: usize, j: usize) -> &mut F {
        debug_assert!(i < self.0.len() && j < self.0.len());
        if i <= j { &mut self.0[i][j] } else { &mut self.0[j][i] }
    }

    pub fn diag(&self) -> Vec<F> {
        (0..self.size()).map(|i| self.at(i, i).clone()).collect()
    }

    pub fn rows(&self) -> &Vec<Vec<F>> { &self.0 }

    pub fn map<T, M>(&self, func: M) -> SymmetricMatrix<T>
        where T: Clone, M: Fn(&F) -> T
    {
        SymmetricMatrix::<T>::from(
            self.rows().into_iter().map(
                |row| row.into_iter().map(&func).collect()
            ).collect::<Vec<Vec<T>>>()
        )
    }
}


// TODO: implement as Mul trait for Vector<R> so that left-multiplication with a scalar is possible
pub type Vector<R> = nalgebra::DVector<R>;

pub fn sample_uniform_mat<R: Scalar + UniformRand, Rng: rand::Rng + ?Sized>(m: usize, n: usize, rng: &mut Rng) -> Matrix<R> {
    Matrix::<R>::from_fn(m, n, |_, _| R::rand(rng))
}

pub fn sample_uniform_mat_symmetric<R: Ring, Rng: rand::Rng + ?Sized>(m: usize, n: usize, rng: &mut Rng) -> Matrix<R> {
    let mut A = Matrix::<R>::zeros(m, n);
    for i in 0..m {
        for j in i..n {
            A[(i, j)] = R::rand(rng);
            A[(j, i)] = A[(i, j)]
        }
    }
    A
}

pub fn sample_uniform_vec<R: Ring, Rng: rand::Rng + ?Sized>(n: usize, rng: &mut Rng) -> Vector<R> {
    Vector::<R>::from_fn(n, |_, _| R::rand(rng))
}

#[allow(non_snake_case)]
#[cfg(test)]
mod tests {
    use ark_std::test_rng;

    use crate::lattice_arithmetic::pow2_cyclotomic_poly_ring::Pow2CyclotomicPolyRing;
    use crate::lattice_arithmetic::ring::Fq;

    use super::*;

    type R = Pow2CyclotomicPolyRing<Fq<3>, 20>;

    #[test]
    fn test_sample_uniform() {
        let m = 10;
        let n = 20;
        let rng = &mut test_rng();
        let A = sample_uniform_mat::<R, _>(m, n, rng);
        assert_eq!(A.nrows(), m);
        assert_eq!(A.ncols(), n);
    }
}