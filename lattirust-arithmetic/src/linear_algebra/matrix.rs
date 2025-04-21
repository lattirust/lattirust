#![allow(non_snake_case)]

use std::ops::Neg;

use ark_std::rand::prelude::SliceRandom;
use ark_std::{rand, UniformRand};
use delegate::delegate;
use nalgebra::{self, ArrayStorage, ComplexField, Dyn, VecStorage};
use num_traits::{One, Zero};
use rayon::prelude::*;

use crate::linear_algebra::generic_matrix::GenericMatrix;
use crate::linear_algebra::{RowVector, Vector};
use crate::linear_algebra::{Scalar, SymmetricMatrix};

pub type Const<const S: usize> = nalgebra::Const<S>;
pub type Matrix<T> = GenericMatrix<T, Dyn, Dyn, VecStorage<T, Dyn, Dyn>>;
pub type SMatrix<T, const R: usize, const C: usize> =
    GenericMatrix<T, Const<R>, Const<C>, ArrayStorage<T, R, C>>;

impl<R: ComplexField> Matrix<R> {
    delegate! {
        to self.0 {
            #[into]
            pub fn symmetric_eigenvalues(&self) -> Vector<R::RealField>;
        }
    }
}

impl<T: Scalar> Matrix<T> {
    pub fn from_vec(m: usize, n: usize, data: Vec<T>) -> Self {
        Self::Inner::from_vec(m, n, data).into()
    }

    pub fn from_fn(m: usize, n: usize, f: impl FnMut(usize, usize) -> T) -> Self {
        Self::Inner::from_fn(m, n, f).into()
    }

    pub fn from_rows(rows: &[RowVector<T>]) -> Self {
        Self::Inner::from_rows(
            rows.iter()
                .cloned()
                .map(|row| row.0)
                .collect::<Vec<_>>()
                .as_slice(),
        )
        .into()
    }

    pub fn from_columns(columns: &[Vector<T>]) -> Self {
        Self::Inner::from_columns(
            columns
                .iter()
                .cloned()
                .map(|row| row.0)
                .collect::<Vec<_>>()
                .as_slice(),
        )
        .into()
    }
}

impl<R: Scalar + Zero> Matrix<R> {
    pub fn zeros(m: usize, n: usize) -> Self {
        Self::Inner::zeros(m, n).into()
    }
    pub fn identity(m: usize, n: usize) -> Self
    where
        R: One,
    {
        Self::Inner::identity(m, n).into()
    }
}

impl<T: Scalar> IntoIterator for Matrix<T>
where
    nalgebra::DMatrix<T>: IntoIterator,
{
    type Item = <nalgebra::DMatrix<T> as IntoIterator>::Item;
    type IntoIter = <nalgebra::DMatrix<T> as IntoIterator>::IntoIter;
    fn into_iter(self) -> Self::IntoIter {
        self.0.into_iter()
    }
}

impl<T: Scalar + UniformRand> Matrix<T> {
    pub fn rand<Rng: rand::Rng + ?Sized>(m: usize, n: usize, rng: &mut Rng) -> Self {
        Self::from_fn(m, n, |_, _| T::rand(rng))
    }

    pub fn par_rand(m: usize, n: usize) -> Self
    where
        T: Send + Sync,
    {
        let data = (0..m * n)
            .into_par_iter()
            .map_init(rand::thread_rng, |mut rng, _| T::rand(&mut rng))
            .collect();
        Self::from_vec(m, n, data)
    }

    pub fn rand_symmetric<Rng: rand::Rng + ?Sized>(n: usize, rng: &mut Rng) -> Self {
        SymmetricMatrix::<T>::rand(n, rng).into()
    }
}

impl<T: Scalar + UniformRand + Zero + One + Neg<Output = T>> Matrix<T> {
    pub fn rand_ternary<Rng: rand::Rng + ?Sized>(m: usize, n: usize, rng: &mut Rng) -> Self {
        Self::from_fn(m, n, |_, _| {
            [-T::one(), T::zero(), T::one()]
                .choose(rng)
                .unwrap()
                .clone()
        })
    }
}

impl<T: Scalar + UniformRand, const R: usize, const C: usize> UniformRand for SMatrix<T, R, C> {
    fn rand<Rng: rand::Rng + ?Sized>(rng: &mut Rng) -> Self {
        nalgebra::SMatrix::<T, R, C>::from_fn(|_, _| T::rand(rng)).into()
    }
}

#[allow(non_snake_case)]
#[cfg(test)]
mod tests {
    use ark_std::test_rng;

    use crate::ring::pow2_cyclotomic_poly_ring::Pow2CyclotomicPolyRing;
    use crate::ring::Zq1;

    use super::*;

    type R = Pow2CyclotomicPolyRing<Zq1<3>, 20>;

    #[test]
    fn test_sample_uniform() {
        let m = 10;
        let n = 20;
        let rng = &mut test_rng();
        let A = Matrix::<R>::rand(m, n, rng);
        assert_eq!(A.nrows(), m);
        assert_eq!(A.ncols(), n);
    }
}
