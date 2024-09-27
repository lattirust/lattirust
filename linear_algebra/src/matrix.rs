use ark_std::{ops::Neg, rand, rand::prelude::SliceRandom, One, UniformRand, Zero};
use delegate::delegate;
use nalgebra::{self, ComplexField, Dyn, VecStorage};
use rayon::prelude::*;

use crate::{generic_matrix::GenericMatrix, RowVector, Scalar, SymmetricMatrix, Vector};

pub type Matrix<T> = GenericMatrix<T, Dyn, Dyn, VecStorage<T, Dyn, Dyn>>;

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
