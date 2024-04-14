use ark_ff::UniformRand;
use ark_std::rand;
use delegate::delegate;
use nalgebra::{
    self, Const, DefaultAllocator, Dim, Dyn, Owned, RawStorage, VecStorage, ViewStorage,
};
use nalgebra::allocator::Allocator;
use num_traits::Zero;

use crate::linear_algebra::generic_matrix::GenericMatrix;
use crate::linear_algebra::Scalar;

pub type GenericVector<T, R, S> = GenericMatrix<T, R, Const<1>, S>;
pub type Vector<T> = GenericVector<T, Dyn, VecStorage<T, Dyn, Const<1>>>;

impl<T: Scalar> From<Vec<T>> for Vector<T> {
    fn from(v: Vec<T>) -> Self {
        nalgebra::DVector::from(v).into()
    }
}

impl<T: Scalar, RStride: Dim, CStride: Dim>
    From<nalgebra::Matrix<T, Dyn, Const<1>, ViewStorage<'_, T, Dyn, Const<1>, RStride, CStride>>>
    for Vector<T>
{
    fn from(
        v: nalgebra::Matrix<T, Dyn, Const<1>, ViewStorage<'_, T, Dyn, Const<1>, RStride, CStride>>,
    ) -> Self {
        v.into_owned().into()
    }
}

impl<T: Scalar> Vector<T> {
    delegate! {
        to self.0 {
            pub fn len(&self) -> usize;
            pub fn as_slice(&self) -> &[T];
        }
    }
}
impl<T, S> GenericVector<T, Dyn, S>
where
    T: Scalar,
    S: RawStorage<T, Dyn, Const<1>>,
    DefaultAllocator: Allocator<T, Dyn, Const<1>>,
{
    pub type VectorBuffer = <DefaultAllocator as Allocator<T, Dyn, Const<1>>>::Buffer;
    pub fn from_fn<F: FnMut(usize, usize) -> T>(
        n: usize,
        f: F,
    ) -> GenericVector<T, Dyn, Self::VectorBuffer> {
        nalgebra::Vector::<T, Dyn, Self::VectorBuffer>::from_fn(n, f).into()
    }

    pub fn from_vec(v: Vec<T>) -> GenericVector<T, Dyn, Self::VectorBuffer> {
        nalgebra::Vector::<T, Dyn, Self::VectorBuffer>::from_vec(v).into()
    }

    pub fn from_element(n: usize, element: T) -> GenericVector<T, Dyn, Self::VectorBuffer> {
        nalgebra::Vector::<T, Dyn, Self::VectorBuffer>::from_element(n, element).into()
    }

    pub fn from_slice(data: &[T]) -> GenericVector<T, Dyn, Owned<T, Dyn, Const<1>>> {
        nalgebra::Vector::<T, Dyn, Owned<T, Dyn, Const<1>>>::from_row_slice(data).into()
    }
}

impl<T: Scalar + Zero> Vector<T> {
    pub fn zeros(n: usize) -> Self {
        Self::Inner::zeros(n).into()
    }
}

pub type GenericRowVector<T, C, S> = GenericMatrix<T, Const<1>, C, S>;
pub type RowVector<T> = GenericRowVector<T, Dyn, VecStorage<T, Const<1>, Dyn>>;

impl<T: Scalar> From<Vec<T>> for RowVector<T> {
    fn from(v: Vec<T>) -> Self {
        nalgebra::RowDVector::from(v).into()
    }
}

impl<T: Scalar> RowVector<T> {
    delegate! {
        to self.0 {
            pub fn len(&self) -> usize;
            pub fn as_slice(&self) -> &[T];
         }
    }
}

impl<T: UniformRand + Scalar> Vector<T> {
    pub fn rand<Rng: rand::Rng + ?Sized>(n: usize, rng: &mut Rng) -> Self {
        Self::from_fn(n, |_, _| T::rand(rng))
    }
}
