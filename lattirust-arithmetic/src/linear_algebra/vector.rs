use ark_ff::UniformRand;
use ark_std::rand;
use ark_std::rand::Rng;
use delegate::delegate;
use nalgebra::allocator::Allocator;
use nalgebra::{
    self, ArrayStorage, Const, DefaultAllocator, Dim, Dyn, Owned, RawStorage, VecStorage,
    ViewStorage,
};
use num_traits::Zero;

use crate::linear_algebra::generic_matrix::GenericMatrix;
use crate::linear_algebra::Scalar;
use crate::ring::SignedRepresentative;
use crate::traits::{WithL2Norm, WithLinfNorm};

pub type GenericVector<T, R, S> = GenericMatrix<T, R, Const<1>, S>;
pub type Vector<T> = GenericVector<T, Dyn, VecStorage<T, Dyn, Const<1>>>;
pub type SVector<T, const N: usize> = GenericVector<T, Const<N>, ArrayStorage<T, N, 1>>;

impl<T: Scalar> From<Vec<T>> for Vector<T> {
    fn from(v: Vec<T>) -> Self {
        nalgebra::DVector::from(v).into()
    }
}

impl<T: Scalar, const N: usize, E> TryFrom<Vec<T>> for SVector<T, N>
where
    [T; N]: TryFrom<Vec<T>, Error = E>,
{
    type Error = E;

    fn try_from(value: Vec<T>) -> Result<Self, Self::Error> {
        Ok(Self::const_from_array(value.try_into()?))
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

impl<T: Scalar, const N: usize> SVector<T, N> {
    pub(crate) const fn const_from_array(array: [T; N]) -> Self {
        Self(nalgebra::SVector::<T, N>::from_array_storage(
            ArrayStorage::<T, { N }, 1>([array; 1]),
        ))
    }
}

impl<T: Scalar + Zero> Vector<T> {
    pub fn zeros(n: usize) -> Self {
        Self::Inner::zeros(n).into()
    }
}

impl<T: UniformRand + Scalar> Vector<T> {
    pub fn rand<Rng: rand::Rng + ?Sized>(n: usize, rng: &mut Rng) -> Self {
        Self::from_fn(n, |_, _| T::rand(rng))
    }

    pub fn rand_vector_with_bounded_norm(n: usize, norm_bound: i128, rng: &mut impl Rng) -> Self
    where
        T: From<SignedRepresentative>,
    {
        loop {
            let mut vec = Vector::<f64>::rand(n, rng);
            let scale = rng.gen_range(0f64..(norm_bound as f64 / vec.0.norm()));
            vec *= scale;

            vec = vec.map(|x| x.round());
            if vec.0.norm() <= norm_bound as f64 {
                return vec.map(|x| T::from(SignedRepresentative(x as i128)));
            }
        }
    }
}

impl<T: Scalar + WithL2Norm, R: Dim, S: RawStorage<T, R, Const<1>>> WithL2Norm
    for GenericVector<T, R, S>
{
    fn l2_norm_squared(&self) -> u128 {
        self.into_iter()
            .cloned()
            .collect::<Vec<_>>()
            .l2_norm_squared()
    }
}

impl<T: Scalar + WithLinfNorm, R: Dim, S: RawStorage<T, R, Const<1>>> WithLinfNorm
    for GenericVector<T, R, S>
{
    fn linf_norm(&self) -> u128 {
        self.into_iter().cloned().collect::<Vec<_>>().linf_norm()
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
