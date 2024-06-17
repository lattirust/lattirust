#![allow(non_snake_case)]

use std::hash::{Hash, Hasher};
use std::iter::Sum;
use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

use delegate::delegate;
use derive_more;
use derive_more::{Display, From, Index, IndexMut, Into};
use nalgebra::allocator::Allocator;
use nalgebra::constraint::{DimEq, ShapeConstraint};
use nalgebra::{
    self, Const, DefaultAllocator, Dim, DimMul, DimProd, DimRange, Owned, RawStorage, Scalar,
    Storage, StorageMut, ViewStorage,
};
use num_traits::{One, Zero};
use rayon::prelude::*;

use crate::linear_algebra::vector::{GenericRowVector, GenericVector};

pub trait ClosedAdd: nalgebra::ClosedAdd {}
impl<T> ClosedAdd for T where T: nalgebra::ClosedAdd {}
#[allow(dead_code)]
pub trait ClosedSub: nalgebra::ClosedSub {}
impl<T> ClosedSub for T where T: nalgebra::ClosedSub {}
pub trait ClosedMul: nalgebra::ClosedMul {}
impl<T> ClosedMul for T where T: nalgebra::ClosedMul {}

#[derive(Clone, Copy, Debug, Display, From, Into, Index, IndexMut)]
pub struct GenericMatrix<T: Scalar, R: Dim, C: Dim, S: RawStorage<T, R, C>>(
    pub(crate) nalgebra::Matrix<T, R, C, S>,
);

impl<T: Scalar, R: Dim, C: Dim, S: RawStorage<T, R, C>> GenericMatrix<T, R, C, S> {
    #[allow(dead_code)]
    pub(crate) type Inner = nalgebra::Matrix<T, R, C, S>;
}

impl<T: Scalar, R: Dim, C: Dim, S: RawStorage<T, R, C>> GenericMatrix<T, R, C, S> {
    delegate! {
        to self.0 {
            pub fn nrows(&self) -> usize;
            pub fn ncols(&self) -> usize;
            pub fn iter(&self) -> impl Iterator<Item=&'_ T>;
        }
    }
}

impl<T: Scalar, R: Dim, C: Dim, S: RawStorage<T, R, C>> GenericMatrix<T, R, C, S> {
    delegate! {
        to self.0 {
            #[into]
            pub fn columns_range<ColRange: DimRange<C>>(&self, range: ColRange) -> GenericMatrix<T, R, ColRange::Size, ViewStorage<'_,T, R, ColRange::Size, S::RStride, S::CStride>>;
        }
    }
}

impl<T: Scalar, R: Dim, C: Dim, S: RawStorage<T, R, C>> GenericMatrix<T, R, C, S> {
    delegate! {
        to self.0 {
            #[into]
            pub fn map<O: Scalar, F: FnMut(T) -> O>(&self, f: F) -> GenericMatrix<O, R, C, Owned<O, R, C>>
            where DefaultAllocator: Allocator<O, R, C>;
        }
    }

    pub fn row(
        &self,
        i: usize,
    ) -> GenericMatrix<
        T,
        Const<1>,
        C,
        ViewStorage<
            '_,
            T,
            Const<1>,
            C,
            <S as RawStorage<T, R, C>>::RStride,
            <S as RawStorage<T, R, C>>::CStride,
        >,
    > {
        self.0.row(i).into()
    }

    pub fn column(
        &self,
        i: usize,
    ) -> GenericMatrix<
        T,
        R,
        Const<1>,
        ViewStorage<
            '_,
            T,
            R,
            Const<1>,
            <S as RawStorage<T, R, C>>::RStride,
            <S as RawStorage<T, R, C>>::CStride,
        >,
    > {
        self.0.column(i).into()
    }

    pub fn transpose(&self) -> GenericMatrix<T, C, R, Owned<T, C, R>>
    where
        DefaultAllocator: Allocator<T, C, R>,
    {
        self.0.transpose().into()
    }
}

impl<T: Scalar, R: Dim, C: Dim, S: RawStorage<T, R, C>> GenericMatrix<T, R, C, S> {
    pub type RowViewStorage<'b> = ViewStorage<'b, T, Const<1>, C, S::RStride, S::CStride>;
    pub fn row_iter(
        &self,
    ) -> impl Iterator<Item = GenericRowVector<T, C, Self::RowViewStorage<'_>>> {
        self.0.row_iter().map(|r| r.into())
    }

    pub fn par_row_iter(
        &self,
    ) -> impl IndexedParallelIterator<Item = GenericRowVector<T, C, Self::RowViewStorage<'_>>>
    where
        T: Send + Sync,
    {
        // TODO: can we implement this without collecting into a vec first?
        self.row_iter().collect::<Vec<_>>().into_par_iter()
    }

    pub type ColumnViewStorage<'b> = ViewStorage<'b, T, R, Const<1>, S::RStride, S::CStride>;
    pub fn column_iter(
        &self,
    ) -> impl Iterator<Item = GenericVector<T, R, Self::ColumnViewStorage<'_>>> {
        self.0.column_iter().map(|c| c.into())
    }

    pub fn par_column_iter(
        &self,
    ) -> impl ParallelIterator<Item = GenericVector<T, R, Self::ColumnViewStorage<'_>>>
    where
        T: Send + Sync,
        S: Sync,
    {
        self.0.par_column_iter().map(|c| c.into()).into_par_iter()
    }
}

impl<T: Scalar + ClosedAdd + ClosedMul + Zero, R: Dim, C: Dim, S: RawStorage<T, R, C>>
    GenericMatrix<T, R, C, S>
{
    pub fn dot<R2: Dim, C2: Dim, S2: RawStorage<T, R2, C2>>(
        &self,
        rhs: &GenericMatrix<T, R2, C2, S2>,
    ) -> T
    where
        ShapeConstraint: DimEq<R, R2> + DimEq<C, C2>,
    {
        self.0.dot(&rhs.0)
    }
}

impl<T: Scalar + ClosedMul, R: Dim, C: Dim, S: Storage<T, R, C>> GenericMatrix<T, R, C, S>
where
    DefaultAllocator: Allocator<T, R, C>,
{
    pub fn component_mul(&self, rhs: &Self) -> GenericMatrix<T, R, C, Owned<T, R, C>> {
        self.0.component_mul(&rhs.0).into()
    }

    pub fn component_mul_assign(&mut self, rhs: &Self)
    where
        S: StorageMut<T, R, C>,
    {
        self.0.component_mul_assign(&rhs.0)
    }
}

impl<T: Scalar, R: Dim, C: Dim, S: RawStorage<T, R, C>> PartialEq for GenericMatrix<T, R, C, S>
where
    nalgebra::Matrix<T, R, C, S>: PartialEq,
{
    fn eq(&self, other: &Self) -> bool {
        self.0.eq(&other.0)
    }
}

impl<T: Scalar, R: Dim, C: Dim, S: RawStorage<T, R, C>> Eq for GenericMatrix<T, R, C, S> where
    nalgebra::Matrix<T, R, C, S>: Eq
{
}

impl<T: Scalar, R: Dim, C: Dim, S: RawStorage<T, R, C>> Hash for GenericMatrix<T, R, C, S>
where
    nalgebra::Matrix<T, R, C, S>: Hash,
{
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.0.hash(state)
    }
}

/// Implement unary operation `GenericMatrix<T>` -> `GenericMatrix<TO>`
macro_rules! impl_unop {
    ($op:ident, $OpTrait:ident) => {
        impl<
                T: Scalar,
                R: Dim,
                C: Dim,
                S: RawStorage<T, R, C>,
                TO: Scalar,
                RO: Dim,
                CO: Dim,
                SO: RawStorage<TO, RO, CO>,
            > $OpTrait for GenericMatrix<T, R, C, S>
        where
            nalgebra::Matrix<T, R, C, S>: $OpTrait<Output = nalgebra::Matrix<TO, RO, CO, SO>>,
        {
            type Output = GenericMatrix<TO, RO, CO, SO>;

            fn $op(self) -> Self::Output {
                self.0.$op().into()
            }
        }
    };
}

/// Implement binary operation `GenericMatrix<T>` x `GenericMatrix<TRhs>` -> `GenericMatrix<TO>`
macro_rules! impl_binop_matrix {
    ($op:ident, $OpTrait:ident) => {
        impl<
                T: Scalar,
                R: Dim,
                C: Dim,
                S: RawStorage<T, R, C>,
                TRhs: Scalar,
                RRhs: Dim,
                CRhs: Dim,
                SRhs: RawStorage<TRhs, RRhs, CRhs>,
                TO: Scalar,
                RO: Dim,
                CO: Dim,
                SO: RawStorage<TO, RO, CO>,
            > $OpTrait<GenericMatrix<TRhs, RRhs, CRhs, SRhs>> for GenericMatrix<T, R, C, S>
        where
            nalgebra::Matrix<T, R, C, S>: $OpTrait<
                nalgebra::Matrix<TRhs, RRhs, CRhs, SRhs>,
                Output = nalgebra::Matrix<TO, RO, CO, SO>,
            >,
        {
            type Output = GenericMatrix<TO, RO, CO, SO>;

            fn $op(self, rhs: GenericMatrix<TRhs, RRhs, CRhs, SRhs>) -> Self::Output {
                self.0.$op(rhs.0).into()
            }
        }

        impl<
                'a,
                'b,
                T: Scalar,
                R: Dim,
                C: Dim,
                S: RawStorage<T, R, C>,
                TRhs: Scalar,
                RRhs: Dim,
                CRhs: Dim,
                SRhs: RawStorage<TRhs, RRhs, CRhs>,
                TO: Scalar,
                RO: Dim,
                CO: Dim,
                SO: RawStorage<TO, RO, CO>,
            > $OpTrait<&'b GenericMatrix<TRhs, RRhs, CRhs, SRhs>> for &'a GenericMatrix<T, R, C, S>
        where
            &'a nalgebra::Matrix<T, R, C, S>: $OpTrait<
                &'b nalgebra::Matrix<TRhs, RRhs, CRhs, SRhs>,
                Output = nalgebra::Matrix<TO, RO, CO, SO>,
            >,
        {
            type Output = GenericMatrix<TO, RO, CO, SO>;

            fn $op(self, rhs: &'b GenericMatrix<TRhs, RRhs, CRhs, SRhs>) -> Self::Output {
                self.0.$op(&rhs.0).into()
            }
        }
    };
}

/// Implement binary operation `GenericMatrix<T>` x `TRhs` -> `GenericMatrix<TO>` and `&GenericMatrix<T>` x `TRhs` -> `GenericMatrix<TO>`
macro_rules! impl_binop {
    ($op:ident, $OpTrait:ident) => {
        impl<T: Scalar, R: Dim, C: Dim, S: RawStorage<T, R, C>, Rhs: Scalar, Output> $OpTrait<Rhs>
            for GenericMatrix<T, R, C, S>
        where
            nalgebra::Matrix<T, R, C, S>: $OpTrait<Rhs, Output = Output>,
        {
            type Output = Output;

            fn $op(self, rhs: Rhs) -> Self::Output {
                self.0.$op(rhs)
            }
        }

        impl<
                'a,
                T: Scalar,
                R: Dim,
                C: Dim,
                S: RawStorage<T, R, C>,
                Rhs: Scalar,
                O: Scalar,
                RO: Dim,
                CO: Dim,
                SO: RawStorage<O, RO, CO>,
            > $OpTrait<Rhs> for &'a GenericMatrix<T, R, C, S>
        where
            &'a nalgebra::Matrix<T, R, C, S>:
                $OpTrait<Rhs, Output = nalgebra::Matrix<O, RO, CO, SO>>,
        {
            type Output = GenericMatrix<O, RO, CO, SO>;

            fn $op(self, rhs: Rhs) -> Self::Output {
                self.0.$op(rhs).into()
            }
        }
    };
}

/// Implement binary assignment operation `GenericMatrix<T>` x `GenericMatrix<TRhs>` -> `GenericMatrix<T>`
macro_rules! impl_binop_assign_matrix {
    ($op:ident, $OpTrait:ident) => {
        impl<
                T: Scalar,
                R: Dim,
                C: Dim,
                S: RawStorage<T, R, C>,
                TRhs: Scalar,
                RRhs: Dim,
                CRhs: Dim,
                SRhs: RawStorage<TRhs, RRhs, CRhs>,
            > $OpTrait<GenericMatrix<TRhs, RRhs, CRhs, SRhs>> for GenericMatrix<T, R, C, S>
        where
            nalgebra::Matrix<T, R, C, S>: $OpTrait<nalgebra::Matrix<TRhs, RRhs, CRhs, SRhs>>,
        {
            fn $op(&mut self, rhs: GenericMatrix<TRhs, RRhs, CRhs, SRhs>) {
                self.0.$op(rhs.0)
            }
        }
    };
}

/// Implement binary assignment operation `GenericMatrix<T>` x `TRhs` -> `GenericMatrix<T>`
macro_rules! impl_binop_assign {
    ($op:ident, $OpTrait:ident) => {
        impl<T: Scalar, R: Dim, C: Dim, S: RawStorage<T, R, C>, Rhs> $OpTrait<Rhs>
            for GenericMatrix<T, R, C, S>
        where
            nalgebra::Matrix<T, R, C, S>: $OpTrait<Rhs>,
        {
            fn $op(&mut self, rhs: Rhs) {
                self.0.$op(rhs)
            }
        }
    };
}

impl_unop!(neg, Neg);
impl_binop_matrix!(add, Add);
impl_binop!(add, Add);
impl_binop_assign_matrix!(add_assign, AddAssign);
impl_binop_assign!(add_assign, AddAssign);
impl_binop_matrix!(sub, Sub);
impl_binop!(sub, Sub);
impl_binop_assign_matrix!(sub_assign, SubAssign);
impl_binop_assign!(sub_assign, SubAssign);
impl_binop_matrix!(mul, Mul);
impl_binop!(mul, Mul);
// No need for impl_binop_assign_matrix!(mul_assign, MulAssign), which is already covered by impl_binop_assign!(mul_assign, MulAssign) below
impl_binop_assign!(mul_assign, MulAssign);

impl<T: Scalar, R: Dim, C: Dim, S: RawStorage<T, R, C>> Sum for GenericMatrix<T, R, C, S>
where
    nalgebra::Matrix<T, R, C, S>: Sum,
{
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.map(|m| m.0).sum::<Self::Inner>().into()
    }
}

impl<T: Scalar, R: Dim, C: Dim, S: RawStorage<T, R, C>, SRhs: RawStorage<T, R, Const<1>>>
    Extend<GenericVector<T, R, SRhs>> for GenericMatrix<T, R, C, S>
where
    nalgebra::Matrix<T, R, C, S>: Extend<nalgebra::Vector<T, R, SRhs>>,
{
    fn extend<I: IntoIterator<Item = GenericVector<T, R, SRhs>>>(&mut self, iter: I) {
        self.0.extend(iter.into_iter().map(|v| v.0));
    }
}

impl<'a, T: Scalar, R: Dim, C: Dim, S: RawStorage<T, R, C>> IntoIterator
    for &'a GenericMatrix<T, R, C, S>
where
    &'a nalgebra::Matrix<T, R, C, S>: IntoIterator,
{
    type Item = <&'a nalgebra::Matrix<T, R, C, S> as IntoIterator>::Item;
    type IntoIter = <&'a nalgebra::Matrix<T, R, C, S> as IntoIterator>::IntoIter;
    fn into_iter(self) -> Self::IntoIter {
        self.0.into_iter()
    }
}

// impl<'a, T: Scalar, R: Dim, C: Dim, S: RawStorage<T, R, C>> IntoIterator
//     for &'a mut GenericMatrix<T, R, C, S>
// where
//     &'a mut nalgebra::Matrix<T, R, C, S>: IntoIterator,
// {
//     type Item = <&'a mut nalgebra::Matrix<T, R, C, S> as IntoIterator>::Item;
//     type IntoIter = <&'a mut nalgebra::Matrix<T, R, C, S> as IntoIterator>::IntoIter;
//     fn into_iter(self) -> Self::IntoIter {
//         self.0.into_iter()
//     }
// }

impl<T: Scalar + Zero + One + ClosedAdd + ClosedMul, R: Dim, C: Dim, S: Storage<T, R, C>>
    GenericMatrix<T, R, C, S>
{
    pub fn kronecker<R2: Dim, C2: Dim, S2>(
        &self,
        rhs: &GenericMatrix<T, R2, C2, S2>,
    ) -> GenericMatrix<T, DimProd<R, R2>, DimProd<C, C2>, Owned<T, DimProd<R, R2>, DimProd<C, C2>>>
    where
        R: DimMul<R2>,
        C: DimMul<C2>,
        S2: Storage<T, R2, C2>,
        DefaultAllocator: Allocator<T, <R as DimMul<R2>>::Output, <C as DimMul<C2>>::Output>,
    {
        self.0.kronecker(&rhs.0).into()
    }
}
