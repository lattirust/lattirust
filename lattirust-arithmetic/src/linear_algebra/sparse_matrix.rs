use std::ops::Mul;

use delegate::delegate;
use derive_more::{From, Index, IndexMut, Into, Mul, MulAssign};
use nalgebra::Dim;
use nalgebra_sparse;
use serde::{Deserialize, Serialize};

use crate::linear_algebra::generic_matrix::GenericMatrix;
use crate::linear_algebra::Matrix;
use crate::linear_algebra::Scalar;

#[derive(Clone, Debug, PartialEq, From, Into, Mul, MulAssign, Index, IndexMut)]
pub struct SparseMatrix<R>(nalgebra_sparse::CscMatrix<R>); // We typically have more rows than columns, hence CSC.

impl<R> SparseMatrix<R> {
    type Inner = nalgebra_sparse::CscMatrix<R>;
}

impl<R: Scalar> SparseMatrix<R> {
    delegate! {
        to self.0 {
            pub fn nrows(&self) -> usize;
            pub fn ncols(&self) -> usize;
            pub fn nnz(&self) -> usize;
            #[into]
            pub fn transpose(&self) -> Self;
        }
    }
}

impl<R> Serialize for SparseMatrix<R>
where
    nalgebra_sparse::CscMatrix<R>: Serialize,
{
    fn serialize<S: serde::Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        self.0.serialize(serializer)
    }
}

impl<'de, R> Deserialize<'de> for SparseMatrix<R>
where
    nalgebra_sparse::CscMatrix<R>: Deserialize<'de>,
{
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        Self::Inner::deserialize(deserializer).map(|x| x.into())
    }
}

impl<
        'a,
    'b,
        Lhs: Scalar,
        Rhs: Scalar,
        RRhs: Dim,
        CRhs: Dim,
        SRhs,
        Out: Scalar,
        ROut: Dim,
        COut: Dim,
        SOut,
    > Mul<&'b GenericMatrix<Rhs, RRhs, CRhs, SRhs>> for &'a SparseMatrix<Lhs>
where
    &'a nalgebra_sparse::CscMatrix<Lhs>: Mul<
        &'b
        nalgebra::Matrix<Rhs, RRhs, CRhs, SRhs>,
        Output = nalgebra::Matrix<Out, ROut, COut, SOut>,
    >,
{
    type Output = GenericMatrix<Out, ROut, COut, SOut>;

    fn mul(self, rhs: &'b GenericMatrix<Rhs, RRhs, CRhs, SRhs>) -> Self::Output {
        self.0.mul(&rhs.0).into()
    }
}

impl<'a, Lhs: Scalar, Rhs: Scalar, Out: Scalar, ROut: Dim, COut: Dim, SOut> Mul<Rhs>
    for &'a SparseMatrix<Lhs>
where
    &'a nalgebra_sparse::CscMatrix<Lhs>: Mul<Rhs, Output = nalgebra::Matrix<Out, ROut, COut, SOut>>,
{
    type Output = GenericMatrix<Out, ROut, COut, SOut>;

    fn mul(self, rhs: Rhs) -> Self::Output {
        self.0.mul(rhs).into()
    }
}

impl<'a, Lhs: Scalar, Rhs: Scalar, Out: Scalar> Mul<&'a SparseMatrix<Rhs>> for &'a Matrix<Lhs>
where
    nalgebra_sparse::CscMatrix<Rhs>: Mul<nalgebra::DMatrix<Lhs>, Output = nalgebra::DMatrix<Out>>,
{
    type Output = Matrix<Out>;

    fn mul(self, rhs: &'a SparseMatrix<Rhs>) -> Self::Output {
        // TODO: nalgebra-sparse only supports sparse-dense multiplication
        let self_transpose = self.0.transpose();
        let rhs_transpose = rhs.0.transpose();
        let dense = rhs_transpose * self_transpose;
        dense.transpose().into()
    }
}
