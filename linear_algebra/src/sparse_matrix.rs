use ark_std::{
    ops::{AddAssign, Mul},
    Zero,
};
use delegate::delegate;
use derive_more::{From, Index, IndexMut, Into, Mul, MulAssign};
use nalgebra::{Dim, RawStorage};
use nalgebra_sparse::{
    csc::{CscCol, CscTripletIter},
    CooMatrix, CscMatrix, SparseFormatError,
};
use serde::{Deserialize, Serialize};

use crate::{generic_matrix::GenericMatrix, Matrix, Scalar};

#[derive(Clone, Debug, PartialEq, From, Into, Mul, MulAssign, Index, IndexMut)]
pub struct SparseMatrix<R>(nalgebra_sparse::CscMatrix<R>); // We typically have more rows than columns, hence CSC.

impl<R: Scalar> SparseMatrix<R> {
    delegate! {
        to self.0 {
            pub fn nrows(&self) -> usize;
            pub fn ncols(&self) -> usize;
            pub fn nnz(&self) -> usize;
            pub fn triplet_iter(&self) -> CscTripletIter<'_, R>;
            pub fn col_offsets(&self) -> &[usize];
            pub fn row_indices(&self) -> &[usize];
            pub fn values(&self) -> &[R];
            pub fn get_col(&self, index: usize) -> Option<CscCol<'_, R>>;
            pub fn col(&self, index: usize) -> CscCol<'_, R>;
            #[into]
            pub fn transpose(&self) -> Self;
        }
    }
}

impl<'a, R: Scalar + AddAssign<R> + Zero> From<&'a [Vec<R>]> for SparseMatrix<R> {
    fn from(matrix: &'a [Vec<R>]) -> Self {
        let mut coo_matrix: CooMatrix<R> = CooMatrix::<R>::new(matrix.len(), matrix[0].len());

        for (i, row) in matrix.iter().enumerate() {
            for (j, value) in row.iter().enumerate() {
                coo_matrix.push(i, j, value.clone());
            }
        }

        SparseMatrix((&coo_matrix).into())
    }
}

impl<R: Scalar> SparseMatrix<R> {
    pub fn try_from_csc_data(
        num_rows: usize,
        num_cols: usize,
        col_offsets: Vec<usize>,
        row_indices: Vec<usize>,
        values: Vec<R>,
    ) -> Result<Self, SparseFormatError> {
        CscMatrix::try_from_csc_data(num_rows, num_cols, col_offsets, row_indices, values)
            .map(|matrix| SparseMatrix(matrix))
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
        nalgebra_sparse::CscMatrix::<R>::deserialize(deserializer).map(|x| x.into())
    }
}

impl<
        'a,
        'b,
        Lhs: Scalar,
        Rhs: Scalar,
        RRhs: Dim,
        CRhs: Dim,
        SRhs: RawStorage<Rhs, RRhs, CRhs>,
        Out: Scalar,
        ROut: Dim,
        COut: Dim,
        SOut: RawStorage<Out, ROut, COut>,
    > Mul<&'b GenericMatrix<Rhs, RRhs, CRhs, SRhs>> for &'a SparseMatrix<Lhs>
where
    &'a nalgebra_sparse::CscMatrix<Lhs>: Mul<
        &'b nalgebra::Matrix<Rhs, RRhs, CRhs, SRhs>,
        Output = nalgebra::Matrix<Out, ROut, COut, SOut>,
    >,
{
    type Output = GenericMatrix<Out, ROut, COut, SOut>;

    fn mul(self, rhs: &'b GenericMatrix<Rhs, RRhs, CRhs, SRhs>) -> Self::Output {
        self.0.mul(&rhs.0).into()
    }
}

impl<
        'a,
        Lhs: Scalar,
        Rhs: Scalar,
        Out: Scalar,
        ROut: Dim,
        COut: Dim,
        SOut: RawStorage<Out, ROut, COut>,
    > Mul<Rhs> for &'a SparseMatrix<Lhs>
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
