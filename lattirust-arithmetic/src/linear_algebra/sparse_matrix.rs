use std::error::Error;
use std::ops::{AddAssign, Mul};

use delegate::delegate;
use derive_more::{From, Index, IndexMut, Into, Mul, MulAssign};
use nalgebra::{Dim, Dyn, RawStorage};
use nalgebra_sparse;
use nalgebra_sparse::CooMatrix;
use num_traits::Zero;
use serde::{Deserialize, Serialize};

use crate::linear_algebra::generic_matrix::GenericMatrix;
use crate::linear_algebra::Scalar;
use crate::linear_algebra::{Matrix, Vector};

#[derive(Clone, Debug, PartialEq, From, Into, Mul, MulAssign, Index, IndexMut)]
pub struct SparseMatrix<R>(nalgebra_sparse::CscMatrix<R>); // We typically have more rows than columns, hence CSC.

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
    pub fn zeros(nrows: usize, ncols: usize) -> Self {
        nalgebra_sparse::CscMatrix::<R>::zeros(nrows, ncols).into()
    }
}

impl<R: Scalar + Copy + Zero + AddAssign> SparseMatrix<R> {
    pub fn try_from_triplets(
        nrows: usize,
        ncols: usize,
        triplets: Vec<(usize, usize, R)>,
    ) -> Result<Self, Box<dyn Error>> {
        let (row_index, col_index, value_index) = itertools::multiunzip(triplets);
        let coo =
            CooMatrix::<R>::try_from_triplets(nrows, ncols, row_index, col_index, value_index)?;
        Ok(SparseMatrix(nalgebra_sparse::CscMatrix::from(&coo)))
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

impl<R: Scalar + Zero> From<SparseMatrix<R>> for Matrix<R> {
    fn from(val: SparseMatrix<R>) -> Self {
        let mut dense = Matrix::<R>::zeros(val.nrows(), val.ncols());
        for (row, col, value) in val.0.triplet_iter() {
            dense[(row, col)] = value.clone();
        }
        dense
    }
}

impl<
        'a,
        'b,
        Lhs: Scalar,
        Rhs: Scalar,
        SRhs: RawStorage<Rhs, Dyn, Dyn>,
        Out: Scalar,
        ROut: Dim,
        COut: Dim,
        SOut: RawStorage<Out, ROut, COut>,
    > Mul<&'b GenericMatrix<Rhs, Dyn, Dyn, SRhs>> for &'a SparseMatrix<Lhs>
where
    &'a nalgebra_sparse::CscMatrix<Lhs>: Mul<
        &'b nalgebra::Matrix<Rhs, Dyn, Dyn, SRhs>,
        Output = nalgebra::Matrix<Out, ROut, COut, SOut>,
    >,
{
    type Output = GenericMatrix<Out, ROut, COut, SOut>;

    fn mul(self, rhs: &'b GenericMatrix<Rhs, Dyn, Dyn, SRhs>) -> Self::Output {
        self.0.mul(&rhs.0).into()
    }
}

/// SparseMatrix * Vector multiplication
impl<'b, Lhs: Scalar, Rhs: Scalar, Out: Scalar> Mul<&'b Vector<Rhs>> for &SparseMatrix<Lhs>
where
    Lhs: Mul<Rhs, Output = Out>,
    Out: Zero + AddAssign,
{
    type Output = Vector<Out>;

    fn mul(self, rhs: &'b Vector<Rhs>) -> Self::Output {
        self.0
            .triplet_iter()
            .fold(Vector::zeros(self.nrows()), |mut acc, (row, col, value)| {
                acc[row] += value.clone() * rhs[col].clone();
                acc
            })
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
        // nalgebra-sparse only supports sparse-dense multiplication
        let self_transpose = self.0.transpose();
        let rhs_transpose = rhs.0.transpose();
        let dense = rhs_transpose * self_transpose;
        dense.transpose().into()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const NUM_ROWS: usize = 1024;
    const NUM_COLS: usize = 2048;
    type R = u128;
    #[test]
    fn test_sparsematrix_vector_mul() {
        let triplets = Vec::from_iter((0..NUM_ROWS).map(|i| (i, i % NUM_COLS, R::from(i as u128))));
        let sparse_mat =
            SparseMatrix::<R>::try_from_triplets(NUM_ROWS, NUM_COLS, triplets).unwrap();
        let dense_mat = Matrix::from_fn(NUM_ROWS, NUM_COLS, |i, j| {
            if j == i % NUM_COLS {
                R::from(i as u128)
            } else {
                0
            }
        });
        let vec = Vector::from_fn(NUM_COLS, |i, _j| R::from((NUM_COLS - i) as u128));

        let res = &sparse_mat * &vec;
        let expected = &dense_mat * &vec;
        assert_eq!(res, expected);
    }
}
