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
    pub fn from_ark_matrix(
        matrix: ark_relations::r1cs::Matrix<R>,
        nrows: usize,
        ncols: usize,
    ) -> Self {
        assert_eq!(nrows, matrix.len());
        let triplets = matrix
            .iter()
            .enumerate()
            .map(|(row_index, row)| {
                row.iter()
                    .map(move |(elem, col_index)| (row_index, *col_index, *elem))
            })
            .flatten()
            .collect();
        let res = Self::try_from_triplets(nrows, ncols, triplets).unwrap();
        res
    }

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
impl<'a, 'b, Lhs: Scalar, Rhs: Scalar, Out: Scalar> Mul<&'b Vector<Rhs>> for &'a SparseMatrix<Lhs>
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
