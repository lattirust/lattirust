use std::ops::{Add, Mul, Sub};

pub mod generic_matrix;
pub mod inner_products;
mod matrix;
pub mod serialization;
mod sparse_matrix;
mod symmetric_matrix;
mod vector;

pub type Vector<T> = vector::Vector<T>;
pub type SVector<T, const N: usize> = vector::SVector<T, N>;
pub type RowVector<T> = vector::RowVector<T>;
pub type SRowVector<T, const N: usize> = vector::SRowVector<T, N>;
pub type Matrix<T> = matrix::Matrix<T>;
pub type SMatrix<T, const R: usize, const C: usize> = matrix::SMatrix<T, R, C>;
pub type SymmetricMatrix<T> = symmetric_matrix::SymmetricMatrix<T>;
pub type SparseMatrix<T> = sparse_matrix::SparseMatrix<T>;
pub trait Scalar = nalgebra::Scalar;

pub trait ClosedAdd: Add<Self, Output = Self> + Sized {}
impl<T> ClosedAdd for T where T: Add<T, Output = T> {}
pub trait ClosedAddAssign: nalgebra::ClosedAddAssign {}
impl<T> ClosedAddAssign for T where T: nalgebra::ClosedAddAssign {}

#[allow(dead_code)]
pub trait ClosedSub: Sub<Self, Output = Self> + Sized {}
impl<T> ClosedSub for T where T: Sub<T, Output = T> {}
pub trait ClosedSubAssign: nalgebra::ClosedSubAssign {}
impl<T> ClosedSubAssign for T where T: nalgebra::ClosedSubAssign {}

pub trait ClosedMul: Mul<Self, Output = Self> + Sized {}
impl<T> ClosedMul for T where T: Mul<Self, Output = Self> {}
pub trait ClosedMulAssign: nalgebra::ClosedMulAssign {}
impl<T> ClosedMulAssign for T where T: nalgebra::ClosedMulAssign {}
