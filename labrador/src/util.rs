#![allow(non_snake_case)]

use ark_relations::r1cs::ConstraintSystem;

use lattirust_arithmetic::linear_algebra::{Matrix, Scalar, SparseMatrix, SymmetricMatrix, Vector};
use lattirust_arithmetic::poly_ring::PolyRing;
use lattirust_arithmetic::ring::Ring;

use crate::binary_r1cs::util::Z2;

pub fn commit<R: Ring>(A: &Matrix<R>, s: &Vector<R>) -> Vector<R> {
    A * s
}

/// Split a vector of size n into ceil(n/new_len) vectors of size new_len
pub fn chunk_pad<T: Ring>(v: &Vector<T>, chunk_size: usize) -> Vec<Vector<T>> {
    let chunk_iter = v.as_slice().chunks_exact(chunk_size);
    let remainder = chunk_iter.remainder();
    let mut res = chunk_iter
        .map(|x| Vector::<T>::from_slice(x))
        .collect::<Vec<_>>();
    // Pad last entry to chunk_size
    let mut last = vec![T::zero(); chunk_size];
    last[..remainder.len()].copy_from_slice(remainder);
    res.push(Vector::<T>::from(last));
    res
}

/// Split a vector of size n into m vectors of size ceil(n/m)
pub fn split<T: Ring>(v: &Vector<T>, m: usize) -> Vec<Vector<T>> {
    let mut res = Vec::<Vector<T>>::with_capacity(m);
    let new_len = v.len().div_ceil(m);
    let mut vals = v.as_slice().to_vec();
    vals.resize(new_len * m, T::zero());
    let vs = vals.as_slice();
    for i in 0..m {
        res.push(Vector::<T>::from_slice(&vs[i * new_len..(i + 1) * new_len]));
    }
    res
}

pub fn flatten_vec_vector<R: Ring>(v: &Vec<Vector<R>>) -> Vector<R> {
    let mut res = Vec::<R>::with_capacity(v.len() * v[0].len());
    for i in 0..v.len() {
        res.extend(v[i].as_slice());
    }
    Vector::<R>::from_vec(res)
}

pub fn flatten_symmetric_matrix<R: Ring>(v: &SymmetricMatrix<R>) -> Vector<R> {
    Vector::<R>::from_vec(v.rows().into_iter().cloned().flatten().collect())
}

pub fn concat<R: Clone + Scalar>(vecs: &[&[R]]) -> Vector<R> {
    let mut vals = Vec::<R>::with_capacity(vecs.into_iter().map(|v| v.len()).sum());
    for v in vecs {
        vals.extend_from_slice(v);
    }
    Vector::<R>::from_vec(vals)
}

pub fn shift_right<R: Ring>(v: &Vec<R>, shift: usize) -> Vec<R> {
    let mut res = vec![R::zero(); shift];
    res.extend(v);
    res
}

pub fn mul_matrix_basescalar<R: PolyRing>(A: &Matrix<R>, x: R::BaseRing) -> Matrix<R> {
    A.map(|a_ij| a_ij * x)
}

pub fn mul_basescalar_vector<R: PolyRing>(s: R::BaseRing, A: &Vector<R>) -> Vector<R> {
    A.map(|a_ij| a_ij * s)
}

/// Compute $\sum_{i,j \in \[r\]} A_ij c_i  c_j$
pub fn linear_combination_symmetric_matrix<R: Ring>(A: &SymmetricMatrix<R>, c: &Vec<R>) -> R {
    let n = A.size();
    debug_assert_eq!(c.len(), n);
    let mut lc = R::zero();
    for i in 0..n {
        for j in 0..n {
            lc += A[(i, j)] * c[i] * c[j];
        }
    }
    lc
}

pub fn ark_sparse_matrices(
    cs: &ConstraintSystem<Z2>,
) -> (SparseMatrix<Z2>, SparseMatrix<Z2>, SparseMatrix<Z2>) {
    todo!()
}
