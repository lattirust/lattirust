#![allow(non_snake_case)]

use std::collections::VecDeque;
use ark_ff::Field;

use nalgebra::Scalar;
use rayon::iter::IntoParallelIterator;
use rayon::iter::ParallelIterator;

use crate::lattice_arithmetic::matrix::{Matrix, Vector};
use crate::lattice_arithmetic::poly_ring::PolyRing;
use crate::lattice_arithmetic::ring::Ring;

pub fn commit<R: Ring>(A: &Matrix<R>, s: &Vector<R>) -> Vector<R> {
    //(A * s.into()).into::<Matrix<R>>()
    A * s
}

/// Split a vector of size n into ceil(n/new_len) vectors of size new_len
pub fn chunk_pad<T: Ring>(v: &Vector<T>, chunk_size: usize) -> Vec<Vector<T>> {
    let chunk_iter = v.data.as_slice().chunks_exact(chunk_size);
    let remainder = chunk_iter.remainder();
    let mut res = chunk_iter.map(|x| Vector::<T>::from_row_slice(x)).collect::<Vec<_>>();
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
    let mut vals = v.data.as_slice().to_vec();
    vals.resize(new_len * m, T::zero());
    let vs = vals.as_slice();
    for i in 0..m {
        res.push(Vector::<T>::from_row_slice(&vs[i * new_len..(i + 1) * new_len]));
    }
    res
}

pub fn flatten_vec_vector<R: Ring>(v: &Vec<Vector<R>>) -> Vector<R> {
    let mut res = Vec::<R>::with_capacity(v.len() * v[0].len());
    for i in 0..v.len() {
        res.extend(v[i].data.as_vec());
    }
    Vector::<R>::from_vec(res)
}

pub fn flatten_symmetric_matrix<R: Ring>(v: &Vec<Vec<R>>) -> Vector<R> {
    Vector::<R>::from_vec(v.into_iter().flatten().cloned().collect())
}

pub fn concat<R: Clone + Scalar>(vecs: &[&Vector<R>]) -> Vector<R> {
    let vals = vecs.into_iter().map(|v| v.data.as_vec()).cloned().flatten().collect::<Vec<R>>();
    Vector::<R>::from_vec(vals)
}

pub fn shift_right<R: Ring>(v: &Vec<R>, shift: usize) -> Vec<R> {
    let mut res = vec![R::zero(); shift];
    res.extend(v);
    res
}

pub fn inner_products_serial<R: PolyRing>(s: &Vec<Vector<R>>) -> Vec<Vec<R>> {
    let mut G = vec![vec![]; s.len()];
    for i in 0..s.len() {
        G[i] = Vec::<R>::with_capacity(i + 1);
        for j in 0..i + 1 {
            G[i].push(&s[i].dot(&s[j]));
        }
    }
    G
}

pub fn mul_matrix_basescalar<R: PolyRing>(A: &Matrix<R>, x: R::BaseRing) -> Matrix<R> {
    A.map(|a_ij| a_ij * x)
}

pub fn mul_basescalar_vector<R: PolyRing>(s: R::BaseRing, A: &Vector<R>) -> Vector<R> {
    A.map(|a_ij| a_ij * s)
}

/// Convert the entries of a lower triangular n x n matrix (in sparse representation) to a vector of length (n*(n+1)) / 2
#[inline(always)]
pub fn vec_from_lowertriang<T>(mut m: VecDeque<VecDeque<T>>) -> Vec<T> {
    debug_assert!(m.len() > 0);
    let mut v = Vec::<T>::with_capacity((m.len() * (m.len() + 1)) / 2);
    for i in 0..m.len() {
        let mut m_i = m.pop_front().unwrap();
        debug_assert_eq!(m_i.len(), i + 1, "representation of lower triangular matrix has wrong dimensions");
        for _ in 0..i + 1 {
            v.push(m_i.pop_front().unwrap()); // repeatedly remove
        }
    }
    v
}

/// Convert a vector of length (n*(n+1)) / 2 to the sparse representation of a lower triangular n x n matrix
#[inline(always)]
pub fn lowertriang_from_vec<T>(mut v: VecDeque<T>, n: usize) -> Vec<Vec<T>> {
    debug_assert_eq!(v.len(), n * (n + 1) / 2);
    (0..n).map(|i|
        (0..i + 1).map(
            |_| v.pop_front().unwrap()
        ).collect()
    ).collect()
}

#[inline(always)]
pub fn lower_triang_indices(n: usize) -> Vec<(usize, usize)> {
    let mut indices = Vec::<(usize, usize)>::with_capacity((n * (n + 1)) / 2);
    for i in 0..n {
        for j in 0..i + 1 {
            indices.push((i, j));
        }
    }
    indices
}

pub fn inner_products<R: PolyRing>(s: &Vec<Vector<R>>) -> Vec<Vec<R>> {
    inner_products2(s, s)
}

/// Compute (<s[i], t[j])_ij, for 0 <= i < n, 0 <= j < n
pub fn inner_products2<R: PolyRing>(s: &Vec<Vector<R>>, t: &Vec<Vector<R>>) -> Vec<Vec<R>> {
    debug_assert_eq!(s.len(), t.len());
    let ranges = lower_triang_indices(s.len());

    lowertriang_from_vec(
        ranges.into_par_iter().map(
            |(i, j)| &s[i].dot(&t[j])
        ).collect::<VecDeque<_>>(),
        s.len(),
    )
}

// Compute sum_{i,j in [r]} A_ij * c_i * c_j
pub fn linear_combination_symmetric_matrix<R: Ring>(A: &Vec<Vec<R>>, c: &Vec<R>) -> R {
    let n = A.len();
    debug_assert_eq!(c.len(), n);
    let mut lc = R::zero();
    for i in 0..n {
        debug_assert_eq!(A[i].len(), i + 1);
        for j in 0..i + 1 {
            lc += A[i][j] * c[i] * c[j];
        }
        for j in i + 1..n {
            // for j >= i+1, get A_ij from A_ji, since A is stored in lower triangular representation
            lc += A[j][i] * c[i] * c[j];
        }
    }
    lc
}

#[cfg(test)]
mod tests {
    use crate::lattice_arithmetic::matrix::sample_uniform_vec;
    use crate::lattice_arithmetic::ntt::ntt_modulus;
    use crate::lattice_arithmetic::pow2_cyclotomic_poly_ring_ntt::Pow2CyclotomicPolyRingNTT;

    use super::*;

    const Q: u64 = ntt_modulus::<64>(32);

    type PR = Pow2CyclotomicPolyRingNTT<Q, 64>;

    #[test]
    fn test_lowertriang_vec() {
        let n = 100;
        let dim = (n * (n + 1)) / 2;
        let x = (0..dim).collect::<VecDeque<_>>();
        let mat = lowertriang_from_vec(x.clone().into(), n);
        let mat_ = mat.clone().into_iter().map(|x| VecDeque::from(x))
            .collect::<VecDeque<_>>();

        assert_eq!(mat_.len(), n);
        for i in 0..mat_.len() {
            assert_eq!(mat_[i].len(), i + 1);
        }
        assert_eq!(x, vec_from_lowertriang(mat_.clone()));


        assert_eq!(mat, lowertriang_from_vec(vec_from_lowertriang(mat_.clone()).into(), n));
    }

    #[test]
    fn test_inner_products() {
        // Test parallelized implementation against a straightforward serial implementation
        let v = vec![sample_uniform_vec::<PR>(2); 3];
        assert_eq!(inner_products_serial(&v), inner_products(&v));
    }
}