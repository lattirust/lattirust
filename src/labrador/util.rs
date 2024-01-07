use std::collections::VecDeque;
use rayon::iter::IntoParallelIterator;
use rayon::iter::ParallelIterator;

use crate::lattice_arithmetic::matrix::{Matrix, Vector};
use crate::lattice_arithmetic::ring::Ring;

pub fn commit<R: Ring>(A: &Matrix<R>, s: &Vector<R>) -> Vector<R> {
    //(A * s.into()).into::<Matrix<R>>()
    A * s
}

pub fn inner_prod<R: Ring>(v: &Vector<R>, w: &Vector<R>) -> R {
    v.dot(w)
}

pub fn inner_products_serial<R: Ring>(s: &Vec<Vector<R>>) -> Vec<Vec<R>> {
    let mut G = vec![vec![]; s.len()];
    for i in 0..s.len() {
        G[i] = Vec::<R>::with_capacity(i + 1);
        for j in 0..i + 1 {
            G[i].push(inner_prod(&s[i], &s[j]));
        }
    }
    G
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
    debug_assert_eq!(v.len(), (n * (n + 1) / 2));
    (0..n).map(|i|
        (0..i + 1).map(
            |j| v.pop_front().unwrap()
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


pub fn inner_products<R: Ring>(s: &Vec<Vector<R>>) -> Vec<Vec<R>> {
    let mut ranges = lower_triang_indices(s.len());

    lowertriang_from_vec(
        ranges.into_par_iter().map(
            |(i, j)| inner_prod(&s[i], &s[j])
        ).collect::<VecDeque<_>>(),
        s.len(),
    )
}

mod tests {
    use super::*;

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
}