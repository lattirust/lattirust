use std::collections::VecDeque;

use num_traits::Zero;
use rayon::prelude::*;

use crate::linear_algebra::{Matrix, Scalar, SymmetricMatrix, Vector};
use crate::linear_algebra::generic_matrix::{ClosedAdd, ClosedMul};
use crate::poly_ring::PolyRing;

/// Convert the entries of a lower triangular n x n matrix (in sparse representation) to a vector of length (n*(n+1)) / 2
#[inline(always)]
pub fn vec_from_lowertriang<T>(mut m: VecDeque<VecDeque<T>>) -> Vec<T> {
    debug_assert!(m.len() > 0);
    let mut v = Vec::<T>::with_capacity((m.len() * (m.len() + 1)) / 2);
    for i in 0..m.len() {
        let mut m_i = m.pop_front().unwrap();
        debug_assert_eq!(
            m_i.len(),
            i + 1,
            "representation of lower triangular matrix has wrong dimensions"
        );
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
    (0..n)
        .map(|i| (0..i + 1).map(|_| v.pop_front().unwrap()).collect())
        .collect()
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

pub fn inner_products_serial<R: PolyRing>(s: &Vec<Vector<R>>) -> SymmetricMatrix<R> {
    let mut symmetric_matrix = vec![vec![]; s.len()];
    for i in 0..s.len() {
        symmetric_matrix[i] = Vec::<R>::with_capacity(i + 1);
        for j in 0..i + 1 {
            symmetric_matrix[i].push(s[i].dot(&s[j]));
        }
    }
    symmetric_matrix.into()
}

/// Compute $(\langle s_{:,i}, s_{:,j}\rangle)_{i, j \in \[n\]}$ for $s \in R^{n \times m}$
pub fn inner_products<R: PolyRing>(s: &Vec<Vector<R>>) -> SymmetricMatrix<R> {
    inner_products2(s, s)
}

/// Compute $(\langle s_{:,i}, s_{:,j}\rangle)_{i, j \in \[n\]}$, where $s \in R^{m \times, n}$
/// This is equivalent to the lower triangular part of the symmetric matrix $s^T \cdot s$.
pub fn inner_products_mat<R: Scalar + ClosedAdd + ClosedMul + Zero + Sync + Send>(
    s: &Matrix<R>,
) -> SymmetricMatrix<R> {
    let ranges = lower_triang_indices(s.ncols());

    lowertriang_from_vec(
        ranges
            .into_par_iter()
            .map(|(i, j)| s.column(i).dot(&s.column(j)))
            .collect::<VecDeque<_>>(),
        s.ncols(),
    )
    .into()
}

/// Compute $(\langle s_i, t_j\rangle)_{i, j \in \[n\]}$ for $s,t \in R^{n \times m}$
pub fn inner_products2<R: PolyRing>(s: &Vec<Vector<R>>, t: &Vec<Vector<R>>) -> SymmetricMatrix<R> {
    debug_assert_eq!(s.len(), t.len());
    let ranges = lower_triang_indices(s.len());

    lowertriang_from_vec(
        ranges
            .into_par_iter()
            .map(|(i, j)| s[i].dot(&t[j]))
            .collect::<VecDeque<_>>(),
        s.len(),
    )
    .into()
}

#[cfg(test)]
mod tests {
    use ark_std::test_rng;

    use crate::linear_algebra::symmetric_matrix::SymmetricMatrix;
    use crate::ntt::ntt_modulus;
    use crate::pow2_cyclotomic_poly_ring_ntt::Pow2CyclotomicPolyRingNTT;
    use crate::ring::Fq;

    use super::*;

    const Q: u64 = ntt_modulus::<64>(32);

    type PR = Pow2CyclotomicPolyRingNTT<Q, 64>;

    #[test]
    fn test_lowertriang_vec() {
        let n = 100;
        let dim = (n * (n + 1)) / 2;
        let x = (0..dim).collect::<VecDeque<_>>();
        let mat = lowertriang_from_vec(x.clone().into(), n);
        let mat_ = mat
            .clone()
            .into_iter()
            .map(|x| VecDeque::from(x))
            .collect::<VecDeque<_>>();

        assert_eq!(mat_.len(), n);
        for i in 0..mat_.len() {
            assert_eq!(mat_[i].len(), i + 1);
        }
        assert_eq!(x, vec_from_lowertriang(mat_.clone()));

        assert_eq!(
            mat,
            lowertriang_from_vec(vec_from_lowertriang(mat_.clone()).into(), n)
        );
    }

    #[test]
    fn test_inner_products() {
        // Test parallelized implementation against a straightforward serial implementation
        let v = vec![Vector::<PR>::rand(2, &mut test_rng()); 3];
        assert_eq!(inner_products_serial(&v), inner_products(&v));
    }

    #[test]
    fn test_inner_products_mat() {
        let rng = &mut test_rng();
        let mat = Matrix::<Fq<Q>>::rand(101, 42, rng);
        let inner_prods = inner_products_mat(&mat);
        let inner_prods_expect: SymmetricMatrix<Fq<Q>> = (mat.transpose() * mat).into();
        assert_eq!(inner_prods, inner_prods_expect);
    }
}
