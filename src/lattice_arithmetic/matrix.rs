#![allow(non_snake_case)]
use ark_std::rand::thread_rng;

use crate::lattice_arithmetic::poly_ring::{PolyRing, SignedRepresentative, UnsignedRepresentative};
use crate::lattice_arithmetic::ring::Ring;
use crate::lattice_arithmetic::traits::{WithL2Norm};

pub type Matrix<R> = nalgebra::DMatrix<R>;

// TODO: implement as Mul trait for Vector<R> so that left-multiplication with a scalar is possible
pub type Vector<R> = nalgebra::DVector<R>;


// TODO: pass rng as param
pub fn sample_uniform_mat<R: PolyRing>(m: usize, n: usize) -> Matrix<R> {
    Matrix::<R>::from_fn(m, n, |_, _| R::rand(&mut thread_rng()))
}

pub fn sample_uniform_mat_symmetric<R: Ring>(m: usize, n: usize) -> Matrix<R> {
    let mut A = Matrix::<R>::zeros(m, n);
    for i in 0..m {
        for j in i..n {
            A[(i, j)] = R::rand(&mut thread_rng());
            A[(j, i)] = A[(i, j)]
        }
    }
    A
}

pub fn sample_uniform_vec<R: Ring>(n: usize) -> Vector<R> {
    Vector::<R>::from_fn(n, |_, _| R::rand(&mut thread_rng()))
}

#[allow(non_snake_case)]
#[cfg(test)]
mod tests {
    use crate::lattice_arithmetic::pow2_cyclotomic_poly_ring::Pow2CyclotomicPolyRing;
    use crate::lattice_arithmetic::ring::Fq;
    use super::*;

    type R = Pow2CyclotomicPolyRing<Fq<3>, 20>;

    #[test]
    fn test_sample_uniform() {
        let m = 10;
        let n = 20;
        let A = sample_uniform_mat::<R>(m, n);
        assert_eq!(A.nrows(), m);
        assert_eq!(A.ncols(), n);
    }
}