#![allow(non_snake_case)]

use ark_std::rand::thread_rng;
use ark_std::UniformRand;
use nalgebra::Scalar;
use rand::Rng;

use crate::lattice_arithmetic::poly_ring::PolyRing;
use crate::lattice_arithmetic::ring::Ring;
use crate::lattice_arithmetic::traits::WithL2Norm;

pub type Matrix<R> = nalgebra::DMatrix<R>;

pub type SparseMatrix<R> = nalgebra_sparse::CscMatrix<R>; // We typically have more rows than columns, hence CSC.

// TODO: implement as Mul trait for Vector<R> so that left-multiplication with a scalar is possible
pub type Vector<R> = nalgebra::DVector<R>;

pub fn sample_uniform_mat<R: Scalar + UniformRand, Rng: rand::Rng + ?Sized>(m: usize, n: usize, rng: &mut Rng) -> Matrix<R> {
    Matrix::<R>::from_fn(m, n, |_, _| R::rand(rng))
}

pub fn sample_uniform_mat_symmetric<R: Ring, Rng: rand::Rng + ?Sized>(m: usize, n: usize, rng: &mut Rng) -> Matrix<R> {
    let mut A = Matrix::<R>::zeros(m, n);
    for i in 0..m {
        for j in i..n {
            A[(i, j)] = R::rand(&mut thread_rng());
            A[(j, i)] = A[(i, j)]
        }
    }
    A
}

pub fn sample_uniform_vec<R: Ring, Rng: rand::Rng + ?Sized>(n: usize, rng: &mut Rng) -> Vector<R> {
    Vector::<R>::from_fn(n, |_, _| R::rand(rng))
}

#[allow(non_snake_case)]
#[cfg(test)]
mod tests {
    use ark_std::test_rng;

    use crate::lattice_arithmetic::pow2_cyclotomic_poly_ring::Pow2CyclotomicPolyRing;
    use crate::lattice_arithmetic::ring::Fq;

    use super::*;

    type R = Pow2CyclotomicPolyRing<Fq<3>, 20>;

    #[test]
    fn test_sample_uniform() {
        let m = 10;
        let n = 20;
        let rng = &mut test_rng();
        let A = sample_uniform_mat::<R, _>(m, n, rng);
        assert_eq!(A.nrows(), m);
        assert_eq!(A.ncols(), n);
    }
}