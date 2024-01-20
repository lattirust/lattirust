#![allow(non_snake_case)]
use ark_std::rand::thread_rng;

use crate::lattice_arithmetic::poly_ring::PolyRing;
use crate::lattice_arithmetic::ring::Ring;
use crate::lattice_arithmetic::traits::{IntegerDiv, Normed};

pub type Matrix<R> = nalgebra::DMatrix<R>;

// TODO: implement as Mul trait for Vector<R> so that left-multiplication with a scalar is possible
pub type Vector<R> = nalgebra::DVector<R>;


// TODO: pass rng as param
pub fn sample_uniform_mat<R: Ring>(m: usize, n: usize) -> Matrix<R> {
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

/// Returns a uniformly random vector r in Z_q^n such that ||r||_2 <= beta
pub fn sample_uniform_vec_bounded<R: PolyRing>(n: usize, beta: f64) -> Vector<R> where u64: From<<R as PolyRing>::BaseRing> {
    let beta_sq = R::BaseRing::from((beta * beta).floor() as u128);
    sample_uniform_vec(n).map(|x_i: R| R::from(
        x_i.coeffs().iter().map(
            |x_ij| x_ij.integer_div(beta_sq)
        ).collect::<Vec<R::BaseRing>>()
    ))
}

pub fn norm_sq_ringelem<R: Ring + Normed<u64>>(x: &R) -> u64 {
    x.norm_squared()
}

pub fn norm_sq_vec_basering<R: PolyRing>(v: &Vector<R::BaseRing>) -> u64
{
    let mut norm_sq: u128 = 0;
    for c in v.as_slice() {
        let c_int: i64 = <<R as PolyRing>::BaseRing as Into<i64>>::into(*c);
        norm_sq += (c_int * c_int) as u128;
    }
    norm_sq as u64
}

pub fn norm_vec_basering<R: PolyRing>(v: &Vector<R::BaseRing>) -> f64
{
    (norm_sq_vec_basering::<R>(v) as f64).sqrt()
}


pub fn norm_sq_vec<R: PolyRing>(v: &Vector<R>) -> u64
{
    norm_sq_vec_basering::<R>(&R::flattened(v))
}

pub fn norm_vec<R: PolyRing>(v: &Vector<R>) -> f64
{
    // TODO: check that this is indeed the norm as used in the paper
    (norm_sq_vec(v) as f64).sqrt()
}


#[allow(non_snake_case)]
#[cfg(test)]
mod tests {
    use crate::lattice_arithmetic::ring::Zq;

    use super::*;

    type Z3 = Zq<3>;

    #[test]
    fn test_sample_uniform() {
        let m = 10;
        let n = 20;
        let A = sample_uniform_mat::<Z3>(m, n);
        assert_eq!(A.nrows(), m);
        assert_eq!(A.ncols(), n);
    }
}