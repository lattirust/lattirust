use lattice_estimator::norms::Norm::L2;
use lattice_estimator::sis::SIS;

use crate::lattice_arithmetic::matrix::{Matrix, sample_uniform_mat};
use crate::lattice_arithmetic::poly_ring::ConvertibleField;

pub type SymmetricMatrix<F, const N: usize> = nalgebra::SMatrix<F, N, N>; // TODO: use lower triangular representation instead

pub const SECPARAM: usize = 128;

pub struct CRS<F: ConvertibleField> {
    pub commitment_mat: Matrix<F>,
    pub norm_bound: f64,
    pub decomposition_basis: u128,
    pub decomposition_length: usize,
    pub n: usize,
}

impl<F: ConvertibleField> CRS<F> {
    pub fn new(n: usize, norm_bound: f64) -> Self {
        let sis = SIS::new(0, F::modulus(), norm_bound, n, L2);
        let h = sis.find_optimal_n(SECPARAM).unwrap();
        let commitment_mat = sample_uniform_mat(h, n, &mut rand::thread_rng());
        let decomposition_basis = norm_bound.sqrt().round_ties_even() as u128; // TODO
        let decomposition_length = norm_bound.log(decomposition_basis as f64).ceil() as usize;
        Self { commitment_mat, norm_bound, decomposition_basis, decomposition_length, n }
    }
}