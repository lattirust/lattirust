#![allow(non_snake_case)]

use std::cmp::max_by;

use serde::Serialize;

use lattice_estimator::msis::MSIS;
use lattice_estimator::norms::Norm;

use crate::labrador::binary_r1cs::util::SECPARAM;
use crate::labrador::prover::Witness;
use crate::lattice_arithmetic::challenge_set::labrador_challenge_set::LabradorChallengeSet;
use crate::lattice_arithmetic::matrix::{Matrix, norm_sq_vec, sample_uniform_mat, sample_uniform_vec, Vector};
use crate::lattice_arithmetic::poly_ring::PolyRing;
use crate::lattice_arithmetic::traits::WithLog2;
use crate::relations::labrador::principal_relation::{ConstantQuadDotProdFunction, PrincipalRelation, QuadDotProdFunction};

#[derive(Clone, Debug, Serialize)]
pub struct CommonReferenceString<R: PolyRing> {
    pub r: usize,
    pub n: usize,
    pub d: usize,
    pub norm_bound_squared: f64,
    pub k: usize,
    pub k1: usize,
    pub k2: usize,
    pub t1: usize,
    pub t2: usize,
    pub num_aggregs: usize,
    pub num_constraints: usize,
    pub num_constant_constraints: usize,
    pub A: Matrix<R>,
    // k x n
    pub B: Vec<Vec<Matrix<R>>>,
    // (r x t1) x k1 x k
    pub C: Vec<Vec<Vec<Vector<R>>>>,
    // (r x r x t2) x k2 x 1
    pub D: Vec<Vec<Vec<Vector<R>>>>,
    // (r x r x t1) x k2 x 1
    pub decomposition_basis: usize,
}

impl<R: PolyRing> CommonReferenceString<R> {
    fn t1(decomposition_basis: usize) -> usize {
        let log2_q: f64 = R::BaseRing::log2_q();
        let log2_b = (decomposition_basis as f64).log2();
        (log2_q / log2_b).round() as usize
    }

    fn t2(r: usize, n: usize, beta_sq: f64, decomposition_basis: usize) -> usize {
        let d = R::dimension();
        let log2_b = (decomposition_basis as f64).log2();
        let s_std_dev_sq: f64 = (beta_sq) / ((r * n * d) as f64); // standard deviation of s vectors = beta / sqrt(r * n * d)
        (f64::log2(f64::sqrt((24 * n * d) as f64) * s_std_dev_sq) / log2_b).round() as usize
    }

    pub fn new(r: usize, n: usize, beta_sq: f64, num_constraints: usize, num_constant_constraints: usize) -> CommonReferenceString<R> {
        let d = R::dimension();
        let q = R::modulus();
        let log2_q: f64 = R::BaseRing::log2_q();

        // Set decomposition basis
        let s = beta_sq.sqrt() / ((r * n * d) as f64).sqrt(); // standard deviation of the Z_q coefficients of the s vectors
        let tau = LabradorChallengeSet::<R>::VARIANCE_SUM_COEFFS;
        let b = ((s * (12. * r as f64 * tau).sqrt()).sqrt()).round() as usize;

        // Set t1 and t2
        let t1 = Self::t1(b);
        let t2 = Self::t2(r, n, beta_sq, b);
        let num_aggregs = (128. / log2_q).ceil() as usize;

        // Compute the norm bound for the next folded instance
        let beta_prime = |k| {
            Self::next_norm_bound_sq(r, n, beta_sq, k, b).sqrt()
        };

        let op_norm = LabradorChallengeSet::<R>::OPERATOR_NORM_THRESHOLD;
        let norm_bound_1 = |kappa| {
            // max(8T (b + 1)β′, 2(b + 1)β′ + 4T sqrt(128/30)β)
            max_by(8. * op_norm * (b + 1) as f64 * beta_prime(kappa),
                   2. * (b + 1) as f64 * beta_prime(kappa) + 4. * op_norm * f64::sqrt(128. / 30.) * beta_prime(kappa),
                   f64::total_cmp,
            )
        };

        // Ensure MSIS_{n=k, d, q, beta_1, m=n} is hard (l_2 norm)
        let msis_1 = MSIS {
            n: 0, // Dummy value, will be set later
            d,
            q,
            length_bound: 0., // Dummy value
            m: n,
            norm: Norm::L2,
        };
        let k = msis_1.find_optimal_n_dynamic(norm_bound_1, SECPARAM).expect("failed to find secure rank for {msis_1}");

        let msis_2 = MSIS {
            n: 0, // Dummy value, will be set later
            d,
            q,
            length_bound: 2. * beta_prime(k),
            m: k,
            norm: Norm::L2,
        };
        let k1 = msis_2.find_optimal_n(SECPARAM).expect("failed to find secure rank for {msis_2}");
        let k2 = k1;

        // TODO: this only gives 125 bits of soundness error for SECPARAM = 128, how do we best document this?
        // TODO: all params should be bigger to account for the slack of sqrt(128/30) per recursion level

        CommonReferenceString {
            r,
            n,
            d,
            norm_bound_squared: beta_sq,
            k,
            k1,
            k2,
            t1,
            t2,
            num_aggregs,
            num_constraints,
            num_constant_constraints,
            A: sample_uniform_mat(k, n),
            B: vec![vec![sample_uniform_mat(k1, k); t1]; r],
            C: (0..r).map(
                |i| (0..i + 1).map(
                    |_| vec![sample_uniform_vec(k2); t2]
                ).collect()
            ).collect(),
            D: (0..r).map(
                |i| (0..i + 1).map(
                    |_| vec![sample_uniform_vec(k2); t1]
                ).collect()
            ).collect(),
            decomposition_basis: b,
        }
    }

    /// Compute the squared norm bound for the next folded instance (cf. Section 5.4 of the Labrador paper)
    pub fn next_norm_bound_sq(r: usize, n: usize, norm_bound_squared: f64, k: usize, decomposition_basis: usize) -> f64 {
        let b_sq = decomposition_basis * decomposition_basis;
        let d = R::dimension();
        let challenge_variance = LabradorChallengeSet::<R>::VARIANCE_SUM_COEFFS;
        let t1 = Self::t1(decomposition_basis);
        let t2 = Self::t2(r, n, norm_bound_squared, decomposition_basis);

        let gamma_sq = norm_bound_squared * challenge_variance;
        let gamma_1_sq = (b_sq * t1) as f64 / 12. * (r * k * d) as f64 + (b_sq * t2) as f64 / 12. * ((r * (r + 1)).div_ceil(2) * d) as f64;
        let gamma_2_sq = (b_sq * t1) as f64 / 12. * ((r * (r + 1)).div_ceil(2) * d) as f64;
        let beta_next_sq: f64 = (2. / b_sq as f64) * gamma_sq + gamma_1_sq * gamma_2_sq;
        beta_next_sq
    }

    pub fn is_wellformed_constraint(&self, c: &QuadDotProdFunction<R>) -> bool {
        match c.A {
            Some(ref A) => A.nrows() == self.r && A.ncols() == self.r && A.transpose() == *A,
            None => true,
        }
    }

    pub fn is_wellformed_const_constraint(&self, c: &ConstantQuadDotProdFunction<R>) -> bool {
        match c.A {
            Some(ref A) => A.nrows() == self.r && A.ncols() == self.r && A.transpose() == *A,
            None => true,
        }
    }

    pub fn is_wellformed_instance(&self, instance: &PrincipalRelation<R>) -> bool {
        instance.quad_dot_prod_funcs.len() == self.num_constraints &&
            instance.quad_dot_prod_funcs.iter().all(|c| self.is_wellformed_constraint(c)) &&
            instance.ct_quad_dot_prod_funcs.len() == self.num_constant_constraints &&
            instance.ct_quad_dot_prod_funcs.iter().all(|c| self.is_wellformed_const_constraint(c))
    }

    pub fn is_wellformed_witness(&self, witness: &Witness<R>) -> bool {
        witness.s.len() == self.r &&
            witness.s.iter().all(|s_i| s_i.len() == self.n) &&
            witness.s.iter().map(|s_i| norm_sq_vec(s_i)).sum::<u64>() as f64 <= self.norm_bound_squared
    }
}