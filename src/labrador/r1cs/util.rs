#![allow(non_snake_case)]

use crate::labrador::binary_r1cs::util::SECPARAM;
use crate::labrador::common_reference_string::CommonReferenceString;
use crate::lattice_arithmetic::matrix::{Matrix, sample_uniform_mat};
use crate::lattice_arithmetic::poly_ring::PolyRing;
use crate::lattice_arithmetic::ring::Ring;
use crate::lattice_arithmetic::traits::Modulus;

pub struct R1CSCRS<R: PolyRing> {
    pub A: Matrix<R>,
    pub B: Matrix<R>,
    pub num_constraints: usize,
    pub num_variables: usize,
    pub m: usize,
    pub m_d: usize,
    pub l: usize,
}

impl<R: PolyRing> R1CSCRS<R> {
    pub fn new(num_constraints: usize, num_variables: usize) -> Self {
        let k = num_constraints;
        let n = num_variables;
        let d = R::dimension();
        // TODO: pad instead
        assert_eq!(k % d, 0, "number of constraints {} must be multiple of d = {}", k, d);
        assert_eq!(n % d, 0, "number of variables {} must be multiple of d = {}", n, d);
        let m = 10usize; // TODO: choose such that MSIS_{m, 2n+6k} is hard with l_inf bound = 1

        let q = R::modulus() as usize;
        assert!(n + 3 * k < q, "n + 3k = {} must be less than q = {} for soundness", n + 3 * k, q);
        assert!(6 * k < q, "6k = {} must be less than q = {} for soundness", 6 * k, q);
        assert!(128 * (n + 3 * k) < 15 * q, "n + 3*k = {} must be less than 15q/128 = {} to be able to compose with Labrador-core", n + 3 * k, 15 * q / 128);

        let m_d = 10usize; // TODO: choose such that MSIS_{m_d, l*k} is hard with l_2 norm bound = beta
        let l = (SECPARAM - 1).div_ceil(18); // protocol has soundness 2*p^-l, where p ~= 2^18 is the smallest prime factor of 2^64 + 1

        Self {
            A: sample_uniform_mat(m.div_ceil(d), (3 * k + n).div_ceil(d)),
            B: sample_uniform_mat(m.div_ceil(d), l * k),
            num_constraints,
            num_variables,
            m,
            m_d,
            l,
        }
    }

    pub fn pr_crs(&self) -> CommonReferenceString<R> {
        let d = R::dimension();
        let r_pr: usize = 8;
        let n_pr = self.num_variables.div_ceil(d);
        let norm_bound = (R::modulus() as f64).sqrt();

        let num_quad_constraints = self.m.div_ceil(d) + 3 * n_pr;
        let num_constant_quad_constraints = 4 + 1 + SECPARAM;

        CommonReferenceString::<R>::new(r_pr, n_pr, norm_bound, num_quad_constraints, num_constant_quad_constraints)
    }
}

pub struct R1CSInstance<R> {
    // TODO: use sparse matrices instead
    pub A: Matrix<R>,
    pub B: Matrix<R>,
    pub C: Matrix<R>,
}