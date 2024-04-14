#![allow(non_snake_case)]

use ark_std::rand::thread_rng;
use num_traits::{One, Zero};

use lattirust_arithmetic::linear_algebra::{Matrix, Scalar};
use lattirust_arithmetic::poly_ring::PolyRing;
use lattirust_arithmetic::ring::Ring;
use lattirust_arithmetic::linear_algebra::Vector;

use crate::binary_r1cs::util::Z2;
use crate::common_reference_string::{CommonReferenceString, SECURITY_PARAMETER};

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
        assert_eq!(
            k % d,
            0,
            "number of constraints {} must be multiple of d = {}",
            k,
            d
        );
        assert_eq!(
            n % d,
            0,
            "number of variables {} must be multiple of d = {}",
            n,
            d
        );
        let m = 10usize; // TODO: choose such that MSIS_{m, 2n+6k} is hard with l_inf bound = 1

        let q = R::modulus() as usize;
        assert!(
            n + 3 * k < q,
            "n + 3k = {} must be less than q = {} for soundness",
            n + 3 * k,
            q
        );
        assert!(
            6 * k < q,
            "6k = {} must be less than q = {} for soundness",
            6 * k,
            q
        );
        assert!(
            128 * (n + 3 * k) < 15 * q,
            "n + 3*k = {} must be less than 15q/128 = {} to be able to compose with Labrador-core",
            n + 3 * k,
            15 * q / 128
        );

        let m_d = 10usize; // TODO: choose such that MSIS_{m_d, l*k} is hard with l_2 norm bound = beta
        let l = (SECURITY_PARAMETER - 1).div_ceil(18); // protocol has soundness 2*p^-l, where p ~= 2^18 is the smallest prime factor of 2^64 + 1

        let rng = &mut ark_std::rand::rngs::OsRng::default();
        Self {
            A: Matrix::<R>::rand(m.div_ceil(d), (3 * k + n).div_ceil(d), rng),
            B: Matrix::<R>::rand(m.div_ceil(d), l * k, rng),
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
        let num_constant_quad_constraints = 4 + 1 + SECURITY_PARAMETER;

        CommonReferenceString::<R>::new(
            r_pr,
            n_pr,
            norm_bound,
            num_quad_constraints,
            num_constant_quad_constraints,
            &mut thread_rng(),
        )
    }
}

pub fn Z2_to_R_vec<R: Ring>(vec: &Vec<Z2>) -> Vec<R> {
    vec.iter()
        .map(|x| if x.is_zero() { R::zero() } else { R::one() })
        .collect()
}

/// Reinterprets a vector of k = k' * d binary coefficients as k' vectors of d binary coefficients, represented as a vector of k' elements of the polynomial ring R with dimension d.
pub fn lift<R: PolyRing>(vec: &Vector<Z2>) -> Vector<R> {
    let d = R::dimension();
    assert_eq!(
        vec.len() % d,
        0,
        "vector length {} must be multiple of dimension {}",
        vec.len(),
        d
    );
    let coeffs = vec
        .as_slice()
        .chunks(d)
        .map(|chunk| R::from(chunk.to_vec().into_iter().map(embed).collect::<Vec<_>>()))
        .collect();
    Vector::<R>::from_vec(coeffs)
}

/// Upcast an element in Z2 to an element in R
pub fn embed<R: Zero + One>(x: Z2) -> R {
    if x.is_zero() {
        R::zero()
    } else {
        debug_assert!(x.is_one());
        R::one()
    }
}

pub fn basis_vector<R: PolyRing>(i: usize, n: usize) -> Vector<R> {
    assert!(i < n, "i = {} must be less than n = {}", i, n);
    let mut coeffs = vec![R::zero(); n];
    coeffs[i] = R::one();
    Vector::<R>::from_vec(coeffs)
}

// pub fn mul_dense_sparse<R: ConvertibleField>(
//     dense: &Matrix<R>,
//     sparse: &SparseMatrix<R>,
// ) -> Matrix<R> {
//     return (&sparse.transpose() * &dense.transpose()).transpose();
// }

pub struct R1CSInstance<R: Scalar> {
    // TODO: use sparse matrices instead
    pub A: Matrix<R>,
    pub B: Matrix<R>,
    pub C: Matrix<R>,
}
