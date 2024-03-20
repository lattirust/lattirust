#![allow(non_snake_case)]

use std::cmp::{max, max_by};

use log::info;
use rand::thread_rng;
use serde::Serialize;

use lattice_estimator::msis2::MSIS;
use lattice_estimator::norms::Norm;

use crate::labrador::binary_r1cs::util::SECPARAM;
use crate::labrador::prover::Witness;
use crate::labrador::shared::BaseTranscript;
use crate::labrador::util::{chunk_pad, concat, flatten_symmetric_matrix, flatten_vec_vector, mul_matrix_basescalar, shift_right};
use crate::lattice_arithmetic::balanced_decomposition::decompose_balanced_vec_polyring;
use crate::lattice_arithmetic::challenge_set::labrador_challenge_set::LabradorChallengeSet;
use crate::lattice_arithmetic::matrix::{Matrix, sample_uniform_mat, sample_uniform_vec, Vector};
use crate::lattice_arithmetic::poly_ring::PolyRing;
use crate::lattice_arithmetic::traits::WithL2Norm;
use crate::relations::labrador::principal_relation::{ConstantQuadDotProdFunction, PrincipalRelation, QuadDotProdFunction};

/// Common reference string for one round of the LaBRADOR protocol
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
    pub b: u128,
    pub b1: u128,
    pub b2: u128,
    pub next_crs: Option<Box<CommonReferenceString<R>>>,
}

fn round_to_even(x: f64) -> u128 {
    if x.floor() as usize % 2 == 0 { x.floor() as u128 } else { x.ceil() as u128 }
}

impl<R: PolyRing> CommonReferenceString<R> {
    fn t1_b1(decomposition_basis: u128) -> (usize, u128) {
        let log2_q: f64 = R::modulus().next_power_of_two().ilog2() as f64;
        let log2_b = (decomposition_basis as f64).log2();
        let t1 = (log2_q / log2_b).round() as usize;
        let b1 = round_to_even((R::modulus() as f64).powf(1. / t1 as f64));
        (t1, b1)
    }

    fn t2_b2(r: usize, n: usize, beta_sq: f64, decomposition_basis: u128) -> (usize, u128) {
        let d = R::dimension();
        let log2_b = (decomposition_basis as f64).log2();
        let s_std_dev_sq: f64 = beta_sq / ((r * n * d) as f64); // standard deviation of s vectors = beta / sqrt(r * n * d)
        let tmp = f64::sqrt((24 * n * d) as f64) * s_std_dev_sq;
        let t2 = (f64::log2(tmp) / log2_b).round() as usize;
        let b2 = round_to_even(tmp.powf(1. / t2 as f64));
        (t2, b2)
    }

    pub fn new<Rng: rand::Rng + ?Sized>(r: usize, n: usize, mut beta_sq: f64, num_constraints: usize, num_constant_constraints: usize, rng: &mut Rng) -> CommonReferenceString<R> {
        let d = R::dimension();
        let q = R::modulus();
        let log2_q: f64 = q.next_power_of_two().ilog2() as f64;
        let mut beta = beta_sq.sqrt();

        info!("Using Z_q[X]/(X^d+1) with q={q} ({} bits), d={d}", q.next_power_of_two().ilog2());
        info!("Setting CRS parameters for n={n}, r={r}, d={d}, beta={beta:.1}, num_constraints={num_constraints}, num_constant_constraints={num_constant_constraints}");

        let MAX_RECURSION_DEPTH = 7;
        beta_sq = beta_sq * f64::sqrt(128./30.).powi(MAX_RECURSION_DEPTH);
        beta = beta_sq.sqrt();
        info!("  Accounting for at most {MAX_RECURSION_DEPTH} recursion levels, use beta={beta:.1}");

        // Checks
        assert!(beta < f64::sqrt(30. / 128.) * (q as f64) / 125.);

        // Set decomposition basis
        let s_sq = beta / ((r * n * d) as f64);
        let s = s_sq.sqrt(); // standard deviation of the Z_q coefficients of the s vectors
        info!("  s={s} (std deviation of coefficients of s-vector)");

        let tau = LabradorChallengeSet::<R>::VARIANCE_SUM_COEFFS;
        let b = round_to_even((s * (12. * r as f64 * tau).sqrt()).sqrt());
        info!("  b={b} (main decomposition basis)");
        assert!(b > 1);

        // Set auxiliary decomposition bases
        let (t1, b1) = Self::t1_b1(b);
        info!("  b1={b1}, t1={t1} (first decomposition basis and decomposition length)");

        let (t2, b2) = Self::t2_b2(r, n, beta_sq, b);
        info!("  b2={b2}, t2={t2} (first decomposition basis and decomposition length)");

        let num_aggregs = (128. / log2_q).ceil() as usize;
        info!("  num_aggregs={num_aggregs} (ceil(128/log(q)))");

        // Compute the norm bound for the next folded instance
        let beta_prime = |k| {
            Self::next_norm_bound_sq(r, n, beta_sq, k, b).sqrt()
        };

        let op_norm = LabradorChallengeSet::<R>::OPERATOR_NORM_THRESHOLD;
        let norm_bound_1 = |kappa| {
            // max(8T(b + 1)β′, 2(b + 1)β′ + 4T sqrt(128/30)β)
            max_by(8. * op_norm * (b + 1) as f64 * beta_prime(kappa),
                   2. * (b + 1) as f64 * beta_prime(kappa) + 4. * op_norm * f64::sqrt(128. / 30.) * beta,
                   f64::total_cmp,
            )
        };

        // Ensure MSIS_{n=k, d, q, beta_1, m=n} is hard (l_2 norm)
        let mut msis_1 = MSIS {
            n: 0, // Dummy value, will be set later
            d,
            q,
            length_bound: 0., // Dummy value
            m: n,
            norm: Norm::L2,
        };
        let k = msis_1.upper_bound_n();
        // let k = msis_1.find_optimal_n_dynamic(norm_bound_1, SECPARAM).expect(format!("failed to find secure rank for {msis_1}. Are there enough constraints in your system?").as_str());
        msis_1 = msis_1.with_n(k).with_length_bound(norm_bound_1(k));
        //info!("  k={k} for the MSIS instance {msis_1}  gives {} bits of security",msis_1.security_level()); // TODO: silently assume that this gives us enough security, which it will for any reasonable parameters
        info!("  Chose largest k={k} for the MSIS instance {msis_1}");

        let mut msis_2 = MSIS {
            n: 0, // Dummy value, will be set later
            d,
            q,
            length_bound: 2. * beta_prime(k),
            m: k,
            norm: Norm::L2,
        };
        let k1 = msis_2.find_optimal_n(SECPARAM).expect(format!("failed to find secure rank for {msis_2}. Are there enough constraints in your system?").as_str());
        let k2 = k1;
        msis_2 = msis_2.with_n(k1).with_length_bound(2. * beta_prime(k));
        info!("  k1=k2={k1} for the MSIS instance {msis_2}  gives {} bits of security", msis_2.security_level());

        // TODO: this only gives 125 bits of soundness error for SECPARAM = 128, how do we best document this?
        // TODO: all params should be bigger to account for the slack of sqrt(128/30) per recursion level

        let mut crs = CommonReferenceString {
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
            A: sample_uniform_mat(k, n, rng),
            B: vec![vec![sample_uniform_mat(k1, k, rng); t1]; r],
            C: (0..r).map(
                |i| (0..i + 1).map(
                    |_| vec![sample_uniform_vec(k2, rng); t2]
                ).collect()
            ).collect(),
            D: (0..r).map(
                |i| (0..i + 1).map(
                    |_| vec![sample_uniform_vec(k2, rng); t1]
                ).collect()
            ).collect(),
            b,
            b1,
            b2,
            next_crs: None,
        };
        crs.next_crs = crs.next_crs().map(|crs| Box::new(crs));
        crs
    }

    /// Compute the squared norm bound for the next folded instance (cf. Section 5.4 of the Labrador paper)
    pub fn next_norm_bound_sq(r: usize, n: usize, norm_bound_squared: f64, k: usize, decomposition_basis: u128) -> f64 {
        let b_sq = decomposition_basis * decomposition_basis;
        let d = R::dimension();
        let challenge_variance = LabradorChallengeSet::<R>::VARIANCE_SUM_COEFFS;

        let (t1, b1) = Self::t1_b1(decomposition_basis);
        let (t2, b2) = Self::t2_b2(r, n, norm_bound_squared, decomposition_basis);

        let gamma_sq = norm_bound_squared * challenge_variance;
        let gamma_1_sq = (b1 * b1 * t1 as u128) as f64 / 12. * (r * k * d) as f64 + (b2 * b2 * t2 as u128) as f64 / 12. * ((r * (r + 1)).div_ceil(2) * d) as f64;
        let gamma_2_sq = (b1 * b1 * t1 as u128) as f64 / 12. * ((r * (r + 1)).div_ceil(2) * d) as f64;
        let beta_next_sq: f64 = (2. / b_sq as f64) * gamma_sq + gamma_1_sq * gamma_2_sq;
        beta_next_sq
    }

    /// Returns true if there should be a next round in the LaBRADOR protocol, and false if `crs' is to be used as the base case.
    fn recurse(&self) -> bool {
        self.last_prover_message_size() > self.proof_size()
    }

    pub fn proof_size(&self) -> usize {
        let log_q = R::modulus().next_power_of_two().ilog2() as usize;
        let beta = self.norm_bound_squared.sqrt();
        // TODO: this assumes that the challenges are sampled from a PRF with a 128-bit seed.
        (self.k1 + self.k2) * self.d * log_q // outer commitments
            + 256 * ((12. * beta) / 2f64.sqrt()).log2().ceil() as usize // JL projection
            + self.num_aggregs * self.d * log_q // JL proof
            + 4 * 128 // challenges
    }

    pub fn last_prover_message_size(&self) -> usize {
        let log_q = R::modulus().next_power_of_two().ilog2() as usize;
        let beta = self.norm_bound_squared.sqrt();
        let tau = LabradorChallengeSet::<R>::VARIANCE_SUM_COEFFS;
        let tmp = ((self.r * (self.r + 1)) / 2) * self.d;
        self.n * self.d * f64::log2(12. * beta * f64::sqrt(tau / (self.n * self.d) as f64)).ceil() as usize // z
            + self.r * self.k * self.d * log_q // t
            + tmp * f64::log2(12. * self.norm_bound_squared * f64::sqrt(2. / (self.r * self.r * self.n * self.d) as f64)).ceil() as usize // g
            + tmp * log_q // h
    }

    fn next_crs(&self) -> Option<CommonReferenceString<R>> {
        if self.recurse() { return None; }

        let crs = self;

        // Generate instance for next iteration of the protocol
        // Set nu, mu such that 2 * n_next ≈ m_next,
        // i.e., set mu = 2 and nu ≈ 2 * n / m
        let m = crs.r * crs.t1 * crs.k + (crs.t1 + crs.t2) * (crs.r * (crs.r + 1)).div_ceil(2);
        let mu = 2;
        let nu = (2. * crs.n as f64 / m as f64).round() as usize;
        let n_next = max(crs.n.div_ceil(nu), m.div_ceil(mu));
        let _m_next = (m as f64 / 2.).round() as usize;
        let r_next = 2 * nu + m;


        let num_quad_constraints = self.k + // k constraints for Az = sum_{i in [r]} c_i t_i (one per row)
            1 + // 1 constraint for <z, z> = sum_{i, j in [r]} g_ij c_i c_j
            1 + // 1 constraint for sum_{i in [r]} <phi_i, z> c_i = sum_{i, j in [r]} h_ij c_i c_j
            1 + // 1 constraint for sum_{i,j in [r]} a_ij * g_ij + sum_{i in [r]} h_ii = b
            crs.k1 + // k1 constraints for u_1
            crs.k2; // k2 constraints for u_2

        let next_norm_bound_squared = CommonReferenceString::<R>::next_norm_bound_sq(crs.r, crs.n, crs.norm_bound_squared, crs.k, crs.b);

        Some(CommonReferenceString::<R>::new(r_next, n_next, next_norm_bound_squared, num_quad_constraints, 0, &mut thread_rng()))
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
            witness.s.iter().map(|s_i| s_i.l2_norm_squared()).sum::<u64>() as f64 <= self.norm_bound_squared
    }
}

pub fn fold_instance<'a, R: PolyRing>(transcript: &BaseTranscript<R>, compute_witness: bool) -> (PrincipalRelation<R>, Option<Witness<R>>) {
    // assert!(transcript.crs.recurse());
    let crs = transcript.crs;

    // Generate instance for next iteration of the protocol
    // Set nu, mu such that 2 * n_next ≈ m_next,
    // i.e., set mu = 2 and nu ≈ 2 * n / m
    let m = crs.r * crs.t1 * crs.k + (crs.t1 + crs.t2) * (crs.r * (crs.r + 1)).div_ceil(2);
    let mu = 2;
    let nu = (2. * crs.n as f64 / m as f64).round() as usize;
    let n_next = max(crs.n.div_ceil(nu), m.div_ceil(mu));
    let _m_next = (m as f64 / 2.).round() as usize;
    let r_next = 2 * nu + m;


    let mut quad_dot_prod_funcs_next = Vec::<QuadDotProdFunction<R>>::with_capacity(crs.k + crs.k1 + crs.k2 + 3);

    // The new witness is the concatenation of the flattened vectors z^(0), z^(1), (t_i^(j))_{i in [r], j in [t1]}, (g_ij)_{i >= j in [r]}, (h_ij)_{i >= j in [r]}

    // │<---------------------- r * t1 * n -------------------->│
    // ┌─────────────┬─────────────┬─────────────┬──────────────┐        ┌─────────────┬─────────────┬─────────────┬──────────────┐
    // │    t_0^(0)  │    t_0^(1)  │     ...     │t_{r-1}^(t1-1)│<-o_tg->│     g_0,0   │     g_0,1   │     ...     │     g_r,r    │
    // └─────────────┴─────────────┴─────────────┴──────────────┘        └─────────────┴─────────────┴─────────────┴──────────────┘
    // ┌────────┬────────┬────────┬                          ┬───────────┬─────────┬
    // │  s_2nu │s_2nu+1 │s_2nu+2 │             ...          │s_2nu+tau-1│s_2nu+tau│   ...
    // └────────┴────────┴────────┘                          └───────────┴─────────┴
    // │<------------------------ tau * n_next ------------------------->│
    // tau = ceil(r * t1 * n / n_next)
    let tau = (crs.r * crs.t1 * crs.n).div_ceil(n_next); // Number of n_next-sized s-vectors over R_q needed to represent all t_i^(j)'s

    let offset_t_g = tau * n_next - crs.r * crs.t1 * crs.n; // Offset to align the t_i^(j)'s with the g_ij's (o_tg in the diagram above)
    let gamma = ((crs.r * (crs.r + 1)).div_ceil(2) * crs.t2 + offset_t_g).div_ceil(n_next); // Number of n_next-sized s-vectors over R_q needed to represent all g_ij's

    let offset_g_h = gamma * n_next - ((crs.r * (crs.r + 1)).div_ceil(2) * crs.t2 + offset_t_g); // Offset to align the g_ij's with the h_ij's
    let eta = ((crs.r * (crs.r + 1)).div_ceil(2) * crs.t1 + offset_g_h).div_ceil(n_next); // Number of n_next-sized s-vectors over R_q needed to represent all h_ij's

    let c = transcript.c.as_ref().expect("c not available");

    let b_ring = R::from(crs.b as u128);
    let mut b_pows = Vec::<R>::with_capacity(max(crs.t1, crs.t2));
    b_pows[0] = R::one();
    for i in 1..max(crs.t1, crs.t2) {
        b_pows.push(b_pows[i - 1] * b_pows[i - 1]);
    }
    let two = R::from(2u128);

    // Constraints for Az = sum_{i in [r]} c_i t_i
    {
        for l in 0..crs.k {
            let mut phis_next = vec![Vector::<R>::zeros(n_next); r_next];
            let A_l = crs.A.row(l).transpose(); // <A_l, z_0>
            phis_next[0..nu].clone_from_slice(crate::labrador::util::split(&A_l, nu).as_slice()); // <A_l * b, z_1>
            phis_next[nu..2 * nu].clone_from_slice(crate::labrador::util::split(&(&A_l * b_ring), nu).as_slice());

            // <(t_j^(k)_l)_{j, k}, (c_j * b^k)_{j, k})>
            let mut c_vec = vec![R::zero(); crs.r * crs.t1 * crs.n];
            for i in 0..crs.r {
                for k in 0..crs.t1 {
                    c_vec[i * crs.r * crs.t1 * crs.n + k * crs.t1 * crs.n + l] = c[i] * b_pows[k];
                }
            }
            let c_vec_split = chunk_pad(&Vector::<R>::from_vec(c_vec), n_next);
            assert_eq!(c_vec_split.len(), tau);
            phis_next[2 * nu..2 * nu + tau].clone_from_slice(c_vec_split.as_slice());

            quad_dot_prod_funcs_next.push(QuadDotProdFunction::<R>::new(Matrix::<R>::zeros(r_next, r_next), phis_next, R::zero()));
        }
    }

    // Constraint for <z, z> = sum_{i, j in [r]} g_ij c_i c_j
    {
        let mut A = Matrix::<R>::zeros(r_next, r_next);
        let b_sq = R::from((crs.b * crs.b) as u128);
        for i in 0..r_next {
            A[(i, i)] = R::one();
            A[(i, i + nu)] = b_ring;
            A[(i + nu, i)] = A[(i, i + nu)];
            A[(i + nu, i + nu)] = b_sq;
        }
        let mut phis_next = vec![Vector::<R>::zeros(n_next); r_next];

        let mut c_prods = Vec::<R>::with_capacity((crs.r * (crs.r + 1)).div_ceil(2) * crs.t2);
        for i in 0..crs.r {
            for j in 0..i + 1 {
                for k in 0..crs.t2 {
                    if i == j {
                        c_prods.push(c[i] * c[i] * b_pows[k]);
                    } else {
                        c_prods.push(c[i] * c[j] * b_pows[k] * two);
                    }
                }
            }
        }
        c_prods = shift_right(&c_prods, offset_t_g);
        let c_prods_flat_split = chunk_pad(&Vector::<R>::from_vec(c_prods), n_next);
        assert_eq!(c_prods_flat_split.len(), gamma);
        phis_next[2 * nu + tau - 1..2 * nu + tau - 1 + gamma].clone_from_slice(c_prods_flat_split.as_slice());

        quad_dot_prod_funcs_next.push(QuadDotProdFunction::<R>::new(A, phis_next, R::zero()));
    }

    // Constraint for sum_{i in [r]} <phi_i, z> c_i = sum_{i, j in [r]} h_ij c_i c_j
    {
        let phi = transcript.phi.as_ref().expect("phi not available");
        debug_assert_eq!(phi.len(), crs.n);
        let mut phis_next = vec![Vector::<R>::zeros(n_next); r_next];

        let mut phi_lc_0 = Vector::<R>::zeros(crs.n);
        let mut phi_lc_1 = Vector::<R>::zeros(crs.n);

        for i in 0..crs.r {
            phi_lc_0 += &phi[i] * c[i];
            phi_lc_1 += &phi[i] * c[i] * b_ring;
        }
        // <sum_{i in [r]} phi_i * c_i, z_0>
        phis_next[0..nu].clone_from_slice(crate::labrador::util::split(&phi_lc_0, nu).as_slice());
        // <sum_{i in [r]} phi_i * c_i * b, z_1>
        phis_next[nu..2 * nu].clone_from_slice(crate::labrador::util::split(&phi_lc_1, nu).as_slice());

        // <(c_i * c_j * b^k)_{i, j, k}, h>
        let mut c_prods = Vec::<R>::with_capacity((crs.r * (crs.r + 1)).div_ceil(2) * crs.t2);
        for i in 0..crs.r {
            for j in 0..i + 1 {
                for k in 0..crs.t1 {
                    if i == j {
                        c_prods.push(c[i] * c[i] * b_pows[k]);
                    } else {
                        c_prods.push(c[i] * c[j] * b_pows[k] * two);
                    }
                }
            }
        }
        c_prods = shift_right(&c_prods, offset_g_h);
        let c_prods_flat_split = chunk_pad(&Vector::<R>::from_vec(c_prods), n_next);
        assert_eq!(c_prods_flat_split.len(), eta);
        phis_next[2 * nu + tau + gamma - 1..2 * nu + tau + gamma - 1 + eta].clone_from_slice(c_prods_flat_split.as_slice());

        quad_dot_prod_funcs_next.push(QuadDotProdFunction::<R>::new(Matrix::<R>::zeros(r_next, r_next), phis_next, R::zero()));
    }

    // Constraint for sum_{i,j in [r]} a_ij * g_ij + sum_{i in [r]} h_ii = b
    {
        let b__ = transcript.b__.as_ref().expect("b'' not available");
        let alpha = transcript.alpha.as_ref().expect("alpha not available");
        let beta = transcript.beta.as_ref().expect("beta not available");
        let psi = transcript.psi.as_ref().expect("psi not available");

        // Compute b = sum_{k in [K]} alpha_k * b^(k) + sum_{k in [K']} beta_k * b''^(k)
        let mut b = R::zero();
        for k in 0..crs.num_constraints {
            b += alpha[k] * transcript.instance.quad_dot_prod_funcs[k].b;
        }
        for k in 0..crs.num_aggregs {
            b += beta[k] * b__[k];
        }
        // Compute a_ij = sum_{k in [K]} alpha_k * a_ij^(k) + sum_{k in [K']} beta_k * a_ij''^(k), where a_ij''^(k) = sum_{l in [L]} psi_l^(k) * a_ij'^(l)
        let mut A = Matrix::<R>::zeros(crs.r, crs.r);
        for k in 0..crs.num_constraints {
            if let Some(ref A_k) = transcript.instance.quad_dot_prod_funcs[k].A {
                A += A_k * alpha[k];
            }
        }
        for k in 0..crs.num_aggregs {
            let mut a_k__ = Matrix::<R>::zeros(crs.r, crs.r);
            for l in 0..crs.num_constant_constraints {
                if let Some(ref A_l_) = transcript.instance.ct_quad_dot_prod_funcs[l].A {
                    a_k__ += mul_matrix_basescalar(A_l_, psi[k][l]);
                }
            }
            A += a_k__ * beta[k];
        }
        let mut phis_next = vec![Vector::<R>::zeros(n_next); r_next];
        // Set phis for a_ij * g_ij
        let mut A_vec = Vec::<R>::with_capacity((crs.r * (crs.r + 1)).div_ceil(2));
        for i in 0..crs.r {
            for j in 0..i {
                A_vec.push(A[(i, j)] + A[(j, i)]);
            }
            A_vec.push(A[(i, i)]);
        }
        phis_next[2 * nu + tau..2 * nu + tau + gamma].clone_from_slice(chunk_pad(&Vector::<R>::from_vec(A_vec), n_next).as_slice());

        // Set phis for 1 * h_ii, i.e.
        let mut indic = Vector::<R>::zeros(crs.r * (crs.r + 1).div_ceil(2));
        for i in 0..crs.r {
            indic[((i + 1) * (i + 2)) / 2] = R::one();
        }
        phis_next[2 * nu + tau + gamma..2 * nu + tau + gamma + eta].clone_from_slice(chunk_pad(&indic, n_next).as_slice());


        quad_dot_prod_funcs_next.push(QuadDotProdFunction::<R>::new(Matrix::<R>::zeros(r_next, r_next), phis_next, b));
    }

    // Constraints for u_1
    let u_1 = transcript.u_1.as_ref().expect("u_1 not available");
    for l in 0..crs.k1 {
        let mut phis = vec![Vector::<R>::zeros(n_next); r_next];
        let mut B_flat = Vec::<R>::with_capacity(crs.r * crs.t1 * crs.n);
        for i in 0..crs.r {
            for k in 0..crs.t1 {
                B_flat.extend(crs.B[i][k].row(l).transpose().data.as_vec());
            }
        }
        let B_flat_split = chunk_pad(&Vector::<R>::from_vec(B_flat), n_next);
        phis[2 * nu..2 * nu + tau].clone_from_slice(B_flat_split.as_slice());

        let mut C_flat = Vec::<R>::with_capacity((crs.r * (crs.r + 1)).div_ceil(2) * crs.t2);
        for i in 0..crs.r {
            for j in 0..i + 1 {
                for k in 0..crs.t2 {
                    C_flat.push(crs.C[i][j][k][l]);
                }
            }
        }
        C_flat = shift_right(&C_flat, offset_t_g);
        let C_flat_split = chunk_pad(&Vector::<R>::from_vec(C_flat), n_next);
        assert_eq!(C_flat_split.len(), gamma);
        phis[2 * nu + tau - 1..2 * nu + tau - 1 + gamma].clone_from_slice(C_flat_split.as_slice());


        quad_dot_prod_funcs_next.push(QuadDotProdFunction::<R>::new(Matrix::<R>::zeros(r_next, r_next), phis, u_1[l]));
    }

    // Constraints for u_2
    let u_2 = transcript.u_2.as_ref().expect("u_2 not available");
    for l in 0..crs.k2 {
        let mut phis = vec![Vector::<R>::zeros(n_next); r_next];

        let mut D_flat = Vec::<R>::with_capacity((crs.r * (crs.r + 1)).div_ceil(2) * crs.t1);
        for i in 0..crs.r {
            for j in 0..i + 1 {
                for k in 0..crs.t1 {
                    D_flat.push(crs.D[i][j][k][l]);
                }
            }
        }
        let D_flat = shift_right(&D_flat, offset_g_h);
        let D_flat_split = chunk_pad(&Vector::<R>::from_vec(D_flat), n_next);
        assert_eq!(D_flat_split.len(), eta);
        phis[2 * nu + tau + gamma - 1..2 * nu + tau + gamma - 1 + eta].clone_from_slice(D_flat_split.as_slice());

        quad_dot_prod_funcs_next.push(QuadDotProdFunction::<R>::new(Matrix::<R>::zeros(r_next, r_next), phis, u_2[l]));
    }


    debug_assert_eq!(crs.next_crs.as_ref().unwrap().num_constraints, quad_dot_prod_funcs_next.len());
    debug_assert_eq!(crs.next_crs.as_ref().unwrap().num_constant_constraints, 0);

    let instance_next = PrincipalRelation::<R> {
        quad_dot_prod_funcs: quad_dot_prod_funcs_next,
        ct_quad_dot_prod_funcs: vec![],
    };

    let witness = if compute_witness {
        let (z, t, G, H) = (transcript.z.as_ref().unwrap(), transcript.t.as_ref().unwrap(), transcript.G.as_ref().unwrap(), transcript.H.as_ref().unwrap());
        let z_decomp = decompose_balanced_vec_polyring(&z, crs.b, Some(2usize));
        assert_eq!(z_decomp.len(), 2);

        let v = concat(&[&flatten_vec_vector(&t).as_slice(), &flatten_symmetric_matrix(&G).as_slice(), &flatten_symmetric_matrix(&H).as_slice()]);
        let z_0_split = crate::labrador::util::split(&z_decomp[0], nu);
        let z_1_split = crate::labrador::util::split(&z_decomp[1], nu);
        let v_split = crate::labrador::util::split(&v, mu);
        Some(Witness::<R> {
            s: vec![z_0_split, z_1_split, v_split].concat(),
        })
    } else {
        None
    };

    (instance_next, witness)
}