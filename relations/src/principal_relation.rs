#![allow(non_snake_case)]

use std::cmp::{max, max_by};

use ark_std::rand;
use ark_std::rand::thread_rng;
use log::info;
use num_traits::{One, ToPrimitive, Zero};
use serde::Serialize;

use lattice_estimator::msis::{MSIS, msis_h_128_l2};
use lattice_estimator::norms::Norm;
use lattirust_arithmetic::challenge_set::labrador_challenge_set::LabradorChallengeSet;
use lattirust_arithmetic::linear_algebra::inner_products::inner_products;
use lattirust_arithmetic::linear_algebra::Matrix;
use lattirust_arithmetic::linear_algebra::Vector;
use lattirust_arithmetic::ring::PolyRing;
use lattirust_arithmetic::traits::WithL2Norm;

use crate::Relation;

#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct Index<R: PolyRing> {
    /// Number of witness vectors
    pub r: usize,
    /// Number of entries in a witness vector
    pub n: usize,
    /// Dimension of the polynomial ring Z_q[X]/(X^d + 1)
    pub d: usize,
    /// Square of the L2-norm bound on the concatenation of witness vectors
    pub norm_bound_squared: f64,
    /// Size of the first-level commitment
    pub k: usize,
    /// Size of the second-level commitment for elements decomposed in basis `b1`
    pub k1: usize,
    /// Size of the second-level commitment for elements decomposed in basis `b2`
    pub k2: usize,
    /// Length of decompositions in basis `b1`
    pub t1: usize,
    /// Length of decompositions in basis `b2`
    pub t2: usize,
    /// Number of aggregated constraints when reducing constant constraints, equal to  `security parameter / log(q)`
    pub num_aggregs: usize,
    /// Number of quadratic-linear constraints
    pub num_constraints: usize,
    /// Number of quadratic-linear constraints on constant coefficients
    pub num_constant_constraints: usize,
    /// First-level commitment matrix, of size k x n
    pub A: Matrix<R>,
    /// Second-level commitment matrices for first decomposition basis, (r x t1) instances, each of size k1 x k
    pub B: Vec<Vec<Matrix<R>>>,
    /// Second-level commitment matrices for second decomposition basis, (r x r x t2) instances, each of size k2 x 1
    pub C: Vec<Vec<Vec<Vector<R>>>>,
    /// Second-level commitment matrices for first decomposition basis, (r x r x t1) instances, each of size k2 x 1
    pub D: Vec<Vec<Vec<Vector<R>>>>,
    /// Decomposition basis for z-vectors, roughly equal to `b1` and `b2`
    pub b: u128,
    /// Decomposition basis for first-level commitments (t-vectors), roughly equal to `b` and `b2`
    pub b1: u128,
    /// Decomposition basis for inner product terms (g), roughly equal to `b1` and `b2`
    pub b2: u128,
    /// A reference to the CRS for the next recursive round, or `None` if this is the CRS for the last round
    pub next_crs: Option<Box<Index<R>>>,
}

impl<R: PolyRing> Index<R> {
    pub fn floor_to_even(x: f64) -> u128 {
        if x.floor() as usize % 2 == 0 {
            x.floor() as u128
        } else {
            x.floor() as u128 - 1
        }
    }

    pub fn round_to_even(x: f64) -> u128 {
        if x.floor() as usize % 2 == 0 {
            x.floor() as u128
        } else if x == x.trunc() {
            // x is an odd integer
            x.ceil() as u128 + 1
        } else {
            x.ceil() as u128
        }
    }

    fn t1_b1(decomposition_basis: u128) -> (usize, u128) {
        let log2_q: f64 = R::modulus().bits() as f64;
        let log2_b = (decomposition_basis as f64).log2();
        let t1 = (log2_q / log2_b).round() as usize;
        let b1 = Self::round_to_even(R::modulus().nth_root(t1 as u32).to_f64().unwrap());
        (t1, b1)
    }

    fn t2_b2(r: usize, n: usize, beta_sq: f64, decomposition_basis: u128) -> (usize, u128) {
        let d = R::dimension();
        let log2_b = (decomposition_basis as f64).log2();
        let s_std_dev_sq: f64 = beta_sq / ((r * n * d) as f64); // standard deviation of s vectors = beta / sqrt(r * n * d)
        let tmp = f64::sqrt((24 * n * d) as f64) * s_std_dev_sq;
        let t2 = (f64::log2(tmp) / log2_b).round() as usize;
        let b2 = Self::round_to_even(tmp.powf(1. / t2 as f64));
        (t2, b2)
    }

    pub fn new<Rng: rand::Rng + ?Sized>(
        r: usize,
        n: usize,
        mut beta_sq: f64,
        num_constraints: usize,
        num_constant_constraints: usize,
        rng: &mut Rng,
    ) -> Index<R> {
        let d = R::dimension();
        let q = R::modulus();
        let log2_q: f64 = q.bits() as f64;
        let mut beta = beta_sq.sqrt();

        info!("Using Z_q[X]/(X^d+1) with q={q} ({} bits), d={d}", log2_q);
        info!("Setting CRS parameters for n={n}, r={r}, d={d}, beta={beta:.1}, num_constraints={num_constraints}, num_constant_constraints={num_constant_constraints}");

        let MAX_RECURSION_DEPTH = 7;
        beta_sq *= f64::sqrt(128. / 30.).powi(MAX_RECURSION_DEPTH);
        beta = beta_sq.sqrt();
        info!(
            "  Accounting for at most {MAX_RECURSION_DEPTH} recursion levels, use beta={beta:.1}"
        );

        // Checks
        assert!(beta < f64::sqrt(30. / 128.) * (q.to_f64().unwrap()) / 125.);

        // Set decomposition basis
        let s_sq = beta / ((r * n * d) as f64);
        let s = s_sq.sqrt();
        // standard deviation of the Z_q coefficients of the s vectors
        info!("  s={s} (std deviation of coefficients of s-vector)");

        let tau = LabradorChallengeSet::<R>::VARIANCE_SUM_COEFFS;
        let b = Self::round_to_even((s * (12. * r as f64 * tau).sqrt()).sqrt());
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
        let beta_prime = |k| Self::next_norm_bound_sq(r, n, beta_sq, k, b).sqrt();

        let op_norm = LabradorChallengeSet::<R>::OPERATOR_NORM_THRESHOLD;
        let norm_bound_1 = |kappa| {
            // max(8T(b + 1)β′, 2(b + 1)β′ + 4T sqrt(128/30)β)
            max_by(
                8. * op_norm * (b + 1) as f64 * beta_prime(kappa),
                2. * (b + 1) as f64 * beta_prime(kappa)
                    + 4. * op_norm * f64::sqrt(128. / 30.) * beta,
                f64::total_cmp,
            )
        };

        // Ensure MSIS_{n=k, d, q, beta_1, m=n} is hard (l_2 norm)
        let mut msis_1 = MSIS {
            h: 0, // Dummy value, will be set later
            d,
            q: q.clone(),
            length_bound: 0., // Dummy value
            w: n,
            norm: Norm::L2,
        };
        let k = msis_1.upper_bound_h();
        // let k = msis_1.find_optimal_h_dynamic(norm_bound_1, SECPARAM).expect(format!("failed to find secure rank for {msis_1}. Are there enough constraints in your system?").as_str());
        msis_1 = msis_1.with_h(k).with_length_bound(norm_bound_1(k));
        //info!("  k={k} for the MSIS instance {msis_1}  gives {} bits of security",msis_1.security_level()); // TODO: silently assume that this gives us enough security, which it will for any reasonable parameters
        info!("  Chose largest k={k} for the MSIS instance {msis_1}");

        let mut msis_2 = MSIS {
            h: 0, // Dummy value, will be set later
            d,
            q: q.clone(),
            length_bound: 2. * beta_prime(k),
            w: k,
            norm: Norm::L2,
        };
        let k1 = msis_h_128_l2(&msis_2).unwrap(); // TODO: Switch back to lattice estimator
                                                  // let k1 = find_optimal_h(&msis_2, SECURITY_PARAMETER).expect(format!("failed to find secure rank for {msis_2}. Are there enough constraints in your system?").as_str());
        let k2 = k1;
        msis_2 = msis_2.with_h(k1).with_length_bound(2. * beta_prime(k));
        info!(
            "  k1=k2={k1} for the MSIS instance {msis_2}  gives {} bits of security",
            msis_2.security_level()
        );

        // TODO: this only gives 125 bits of soundness error for SECPARAM = 128, how do we best document this?
        // TODO: all params should be bigger to account for the slack of sqrt(128/30) per recursion level

        let mut crs = Index {
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
            A: Matrix::<R>::par_rand(k, n),
            B: vec![vec![Matrix::<R>::rand(k1, k, rng); t1]; r],
            C: (0..r)
                .map(|i| {
                    (0..i + 1)
                        .map(|_| vec![Vector::<R>::rand(k2, rng); t2])
                        .collect()
                })
                .collect(),
            D: (0..r)
                .map(|i| {
                    (0..i + 1)
                        .map(|_| vec![Vector::<R>::rand(k2, rng); t1])
                        .collect()
                })
                .collect(),
            b,
            b1,
            b2,
            next_crs: None,
        };
        crs.next_crs = crs.next_crs().map(Box::new);
        crs
    }

    /// Compute the squared norm bound for the next folded instance (cf. Section 5.4 of the Labrador paper)
    pub fn next_norm_bound_sq(
        r: usize,
        n: usize,
        norm_bound_squared: f64,
        k: usize,
        decomposition_basis: u128,
    ) -> f64 {
        let b_sq = decomposition_basis * decomposition_basis;
        let d = R::dimension();
        let challenge_variance = LabradorChallengeSet::<R>::VARIANCE_SUM_COEFFS;

        let (t1, b1) = Self::t1_b1(decomposition_basis);
        let (t2, b2) = Self::t2_b2(r, n, norm_bound_squared, decomposition_basis);

        let gamma_sq = norm_bound_squared * challenge_variance;
        let gamma_1_sq = (b1 * b1 * t1 as u128) as f64 / 12. * (r * k * d) as f64
            + (b2 * b2 * t2 as u128) as f64 / 12. * ((r * (r + 1)).div_ceil(2) * d) as f64;
        let gamma_2_sq =
            (b1 * b1 * t1 as u128) as f64 / 12. * ((r * (r + 1)).div_ceil(2) * d) as f64;
        let beta_next_sq: f64 = (2. / b_sq as f64) * gamma_sq + gamma_1_sq * gamma_2_sq;
        beta_next_sq
    }

    /// Returns true if there should be a next round in the LaBRADOR protocol, and false if `crs' is to be used as the base case.
    fn recurse(&self) -> bool {
        self.last_prover_message_size() > self.proof_size()
    }

    pub fn proof_size(&self) -> usize {
        let log_q = R::modulus().bits() as usize;
        let beta = self.norm_bound_squared.sqrt();
        // TODO: this assumes that the challenges are sampled from a PRF with a 128-bit seed.
        (self.k1 + self.k2) * self.d * log_q // outer commitments
            + 256 * ((12. * beta) / 2f64.sqrt()).log2().ceil() as usize // JL projection
            + self.num_aggregs * self.d * log_q // JL proof
            + 4 * 128 // challenges
    }

    pub fn last_prover_message_size(&self) -> usize {
        let log_q = R::modulus().bits() as usize;
        let beta = self.norm_bound_squared.sqrt();
        let tau = LabradorChallengeSet::<R>::VARIANCE_SUM_COEFFS;
        let tmp = ((self.r * (self.r + 1)) / 2) * self.d;
        self.n * self.d * f64::log2(12. * beta * f64::sqrt(tau / (self.n * self.d) as f64)).ceil() as usize // z
            + self.r * self.k * self.d * log_q // t
            + tmp * f64::log2(12. * self.norm_bound_squared * f64::sqrt(2. / (self.r * self.r * self.n * self.d) as f64)).ceil() as usize // g
            + tmp * log_q // h
    }

    fn next_crs(&self) -> Option<Index<R>> {
        if self.recurse() {
            return None;
        }

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

        let next_norm_bound_squared =
            Index::<R>::next_norm_bound_sq(crs.r, crs.n, crs.norm_bound_squared, crs.k, crs.b);

        Some(Index::<R>::new(
            r_next,
            n_next,
            next_norm_bound_squared,
            num_quad_constraints,
            0,
            &mut thread_rng(),
        ))
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

    pub fn is_wellformed_instance(&self, instance: &Instance<R>) -> bool {
        instance.quad_dot_prod_funcs.len() == self.num_constraints
            && instance
                .quad_dot_prod_funcs
                .iter()
                .all(|c| self.is_wellformed_constraint(c))
            && instance.ct_quad_dot_prod_funcs.len() == self.num_constant_constraints
            && instance
                .ct_quad_dot_prod_funcs
                .iter()
                .all(|c| self.is_wellformed_const_constraint(c))
    }

    pub fn is_wellformed_witness(&self, witness: &Witness<R>) -> bool {
        witness.s.len() == self.r
            && witness.s.iter().all(|s_i| s_i.len() == self.n)
            && witness
                .s
                .iter()
                .map(|s_i| s_i.l2_norm_squared())
                .sum::<u128>() as f64
                <= self.norm_bound_squared
    }
}

#[derive(Clone, Eq, PartialEq, Debug)]
pub struct Instance<R: PolyRing> {
    pub quad_dot_prod_funcs: Vec<QuadDotProdFunction<R>>,
    pub ct_quad_dot_prod_funcs: Vec<ConstantQuadDotProdFunction<R>>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct QuadDotProdFunction<R: PolyRing> {
    // TODO: A is always symmetric, so we could at least use a symmetric matrix type. A is also very sparse in some cases.
    pub A: Option<Matrix<R>>,
    // TODO: phi can be quite sparse
    pub phi: Vec<Vector<R>>,
    pub b: R,
    _private: (), // Forbid direct initialization, force users to use new(), which does some basis debug_asserts
}

impl<R: PolyRing> QuadDotProdFunction<R> {
    pub fn new(A: Matrix<R>, phi: Vec<Vector<R>>, b: R) -> Self {
        let (r, n) = (A.nrows(), phi[0].len());
        debug_assert_eq!(A.ncols(), r, "A should be square");
        debug_assert_eq!(A.transpose(), A, "A should be symmetric");

        debug_assert_eq!(
            phi.len(),
            r,
            "phi should have the same length as the dimensions of A"
        );
        debug_assert!(
            phi.iter().all(|phi_i| phi_i.len() == n),
            "each phi_i should have the same length"
        );
        Self {
            A: Some(A),
            phi,
            b,
            _private: (),
        }
    }

    pub fn new_linear(phi: Vec<Vector<R>>, b: R) -> Self {
        let n = phi[0].len();
        debug_assert!(
            phi.iter().all(|phi_i| phi_i.len() == n),
            "each phi_i should have the same length"
        );
        Self {
            A: None,
            phi,
            b,
            _private: (),
        }
    }

    pub fn new_dummy<Rng: rand::Rng + ?Sized>(r: usize, n: usize, rng: &mut Rng) -> Self {
        Self::new(
            Matrix::<R>::rand_symmetric(r, rng),
            vec![Vector::<R>::rand(n, rng); r],
            R::zero(),
        )
    }

    pub fn new_empty(r: usize, n: usize) -> Self {
        Self::new(
            Matrix::<R>::zeros(r, r),
            vec![Vector::<R>::zeros(n); r],
            R::zero(),
        )
    }

    pub fn is_valid_witness(&self, witness: &Witness<R>) -> bool {
        let inner_prods = inner_products(&witness.s);

        let mut res = R::zero();
        if let Some(A) = &self.A {
            let r = A.nrows();
            for i in 0..r {
                for j in 0..r {
                    res += A[(i, j)] * inner_prods[(i, j)];
                }
            }
        }

        for i in 0..self.phi.len() {
            res += self.phi[i].dot(&witness.s[i]);
        }

        res == self.b
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ConstantQuadDotProdFunction<R: PolyRing> {
    pub A: Option<Matrix<R>>,
    pub phi: Vec<Vector<R>>,
    pub b: R::BaseRing,
    _private: (), // Forbid direct initialization, force users to use new(), which does some basis debug_asserts
}

impl<R: PolyRing> ConstantQuadDotProdFunction<R> {
    pub fn new(A: Matrix<R>, phi: Vec<Vector<R>>, b: R::BaseRing) -> Self {
        let (r, n) = (A.nrows(), phi[0].len());
        debug_assert_eq!(A.ncols(), r, "A should be square");
        debug_assert_eq!(A.transpose(), A, "A should be symmetric");

        debug_assert_eq!(
            phi.len(),
            r,
            "phi should have the same length as the dimensions of A"
        );
        debug_assert!(
            phi.iter().all(|phi_i| phi_i.len() == n),
            "each phi_i should have the same length"
        );
        Self {
            A: Some(A),
            phi,
            b,
            _private: (),
        }
    }

    pub fn new_linear(phi: Vec<Vector<R>>, b: R::BaseRing) -> Self {
        let n = phi[0].len();
        debug_assert!(
            phi.iter().all(|phi_i| phi_i.len() == n),
            "each phi_i should have the same length"
        );
        Self {
            A: None,
            phi,
            b,
            _private: (),
        }
    }

    pub fn new_dummy<Rng: rand::Rng + ?Sized>(r: usize, n: usize, rng: &mut Rng) -> Self {
        Self::new(
            Matrix::<R>::rand_symmetric(r, rng),
            vec![Vector::<R>::rand(n, rng); r],
            R::BaseRing::zero(),
        )
    }

    pub fn is_valid_witness(&self, witness: &Witness<R>) -> bool {
        let inner_prods = inner_products(&witness.s);

        let mut res = R::zero();
        if let Some(A) = &self.A {
            let r = A.nrows();
            for i in 0..r {
                for j in 0..i + 1 {
                    res += A[(i, j)] * inner_prods[(i, j)];
                }
                for j in i + 1..r {
                    res += A[(i, j)] * inner_prods[(j, i)];
                }
            }
        }

        for i in 0..self.phi.len() {
            res += self.phi[i].dot(&witness.s[i]);
        }

        res.coeffs()[0] == self.b
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Witness<R: PolyRing> {
    pub s: Vec<Vector<R>>,
}

impl<R: PolyRing> Witness<R> {
    pub fn new_dummy(rank: usize, multiplicity: usize, _norm_bound: u64) -> Self {
        Self {
            // Guaranteed to have low norm, but not very interesting
            s: vec![Vector::<R>::from_fn(multiplicity, |_, _| R::zero()); rank],
        }
    }
}

impl<R: PolyRing> Instance<R> {
    pub fn is_valid_witness(&self, witness: &Witness<R>) -> bool {
        self.quad_dot_prod_funcs
            .iter()
            .all(|c| c.is_valid_witness(witness))
            && self
                .ct_quad_dot_prod_funcs
                .iter()
                .all(|c| c.is_valid_witness(witness))
    }
}

pub struct Size {
    pub num_witnesses: usize,
    pub witness_size: usize,
    pub norm_bound_sq: f64,
    pub num_constraints: usize,
    pub num_constant_constraints: usize,
}

pub struct PrincipalRelation<R: PolyRing> {
    _marker: std::marker::PhantomData<R>,
}

impl<R: PolyRing> Relation for PrincipalRelation<R> {
    type Size = Size;
    type Index = Index<R>;
    type Instance = Instance<R>;
    type Witness = Witness<R>;

    fn is_well_defined(pp: &Self::Index, x: &Self::Instance, w: Option<&Self::Witness>) -> bool {
        pp.is_wellformed_instance(x)
            && match w {
                Some(w) => pp.is_wellformed_witness(w),
                None => true,
            }
    }

    fn is_satisfied(pp: &Self::Index, x: &Self::Instance, w: &Self::Witness) -> bool {
        Self::is_well_defined(pp, x, Some(w)) && x.is_valid_witness(w)
    }

    fn generate_satisfied_instance(
        size: &Self::Size,
    ) -> (Self::Index, Self::Instance, Self::Witness) {
        assert!(size.num_witnesses > 0, "Need at least one witness");
        assert!(size.witness_size > 0, "Need positive witness size");
        assert!(size.num_constraints > 0 || size.num_constant_constraints > 0, "Need at least one constraint");

        // TODO: implement a more interesting satisfied relation
        let index = Index::<R>::new(
            size.num_witnesses,
            size.witness_size,
            size.norm_bound_sq,
            size.num_constraints,
            size.num_constant_constraints,
            &mut thread_rng(),
        );
        let instance = Instance::<R> {
            quad_dot_prod_funcs: (0..index.num_constraints)
                .map(|_| {
                    QuadDotProdFunction::<R>::new_dummy(
                        size.num_witnesses,
                        size.witness_size,
                        &mut thread_rng(),
                    )
                })
                .collect(),
            ct_quad_dot_prod_funcs: (0..index.num_constant_constraints)
                .map(|_| {
                    ConstantQuadDotProdFunction::<R>::new_dummy(
                        size.num_witnesses,
                        size.witness_size,
                        &mut thread_rng(),
                    )
                })
                .collect(),
        };

        let witness = Witness::new_dummy(size.num_witnesses, size.witness_size, 0);
        (index, instance, witness)
    }

    fn generate_unsatisfied_instance(
        size: &Self::Size,
    ) -> (Self::Index, Self::Instance, Self::Witness) {
        assert!(size.num_witnesses > 0, "Need at least one witness");
        assert!(size.witness_size > 0, "Need positive witness size");
        assert!(size.num_constraints > 0 || size.num_constant_constraints > 0, "Need at least one constraint");

        // TODO: implement a more interesting satisfied relation
        let index = Index::<R>::new(
            size.num_witnesses,
            size.witness_size,
            size.norm_bound_sq,
            size.num_constraints,
            size.num_constant_constraints,
            &mut thread_rng(),
        );
        let mut instance = Instance::<R> {
            quad_dot_prod_funcs: (0..index.num_constraints)
                .map(|_| {
                    QuadDotProdFunction::<R>::new_dummy(
                        size.num_witnesses,
                        size.witness_size,
                        &mut thread_rng(),
                    )
                })
                .collect(),
            ct_quad_dot_prod_funcs: (0..index.num_constant_constraints)
                .map(|_| {
                    ConstantQuadDotProdFunction::<R>::new_dummy(
                        size.num_witnesses,
                        size.witness_size,
                        &mut thread_rng(),
                    )
                })
                .collect(),
        };
        
        // Add single unsatisfied constraint
        instance.ct_quad_dot_prod_funcs.last_mut().unwrap().b = R::BaseRing::one();

        let witness = Witness::new_dummy(size.num_witnesses, size.witness_size, 0);
        (index, instance, witness)
    }
}

#[cfg(test)]
mod test {
    use lattirust_arithmetic::ntt::ntt_modulus;
    use lattirust_arithmetic::ring::Pow2CyclotomicPolyRingNTT;

    use crate::{test_generate_satisfied_instance, test_generate_unsatisfied_instance};
    use crate::principal_relation::{Instance, PrincipalRelation, Size};
    use crate::Relation;

    const Q: u64 = ntt_modulus::<64>(32);
    const D: usize = 64;

    type R = Pow2CyclotomicPolyRingNTT<Q, D>;
    type RELATION = PrincipalRelation<R>;

    const TEST_SIZE: Size = Size {
        num_witnesses: 10,
        witness_size: 128,
        norm_bound_sq: (Q as f64) / 10.,
        num_constraints: 32,
        num_constant_constraints: 64,
    };

    test_generate_satisfied_instance!(RELATION, TEST_SIZE);

    test_generate_unsatisfied_instance!(RELATION, TEST_SIZE);
}
