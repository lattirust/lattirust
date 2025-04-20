#![allow(non_snake_case)]

use ark_std::rand;
use ark_std::rand::thread_rng;
use num_traits::{Float, ToPrimitive};
use std::cmp::{max, max_by};
use std::fmt::Debug;

use lattice_estimator::msis;
use lattice_estimator::msis::{msis_h_128_l2, MSIS};
use lattice_estimator::norms::Norm;
use serde::{Deserialize, Serialize};
use tracing::info;

use lattirust_arithmetic::challenge_set::labrador_challenge_set::LabradorChallengeSet;
use lattirust_arithmetic::linear_algebra::Matrix;
use lattirust_arithmetic::ring::PolyRing;
use relations::principal_relation::Size;

/// Common reference string for one round of the LaBRADOR protocol
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct CommonReferenceString<R: PolyRing> {
    pub sec_param: usize,
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
    /// Second-level commitment matrix for first decomposition basis, of size k1 x (t1 * r * k)
    pub B: Matrix<R>,
    /// Second-level commitment matrix for second decomposition basis, of size k2 x (t2 * ((r * (r+1))/2))
    pub C: Matrix<R>,
    /// Second-level commitment matrix for first decomposition basis, of size k1 x (t1 * ((r * (r+1))/2))
    pub D: Matrix<R>,
    /// Decomposition basis for z-vectors, roughly equal to `b1` and `b2`
    pub b: u128,
    /// Decomposition basis for first-level commitments (t-vectors), roughly equal to `b` and `b2`
    pub b1: u128,
    /// Decomposition basis for inner product terms (g), roughly equal to `b1` and `b2`
    pub b2: u128,
    /// A reference to the CRS for the next recursive round, or `None` if this is the CRS for the last round
    pub next_crs: Option<Box<CommonReferenceString<R>>>,
}

#[derive(Copy, Clone, Debug, PartialEq)]
pub struct FoldedSize {
    pub size: Size,
    pub size_z: usize, // Number of Rq-elements in the decomposed z-witness
    pub size_t: usize, // Number of Rq-elements in the decomposed t-witness
    pub size_g: usize, // Number of Rq-elements in the decomposed g-witness
    pub size_h: usize, // Number of Rq-elements in the decomposed h-witness
    pub nu: usize, 
    pub mu: usize,
}

impl FoldedSize {
    pub fn num_witnesses_z(&self) -> usize {
        self.size_z.div_ceil(self.size.witness_len)
    }

    pub fn num_witnesses_t_g_h(&self) -> usize {
        (self.size_t + self.size_g + self.size_h).div_ceil(self.size.witness_len)
    }
}

impl<R: PolyRing> CommonReferenceString<R> {
    pub fn floor_to_even(x: f64) -> u128 {
        let floor = x.floor() as u128;
        if floor % 2 == 0 {
            floor
        } else {
            floor - 1
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

    pub fn ceil_to_even(x: f64) -> u128 {
        let ceil = x.ceil() as u128;
        if ceil % 2 == 0 {
            ceil
        } else {
            ceil + 1
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

    pub fn new_for_size(size: Size) -> CommonReferenceString<R> {
        let mut rng = thread_rng();
        Self::new(
            size.num_witnesses,
            size.witness_len,
            size.norm_bound_sq,
            size.num_constraints,
            size.num_constant_constraints,
            &mut rng,
        )
    }

    pub fn new<Rng: rand::Rng + ?Sized>(
        r: usize,
        n: usize,
        mut beta_sq: f64,
        num_constraints: usize,
        num_constant_constraints: usize,
        rng: &mut Rng,
    ) -> CommonReferenceString<R> {
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
        // let tau = LabradorChallengeSet::<R>::VARIANCE_SUM_COEFFS;

        // ||z||_∞ <= r * ||c_i||_∞ * ||s_i||_∞  <= r * ||c_i||_∞ * ||s_i||_2
        // b^2 must be < ||z||_∞ to be able to decompose each entry in z into 2 parts using standard (non-balanced) decomposition
        let linf_norm_z =
            r as u128 * LabradorChallengeSet::<R>::LINF_NORM * beta_sq.sqrt().floor() as u128;
        let b = Self::ceil_to_even((linf_norm_z as f64).sqrt());
        // let b = Self::round_to_even((s * (12. * r as f64 * tau).sqrt()).sqrt());
        info!("  b={b} (main decomposition basis)");
        assert!(b > 1);

        // Set auxiliary decomposition bases
        let (t1, b1) = Self::t1_b1(b);
        info!("  b1={b1}, t1={t1} (first decomposition basis and decomposition length)");

        let (t2, b2) = Self::t2_b2(r, n, beta_sq, b);
        info!("  b2={b2}, t2={t2} (second decomposition basis and decomposition length)");

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
        // let k = msis_1.upper_bound_h();
        let k = msis::lattice_estimator::find_optimal_h_dynamic(&msis_1, norm_bound_1, 128).expect(format!("failed to find secure rank for {msis_1}. Are the witness vectors long enough in your system?").as_str());
        // let k = msis::find_optimal_h(&msis_1, 128).expect(format!("failed to find secure rank for {msis_1}. Are there enough constraints in your system?").as_str());
        // let k = msis_1.find_optimal_h_dynamic(norm_bound_1, SECPARAM).expect(format!("failed to find secure rank for {msis_1}. Are there enough constraints in your system?").as_str());
        msis_1 = msis_1.with_h(k).with_length_bound(norm_bound_1(k));
        info!(
            "  k={k} for the MSIS instance {msis_1}  gives {} bits of security",
            msis_1.security_level()
        ); // TODO: silently assume that this gives us enough security, which it will for any reasonable parameters
        info!("  Chose largest k={k} for the MSIS instance {msis_1}");
        assert!(k > 0);

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
        assert!(k > 1);

        // TODO: this only gives 125 bits of soundness error for SECPARAM = 128, how do we best document this?
        // TODO: all params should be bigger to account for the slack of sqrt(128/30) per recursion level

        let mut crs = Self {
            sec_param: 128,
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
            B: Matrix::<R>::rand(k1, t1 * r * k, rng),
            C: Matrix::<R>::rand(k2, t2 * ((r * (r + 1)) / 2), rng),
            D: Matrix::<R>::rand(k1, t1 * ((r * (r + 1)) / 2), rng),
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

    /// Generate optimal sizes for the next round of the LaBRADOR protocol
    pub(crate) fn next_size(&self) -> FoldedSize {
        // Generate instance for next iteration of the protocol
        let size_zi = self.n;
        let z_decomp_len = 2;
        let size_z = z_decomp_len * size_zi;
        let size_t = self.r * self.t1 * self.k;
        let size_g = (self.r * (self.r + 1)).div_ceil(2) * self.t2;
        let size_h = (self.r * (self.r + 1)).div_ceil(2) * self.t1;
        let size_t_g_h = size_t + size_g + size_h;
        info!("Choosing splitting parameters for next round: size_z = {size_z}, size_t = {size_t}, size_g = {size_g}, size_h = {size_h}");
        info!(
            "Incoming constraints: standard: {}, constant-coeff: {}",
            self.num_constraints, self.num_constant_constraints
        );

        // Find nu, mu such that m * nu ≈ mu * n
        // TODO: Can we find these more efficiently?
        // TODO: Also, this is optimizing for proof size, maybe we want to only consider power-of-two values for n_next?
        let mut min_diff = f64::infinity();
        let mut best_mu = 0;
        let mut best_nu = 0;
        for mu in 1..size_t_g_h {
            for nu in 1..size_z {
                let diff = (((size_t_g_h as f64) / mu as f64)
                    - ((size_z as f64) / ((z_decomp_len * nu) as f64)))
                    .abs();
                if diff < min_diff {
                    min_diff = diff;
                    best_mu = mu;
                    best_nu = nu;
                }
            }
        }
        info!("Best splitting parameters for next round: nu = {best_nu}, mu = {best_mu}");

        let n_next = max(
            size_z.div_ceil(z_decomp_len * best_nu),
            size_t_g_h.div_ceil(best_mu),
        );
        let r_next = z_decomp_len * best_nu + best_mu;
        info!("Setting r_next = {r_next}, n_next = {n_next}");

        let num_quad_constraints = self.k + // k constraints for Az = sum_{i in [r]} c_i t_i (one per row)
            1 + // 1 constraint for <z, z> = sum_{i, j in [r]} g_ij c_i c_j
            1 + // 1 constraint for sum_{i in [r]} <phi_i, z> c_i = sum_{i, j in [r]} h_ij c_i c_j
            1 + // 1 constraint for sum_{i,j in [r]} a_ij * g_ij + sum_{i in [r]} h_ii = b
            self.k1 + // k1 constraints for u_1
            self.k2; // k2 constraints for u_2

        let next_norm_bound_squared = CommonReferenceString::<R>::next_norm_bound_sq(
            self.r,
            self.n,
            self.norm_bound_squared,
            self.k,
            self.b,
        );

        let next_size = Size {
            num_witnesses: r_next,
            witness_len: n_next,
            norm_bound_sq: next_norm_bound_squared,
            num_constraints: num_quad_constraints,
            num_constant_constraints: 0,
        };

        info!(
            "Outgoing constraints: standard: {}, constant-coeff: {}",
            next_size.num_constraints, next_size.num_constant_constraints
        );

        let folded_size = FoldedSize {
            size: next_size,
            size_z,
            size_t,
            size_g,
            size_h,
            nu: best_nu,
            mu: best_mu,
        };

        // There might be several empty witness vectors for padding (for example, if n_next is small)
        debug_assert!(best_nu * z_decomp_len >= folded_size.num_witnesses_z()); 
        debug_assert!(best_mu >= folded_size.num_witnesses_t_g_h());

        folded_size
    }

    fn next_crs(&self) -> Option<CommonReferenceString<R>> {
        if self.recurse() {
            return None;
        }
        let next_size = self.next_size();

        Some(CommonReferenceString::<R>::new_for_size(next_size.size))
    }
}
