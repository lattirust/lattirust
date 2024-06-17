use std::marker::PhantomData;

use ark_std::rand::thread_rng;
use ark_std::{One, UniformRand};
use derive_more::Display;
use log::{debug, info};
use nimue::{DuplexHash, IOPattern};
use num_bigint::BigUint;
use num_traits::cast::ToPrimitive;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use labrador::common_reference_string::floor_to_even;
use lattice_estimator::norms::Norm::L2;
use lattice_estimator::sis::SIS;
use lattirust_arithmetic::balanced_decomposition::balanced_decomposition_max_length;
use lattirust_arithmetic::challenge_set::ternary::{TernaryChallengeSet, Trit};
use lattirust_arithmetic::linear_algebra::inner_products::inner_products_mat;
use lattirust_arithmetic::linear_algebra::SymmetricMatrix;
use lattirust_arithmetic::linear_algebra::{Matrix, Scalar, Vector};
use lattirust_arithmetic::nimue::iopattern::{
    RatchetIOPattern, SerIOPattern, SqueezeFromRandomBytes,
};
use lattirust_arithmetic::ring::{ConvertibleRing, SignedRepresentative};
use lattirust_arithmetic::traits::WithL2Norm;
use relations::traits::Relation;

const COMPLETENESS_ERROR: usize = 128;

#[derive(Clone, Copy, Display, Debug, Serialize, Deserialize)]
pub enum OptimizationMode {
    OptimizeForSpeed,
    OptimizeForSize,
    OptimizeForSpeedWithCompletenessError,
}

#[derive(Clone, Debug, Serialize)]
pub struct PublicParameters<F: ConvertibleRing> {
    pub commitment_mat: Matrix<F>,
    pub norm_bound: f64,
    pub decomposition_basis: u128,
    pub decomposition_length: usize,
    pub mode: OptimizationMode,
    pub inner_security_parameter: usize,
    pub commitment_length: usize,
}

fn pretty_print(param: f64) -> String {
    format!("{param} = 2^{}", param.log2())
}
impl<F: ConvertibleRing> PublicParameters<F> {
    /// Return the lowest norm bound such that folding is norm-preserving for a given `basis`, witness length `n`, and inner security parameter `security_parameter` (t in the paper).
    fn norm_lower_bound(basis: u128, n: usize, security_parameter: usize) -> f64 {
        // beta >= 2 * security_parameter * k * norm_decomp = 2 * security_parameter * k * b/2 * sqrt(n)
        // => beta >= security_parameter * (log_b(beta) + 2) * b * sqrt(n)
        // => beta / log_b(beta) >= security_parameter * b * sqrt(n) + 2 * security_parameter * b * sqrt(n)
        // => beta / ln(beta) >= security_parameter * b/ln(b) * sqrt(n) + 2 * security_parameter * b * sqrt(n)
        // To enforce beta / ln(beta) >= security_parameter * b/ln(b) * sqrt(n), we set beta = e^-W_{-1}(-1 / (security_parameter * b/ln(b) * sqrt(n))),
        // to which we then add 2 * security_parameter * b * sqrt(n).
        let tmp = security_parameter as f64 * (n as f64).sqrt();
        let threshold = tmp * (basis as f64) / (basis as f64).ln();
        let lambert_val = lambert_w_min1(-1. / threshold);

        let norm_bound = (-lambert_val).exp();
        let norm_bound = norm_bound.ceil(); // Round up to remove any floating point errors from this numerical estimate.
        let ln_norm_bound = norm_bound.ln();
        assert!(
            norm_bound / ln_norm_bound
                >= threshold,
            "norm_bound / ln(norm_bound) = {norm_bound}/{ln_norm_bound} = {}  must be at least security_parameter * sqrt(n) * b/ln(b) = {}",
            norm_bound / ln_norm_bound,
            threshold
        );
        norm_bound + 2. * basis as f64 * security_parameter as f64 * (n as f64).sqrt()
    }

    /// Set all parameters for a given witness length `n`and optimization mode `mode`, but do not generate expensive values (e.g., the commitment matrix).
    fn set_parameters(
        instance_length: usize,
        mode: OptimizationMode,
        security_parameter: usize,
        log_fiat_shamir: usize,
    ) -> Self {
        /*
        For inputs of length n we need to fulfill two conditions:
        1. The commitment must be binding, i.e., SIS[h, w=n, q, beta, l2] is hard
        2. For perfect completeness, the columns of the folded witness must have norm at most norm_bound, i.e., 2 * security_parameter * decomposition_length * norm_decomp <= norm_bound, where norm_decomp = floor(decomposition_basis/2) * sqrt(n)
         */
        let inner_security_parameter = (((security_parameter + 1 + log_fiat_shamir) as f64
            + f64::log2(5.))
            / (f64::log2(3.) - 1.))
            .ceil() as usize;
        info!("Setting parameters for security parameter = {security_parameter}, Fiat-Shamir log loss = {log_fiat_shamir}, witness size n = {instance_length}, modulus  q = {}, mode = {:?}", F::modulus(), mode);
        info!(" Inner security parameter = {inner_security_parameter}");
        assert!(
            5f64 * (2f64 / 3f64).powf(inner_security_parameter as f64)
                <= 2f64.powf(-((security_parameter + 1) as f64))
        );
        let sis_security_parameter = security_parameter + 1 + log_fiat_shamir;
        info!(" SIS security parameter   = {sis_security_parameter}");

        let norm_bound: f64 = match mode {
            OptimizationMode::OptimizeForSpeed => {
                // TODO: Do binary search for basis in [2, sqrt(norm_bound)] to find the largest b such that there exists an h for which SIS[h, w=n, q, beta, l2] is security-parameter-hard.
                let target_decomposition_length = 4;
                let norm_bound_hi = ((target_decomposition_length * inner_security_parameter)
                    * (target_decomposition_length * inner_security_parameter)
                    * instance_length) as f64;
                let _norm_bound_lo =
                    Self::norm_lower_bound(2, instance_length, inner_security_parameter);
                let _basis_hi = floor_to_even(norm_bound_hi.sqrt());
                let _basis_lo = 2;

                norm_bound_hi
            }
            OptimizationMode::OptimizeForSpeedWithCompletenessError => {
                let norm_bound = (instance_length as f64)
                    * ((8f64 * inner_security_parameter as f64) / 3f64
                        + f64::sqrt(COMPLETENESS_ERROR as f64 * f64::ln(2.)) / 2f64)
                        .powf(2.);
                assert!(norm_bound < F::modulus().to_f64().unwrap());
                norm_bound
            }
            OptimizationMode::OptimizeForSize => {
                Self::norm_lower_bound(2, instance_length, inner_security_parameter)
            }
        };
        info!("  Setting norm_bound = {}", pretty_print(norm_bound));

        log::logger().flush();
        assert!(norm_bound < F::modulus().to_f64().unwrap());

        let decomposition_basis: u128 = match mode {
            OptimizationMode::OptimizeForSpeedWithCompletenessError => {
                let basis = floor_to_even(norm_bound.sqrt());
                // Ensure that the constant 4 used to choose norm_bound above is correct.
                debug_assert_eq!(
                    balanced_decomposition_max_length(basis, norm_bound.floor() as u128),
                    4,
                    "assumption on decomposition length used to set norm_bound is violated"
                );
                info!(
                    "  Setting decomposition_basis = floor_to_even(sqrt(norm_bound)) = {}",
                    pretty_print(basis as f64)
                );
                basis
            }
            OptimizationMode::OptimizeForSpeed => {
                // For speed, we want to set decomposition_basis to be as large as possible while ensuring perfect completeness.
                let basis = floor_to_even(norm_bound.sqrt());
                // Ensure that the constant 4 used to choose norm_bound above is correct.
                debug_assert_eq!(
                    balanced_decomposition_max_length(basis, norm_bound.floor() as u128),
                    4,
                    "assumption on decomposition length used to set norm_bound is violated"
                );
                info!(
                    "  Setting decomposition_basis = floor_to_even(sqrt(norm_bound)) = {}",
                    pretty_print(basis as f64)
                );
                basis
            }
            OptimizationMode::OptimizeForSize => {
                /*
                The proof size is composed as follows
                    commitment: h * 2*decomposition_length*security_parameter * log2(q) bits
                    inner products: (2*decomposition_length*security_parameter)^2 * floor(decomposition_basis/2)^2 bits
                The minimum of the sum of these (for decomposition_basis >= 2) is attained for decomposition_basis ≈ e, which we round down to 2.
                */
                let basis = 2;
                info!("  Setting decomposition_basis = 2");
                basis
            }
        };

        let decomposition_length: usize =
            balanced_decomposition_max_length(decomposition_basis, norm_bound as u128);
        info!(
            "    decomposition_length = {}",
            pretty_print(decomposition_length as f64)
        );
        log::logger().flush();

        let norm_decomp = Self::decomposed_norm_max(decomposition_basis, instance_length);
        info!(
            "    decomposed witnesses will have column norm <= {}",
            pretty_print(norm_decomp)
        );
        log::logger().flush();

        let folded_norm = (2 * security_parameter * decomposition_length) as f64 * norm_decomp;
        info!(
            "    folded witness will have column norm <= {}",
            pretty_print(folded_norm)
        );
        log::logger().flush();

        debug_assert!(folded_norm <= norm_bound);

        info!("  Setting SIS parameters for commitment...");
        let sis = SIS::new(0, F::modulus(), norm_bound, instance_length, L2);
        let commitment_length = sis.find_optimal_h(sis_security_parameter).unwrap();

        info!(
            "    using SIS parameters {} for commitment, achieving {} bits of security",
            sis.with_h(commitment_length),
            sis.with_h(commitment_length).security_level()
        );
        log::logger().flush();

        Self {
            commitment_mat: Matrix::<F>::zeros(0, 0),
            norm_bound,
            decomposition_basis,
            decomposition_length,
            mode,
            inner_security_parameter,
            commitment_length,
        }
    }

    pub fn new(
        n: usize,
        mode: OptimizationMode,
        security_parameter: usize,
        log_fiat_shamir: usize,
    ) -> Self {
        let mut pp = Self::set_parameters(n, mode, security_parameter, log_fiat_shamir);

        info!("  Generating commitment matrix...");
        let commitment_mat = Matrix::<F>::par_rand(pp.commitment_length, n);
        info!("    done");

        pp.commitment_mat = commitment_mat;
        pp
    }

    pub fn witness_len(&self) -> usize {
        self.commitment_mat.ncols()
    }

    pub fn decomposed_norm_max(decomposition_basis: u128, witness_len: usize) -> f64 {
        decomposition_basis.div_floor(2) as f64 * (witness_len as f64).sqrt()
    }

    fn signed_bits(bound: usize) -> usize {
        bound.ilog2() as usize + 1 + 1
    }

    pub fn proof_size_bytes_with_mode(
        n: usize,
        mode: OptimizationMode,
        security_parameter: usize,
        log_fiat_shamir: usize,
    ) -> usize {
        Self::set_parameters(n, mode, security_parameter, log_fiat_shamir).proof_size_bytes()
    }

    pub fn proof_size_bytes(&self) -> usize {
        let modulus_bits = (F::modulus() - BigUint::from(1u32)).bits() as usize;

        // number of bits needed to represent integers in the range [-norm_bound^2, norm_bound^2]
        let norm_bound_sq_bits =
            Self::signed_bits((self.norm_bound * self.norm_bound).floor() as usize);

        // number of bits needed to represent integers in the range [-(decomposition_basis/2)^2, (decomposition_basis/2)^2]
        let decomp_upper_bound = (self.decomposition_basis / 2) as usize;
        let decomp_basis_sq_bits = Self::signed_bits(decomp_upper_bound * decomp_upper_bound);

        // lambda * lambda matrix, with entries w_2,i^T * w_1,j in [-norm_bound^2, norm_bound^2]
        let merge_proof_size =
            self.inner_security_parameter * self.inner_security_parameter * norm_bound_sq_bits;

        // 2 * k * lambda symmetric matrix, with signed entries in [-(b/2)^2, (b/2)^2]
        let inner_prods_size = ((2 * self.decomposition_length * self.inner_security_parameter)
            * (2 * self.decomposition_length * self.inner_security_parameter))
            / 2
            * decomp_basis_sq_bits;

        // m x 2*k*lambda matrix, with entries <= F
        let commitment_size = self.commitment_mat.nrows()
            * (2 * self.decomposition_length * self.inner_security_parameter)
            * modulus_bits;

        // 2*k*lambda x lambda with entries in {-1, 0, 1}
        let challenge_size = (2 * self.decomposition_length * self.inner_security_parameter)
            * self.inner_security_parameter
            * Self::signed_bits(1);

        let proof_size_bits =
            merge_proof_size + inner_prods_size + commitment_size + challenge_size;
        (proof_size_bits + 7) / 8
    }
}

impl<F: ConvertibleRing> PublicParameters<F> {
    pub fn powers_of_basis(&self) -> Vec<F> {
        self.powers_of_basis_int()
            .iter()
            .map(|&x| F::from(x))
            .collect()
    }

    pub fn powers_of_basis_int(&self) -> Vec<SignedRepresentative> {
        let mut pows = Vec::<SignedRepresentative>::with_capacity(self.decomposition_length);
        pows.push(SignedRepresentative::one());
        let basis = SignedRepresentative(self.decomposition_basis as i128);
        assert!(
            F::modulus()
                .to_f64()
                .unwrap()
                .log(self.decomposition_basis as f64)
                > (self.decomposition_length - 1) as f64,
            "ring must be large enough to support basis^i for i = 0..decomposition_length"
        );
        for i in 1..self.decomposition_length {
            pows.push(pows[i - 1] * basis);
        }
        pows
    }
}

#[derive(Clone, Debug, Serialize, PartialEq)]
pub struct Instance<F: ConvertibleRing> {
    pub commitment: Matrix<F>,
    pub inner_products: SymmetricMatrix<SignedRepresentative>,
}

impl<F: ConvertibleRing> Instance<F> {
    pub fn new(pp: &PublicParameters<F>, w: &Witness<F>) -> Instance<F> {
        Instance {
            commitment: &pp.commitment_mat * w,
            inner_products: inner_products_mat(&to_integers(&w)),
        }
    }
}

pub type Witness<F> = Matrix<F>;

pub struct BaseRelation<F: ConvertibleRing> {
    _marker: PhantomData<F>,
}

impl<F: ConvertibleRing> Relation for BaseRelation<F> {
    type PublicParameters = PublicParameters<F>;
    type Instance = Instance<F>;
    type Witness = Witness<F>;

    fn is_well_defined(
        pp: &Self::PublicParameters,
        x: &Self::Instance,
        w: Option<&Self::Witness>,
    ) -> bool {
        let mut is_well_defined = pp.commitment_mat.nrows() == x.commitment.nrows();
        if let Some(w) = w {
            is_well_defined = is_well_defined
                && pp.commitment_mat.ncols() == w.nrows()
                && x.commitment.ncols() == w.ncols()
                && x.inner_products.size() == w.ncols();
        }
        is_well_defined
    }

    fn is_satisfied(pp: &Self::PublicParameters, x: &Self::Instance, w: &Self::Witness) -> bool {
        let norm_bound_sq = pp.norm_bound * pp.norm_bound;
        let w_int = to_integers(w);

        let commitments_match = &pp.commitment_mat * w == x.commitment; // mod q
        debug!("  commitments match:      {commitments_match} (pp.commitment_mat * w == x.commitment (mod q))");

        let inner_products_match = inner_products_mat(&w_int) == x.inner_products; // over the integers
        debug!("  inner products match:   {inner_products_match} (w_int.transpose() * w_int == x.inner_products (over the integers))");

        let inner_products_are_small = x
            .inner_products
            .diag()
            .into_iter()
            .all(|inner_product_ii| inner_product_ii.0 as f64 <= norm_bound_sq);
        debug!("  inner products bounded: {inner_products_are_small} (x.inner_products()[i][i] <= norm_bound_sq)");

        Self::is_well_defined(pp, x, Some(w))
            && commitments_match
            && inner_products_match
            && inner_products_are_small
    }
}

pub fn norm_l2_squared_columnwise<F: ConvertibleRing>(mat: &Matrix<F>) -> Vec<u128> {
    mat.column_iter().map(|c| c.l2_norm_squared()).collect()
}

pub fn norm_l2_columnwise<F: ConvertibleRing>(mat: &Matrix<F>) -> Vec<f64> {
    mat.column_iter().map(|c| c.l2_norm()).collect()
}

pub fn to_integers<F: ConvertibleRing>(mat: &Matrix<F>) -> Matrix<SignedRepresentative> {
    mat.map(|x| Into::<SignedRepresentative>::into(x))
}

pub trait LovaIOPattern
where
    Self: SerIOPattern + SqueezeFromRandomBytes + RatchetIOPattern,
{
    fn merge<F: ConvertibleRing>(self, pp: &PublicParameters<F>) -> Self {
        self.absorb_matrix::<SignedRepresentative>(
            pp.inner_security_parameter,
            pp.inner_security_parameter,
            "cross inner products",
        )
        .ratchet()
    }

    fn reduce<F: ConvertibleRing>(self, pp: &PublicParameters<F>) -> Self {
        self.absorb_matrix::<F>(
            pp.commitment_mat.nrows(),
            2 * pp.inner_security_parameter * pp.decomposition_length,
            "commitment",
        )
        .ratchet()
        .absorb_symmetric_matrix::<SignedRepresentative>(
            2 * pp.inner_security_parameter * pp.decomposition_length,
            "inner products",
        )
        .ratchet()
        .squeeze_matrix::<Trit, TernaryChallengeSet<Trit>>(
            2 * pp.inner_security_parameter * pp.decomposition_length,
            pp.inner_security_parameter,
            "challenge",
        )
        .ratchet()
    }

    fn fold<F: ConvertibleRing>(self, pp: &PublicParameters<F>) -> Self {
        self.merge(pp).reduce(pp)
    }
}

impl<H: DuplexHash<u8>> LovaIOPattern for IOPattern<H> {}

pub fn rand_matrix_with_bounded_column_norms<
    F: UniformRand + Scalar + From<SignedRepresentative> + Send + Sync,
>(
    nrows: usize,
    ncols: usize,
    norm_bound: i128,
) -> Matrix<F> {
    Matrix::<F>::from_columns(
        (0..ncols)
            .into_par_iter()
            .map(|_| {
                Vector::<F>::rand_vector_with_bounded_norm(nrows, norm_bound, &mut thread_rng())
            })
            .collect::<Vec<_>>()
            .as_slice(),
    )
}

/// For $z \in (-0.25, 0)$, computes an approximation of $W_{-1}(z)$.
/// We use a Newton approximation with starting point as described in Lóczi, Lajos. (2022). Guaranteed- and high-precision evaluation of the Lambert W function. Applied Mathematics and Computation. 433. 10.1016/j.amc.2022.127406.
pub fn lambert_w_min1(z: f64) -> f64 {
    assert!(-0.25 < z && z < 0.);
    let ln_min_z = (-z).ln();
    let mut w = ln_min_z - (-ln_min_z).ln();
    const NUM_STEPS: usize = 1 << 10;
    for _ in 0..NUM_STEPS {
        let exp_w = w.exp();
        w = w - (w * exp_w - z) / (exp_w + w * exp_w);
    }
    w
}

pub fn header() -> String {
    r" 
 █████                                   
░░███                                    
 ░███         ██████  █████ ███████████  
 ░███        ███░░███░░███ ░░███░░░░░███ 
 ░███       ░███ ░███ ░███  ░███ ███████ 
 ░███      █░███ ░███ ░░███ ███ ███░░███ 
 ███████████░░██████   ░░█████ ░░████████
░░░░░░░░░░░  ░░░░░░     ░░░░░   ░░░░░░░░ 

                                         "
    .to_string()
}
