use std::marker::PhantomData;

use ark_std::rand::Rng;
use ark_std::{rand, One, UniformRand};
use log::{debug, info};
use nimue::{DuplexHash, IOPattern};
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

#[derive(Clone, Debug, Serialize, Deserialize)]
pub enum OptimizationMode {
    OptimizeForSpeed,
    OptimizeForSize,
}

#[derive(Clone, Debug, Serialize)]
pub struct PublicParameters<F: ConvertibleRing> {
    pub commitment_mat: Matrix<F>,
    pub norm_bound: f64,
    pub decomposition_basis: u128,
    pub decomposition_length: usize,
    pub mode: OptimizationMode,
    pub security_parameter: usize,
}

fn pretty_print(param: f64) -> String {
    format!("{param} = 2^{}", param.log2())
}
impl<F: ConvertibleRing> PublicParameters<F> {
    pub fn new(n: usize, mode: OptimizationMode) -> Self {
        /*
        For inputs of length n we need to fulfill two conditions:
        1. The commitment must be binding, i.e., SIS[h, w=n, q, beta, l2] is hard
        2. For perfect completeness, the columns of the folded witness must have norm at most norm_bound, i.e., 2 * security_parameter * decomposition_length * norm_decomp <= norm_bound, where norm_decomp = floor(decomposition_basis/2) * sqrt(n)
         */
        let security_parameter = 128;
        info!("Setting parameters for security parameter = {security_parameter}, witness size n = {n}, modulus  q = {}, mode = {:?}", F::modulus(), mode);
        let norm_bound: f64 = match mode {
            OptimizationMode::OptimizeForSpeed => {
                let target_decomposition_length = 4;
                let norm_bound = ((target_decomposition_length * security_parameter)
                    * (target_decomposition_length * security_parameter)
                    * n) as f64;
                info!("  Setting norm_bound = ({target_decomposition_length} * security_parameter)^2 * n = {}", pretty_print(norm_bound));
                norm_bound
            }
            OptimizationMode::OptimizeForSize => {
                // beta >= 2 * security_parameter * k * norm_decomp = 2 * security_parameter * k * b/2 * sqrt(n)
                // => beta >= security_parameter * (log_b(beta) + 2) * b * sqrt(n)
                // => beta / log_b(beta) >= security_parameter * b * sqrt(n) + 2 * security_parameter * b * sqrt(n)
                // => beta / ln(beta) >= security_parameter * b/ln(b) * sqrt(n) + 2 * security_parameter * b * sqrt(n)
                // To enforce beta / ln(beta) >= security_parameter * b/ln(b) * sqrt(n), we set beta = e^-W_{-1}(-1 / (security_parameter * b/ln(b) * sqrt(n))),
                // to which we then add 2 * security_parameter * b * sqrt(n).
                let basis = 2;
                let tmp = security_parameter as f64 * (n as f64).sqrt();
                let threshold = tmp * (basis as f64) / (basis as f64).ln();
                let lambert_val = lambert_w_min1(-1. / threshold);

                let norm_bound = (-lambert_val).exp();
                let norm_bound = norm_bound.ceil(); // Round up to remove any floating point errors from this numerical estimate.
                info!("  Setting norm_bound = {}", pretty_print(norm_bound));

                let ln_norm_bound = norm_bound.ln();
                assert!(
                    norm_bound / ln_norm_bound
                        >= threshold,
                    "norm_bound / ln(norm_bound) = {norm_bound}/{ln_norm_bound} = {}  must be at least security_parameter * sqrt(n) * b/ln(b) = {}",
                    norm_bound / ln_norm_bound,
                    threshold
                );
                norm_bound
                    + 2. * basis as f64
                        * security_parameter as f64
                        * (basis as f64)
                        * (n as f64).sqrt()
            }
        };
        log::logger().flush();
        assert!((norm_bound as u128) < F::modulus());

        let decomposition_basis: u128 = match mode {
            OptimizationMode::OptimizeForSpeed => {
                // For speed, we want to set decomposition_basis to be as large as possible while ensuring perfect completeness.
                let basis = floor_to_even(norm_bound.sqrt());
                // Ensure that the constant 3 used to choose norm_bound above is correct.
                debug_assert_eq!(
                    balanced_decomposition_max_length(basis, norm_bound.floor() as u128),
                    4,
                    "assumption on decomposition length used to set norm_bound is violated"
                );
                info!("  Setting decomposition_basis = floor(sqrt(norm_bound) + 1) = {basis}");
                basis
            }
            OptimizationMode::OptimizeForSize => {
                /*
                The proof size is composed as follows
                    commitment: h * 2*decomposition_length*security_parameter * log2(q) bits
                    inner products: (2*decomposition_length*security_parameter)^2 * floor(decomposition_basis/2)^2 bits
                    challenge: 2*decomposition_length*security_parameter^2 * 2 bits
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
            "  decomposition_length = {}",
            pretty_print(decomposition_length as f64)
        );
        log::logger().flush();

        let norm_decomp = Self::decomposed_norm_max(decomposition_basis, n);
        info!(
            "  decomposed witnesses will have column norm <= {}",
            pretty_print(norm_decomp)
        );
        log::logger().flush();

        let folded_norm = (2 * security_parameter * decomposition_length) as f64 * norm_decomp;
        info!(
            "  folded witness will have column norm <= {}",
            pretty_print(folded_norm)
        );
        log::logger().flush();

        debug_assert!(folded_norm <= norm_bound);

        let sis = SIS::new(0, F::modulus(), norm_bound, n, L2);
        let h = sis.find_optimal_h(security_parameter).unwrap();
        info!(
            "  Using SIS parameters {} for commitment, achieving {} bits of security",
            sis.with_h(h),
            sis.with_h(h).security_level()
        );
        log::logger().flush();

        let commitment_mat = Matrix::<F>::rand(h, n, &mut rand::thread_rng());

        Self {
            commitment_mat,
            norm_bound,
            decomposition_basis,
            decomposition_length,
            mode,
            security_parameter,
        }
    }

    pub fn witness_len(&self) -> usize {
        self.commitment_mat.ncols()
    }

    pub fn decomposed_norm_max(decomposition_basis: u128, witness_len: usize) -> f64 {
        decomposition_basis.div_floor(2) as f64 * (witness_len as f64).sqrt()
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
        for i in 1..self.decomposition_length {
            pows.push(pows[i - 1] * basis);
            assert!(
                pows[i].0 < F::modulus() as i128,
                "basis^i = {basis}^{i} = {} must be less than modulus = {}",
                pows[i],
                F::modulus()
            );
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
        let inner_products_match =
            Into::<SymmetricMatrix<SignedRepresentative>>::into(w_int.transpose() * w_int)
                == x.inner_products; // over the integers
        let inner_products_are_small = x
            .inner_products
            .diag()
            .into_iter()
            .all(|inner_product_ii| inner_product_ii.0 as f64 <= norm_bound_sq);

        debug!("  commitments match:      {commitments_match} (pp.commitment_mat * w == x.commitment (mod q))");
        debug!("  inner products match:   {inner_products_match} (w_int.transpose() * w_int == x.inner_products (over the integers))");
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
            pp.security_parameter,
            pp.security_parameter,
            "cross inner products",
        )
        .ratchet()
    }

    fn reduce<F: ConvertibleRing>(self, pp: &PublicParameters<F>) -> Self {
        self.absorb_matrix::<F>(
            pp.commitment_mat.nrows(),
            2 * pp.security_parameter * pp.decomposition_length,
            "commitment",
        )
        .ratchet()
        .absorb_symmetric_matrix::<SignedRepresentative>(
            2 * pp.security_parameter * pp.decomposition_length,
            "inner products",
        )
        .ratchet()
        .squeeze_matrix::<Trit, TernaryChallengeSet<Trit>>(
            2 * pp.security_parameter * pp.decomposition_length,
            pp.security_parameter,
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
    F: UniformRand + Scalar + From<SignedRepresentative>,
>(
    nrows: usize,
    ncols: usize,
    rng: &mut impl Rng,
    norm_bound: i128,
) -> Matrix<F> {
    Matrix::<F>::from_columns(
        (0..ncols)
            .map(|_| Vector::<F>::rand_vector_with_bounded_norm(nrows, norm_bound, rng))
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
    const NUM_STEPS: usize = 1<<10;
    for _ in 0..NUM_STEPS {
        let exp_w = w.exp();
        w = w - (w * exp_w - z) / (exp_w + w * exp_w);
    }
    w
}
