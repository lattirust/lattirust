use std::marker::PhantomData;

use ark_ff::Field;
use nalgebra::Scalar;
use nimue::{DuplexHash, IOPattern};
use serde::Serialize;

use lattice_estimator::norms::Norm::L2;
use lattice_estimator::sis::SIS;

use crate::labrador::common_reference_string::round_to_even;
use crate::labrador::util::inner_products_mat;
use crate::lattice_arithmetic::challenge_set::ternary::{TernaryChallengeSet, Trit};
use crate::lattice_arithmetic::matrix::{Matrix, sample_uniform_mat, SymmetricMatrix};
use crate::lattice_arithmetic::poly_ring::{ConvertibleField, SignedRepresentative};
use crate::lattice_arithmetic::traits::WithL2Norm;
use crate::nimue::iopattern::{RatchetIOPattern, SerIOPattern, SqueezeFromRandomBytes};
use crate::relations::traits::Relation;

pub const SECPARAM: usize = 128;

impl<F: ConvertibleField> PublicParameters<F> {
    pub fn new(n: usize, norm_bound: f64) -> Self {
        let sis = SIS::new(0, F::modulus(), norm_bound, n, L2);
        let h = sis.find_optimal_n(SECPARAM).unwrap();
        let commitment_mat = sample_uniform_mat(h, n, &mut rand::thread_rng());
        let decomposition_basis = round_to_even(norm_bound.sqrt()); // TODO
        let decomposition_length = norm_bound.log(decomposition_basis as f64).ceil() as usize;
        Self { commitment_mat, norm_bound, decomposition_basis, decomposition_length }
    }
}

#[derive(Clone, Debug, Serialize)]
pub struct PublicParameters<F: ConvertibleField> {
    pub commitment_mat: Matrix<F>,
    pub norm_bound: f64,
    pub decomposition_basis: u128,
    pub decomposition_length: usize,
}

impl<F: ConvertibleField> PublicParameters<F> {
    pub fn powers_of_basis(&self) -> Vec<F> {
        let mut pows = Vec::<F>::with_capacity(self.decomposition_length);
        pows.push(F::one());
        let basis = F::from(self.decomposition_basis);
        for i in 1..self.decomposition_length {
            pows.push(pows[i - 1] * basis)
        }
        pows
    }
}

#[derive(Clone, Debug, Serialize, PartialEq)]
pub struct Instance<F: ConvertibleField> {
    pub commitment: Matrix<F>,
    pub inner_products: SymmetricMatrix<i128>,
}

impl<F: ConvertibleField> Instance<F> {
    pub fn new(pp: &PublicParameters<F>, w: &Witness<F>) -> Instance<F> {
        Instance {
            commitment: &pp.commitment_mat * w,
            inner_products: inner_products_mat(&to_integers(&w)).into(),
        }
    }
}

pub type Witness<F> = Matrix<F>;

pub struct BaseRelation<F: ConvertibleField> {
    _marker: PhantomData<F>,
}

impl<F: ConvertibleField> Relation for BaseRelation<F> {
    type PublicParameters = PublicParameters<F>;
    type Instance = Instance<F>;
    type Witness = Witness<F>;

    fn is_well_defined(pp: &Self::PublicParameters, x: &Self::Instance, w: Option<&Self::Witness>) -> bool {
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
        Self::is_well_defined(pp, x, Some(w))
            && &pp.commitment_mat * w == x.commitment
            && Into::<SymmetricMatrix<i128>>::into(w_int.transpose() * w_int) == x.inner_products
            && x.inner_products.diag().into_iter().all(|inner_product_ii| inner_product_ii as f64 <= norm_bound_sq)
    }
}

pub fn norm_l2_squared_columnwise<F: Scalar + WithL2Norm + Copy>(mat: &Matrix<F>) -> Vec<u64> {
    mat.column_iter().map(|c| c.into_iter().map(|v| *v).collect::<Vec<F>>().l2_norm_squared()).collect()
}

pub fn to_integers<F: ConvertibleField>(mat: &Matrix<F>) -> Matrix<i128> {
    mat.map(|x| Into::<SignedRepresentative>::into(x).0)
}

pub trait LovaIOPattern
    where Self: SerIOPattern + SqueezeFromRandomBytes + RatchetIOPattern
{
    fn folding_round<F: ConvertibleField>(self, pp: &PublicParameters<F>) -> Self {
        self.absorb_matrix::<F>(pp.commitment_mat.nrows(), 2 * SECPARAM * pp.decomposition_length, "commitment")
            .ratchet()
            .absorb_symmetric_matrix::<i128>(2 * SECPARAM * pp.decomposition_length, "inner products")
            .ratchet()
            .squeeze_matrix::<Trit, TernaryChallengeSet<Trit>>(2 * SECPARAM * pp.decomposition_length, SECPARAM, "challenge")
    }
}

impl<H: DuplexHash<u8>> LovaIOPattern for IOPattern<H> {}