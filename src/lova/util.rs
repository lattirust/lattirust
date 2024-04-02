use std::marker::PhantomData;
use ark_ff::Field;
use nalgebra::Scalar;
use serde::Serialize;

use lattice_estimator::norms::Norm::L2;
use lattice_estimator::sis::SIS;

use crate::lattice_arithmetic::matrix::{Matrix, sample_uniform_mat, SymmetricMatrix};
use crate::lattice_arithmetic::poly_ring::{ConvertibleField, UnsignedRepresentative};
use crate::lattice_arithmetic::traits::WithL2Norm;
use crate::relations::traits::Relation;

pub const SECPARAM: usize = 128;

impl<F: ConvertibleField> PublicParameters<F> {
    pub fn new(n: usize, norm_bound: f64) -> Self {
        let sis = SIS::new(0, F::modulus(), norm_bound, n, L2);
        let h = sis.find_optimal_n(SECPARAM).unwrap();
        let commitment_mat = sample_uniform_mat(h, n, &mut rand::thread_rng());
        let decomposition_basis = norm_bound.sqrt().round_ties_even() as u128; // TODO
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
    pub(crate) commitment: Matrix<F>,
    pub(crate) inner_products: SymmetricMatrix<F>,
}

pub type Witness<F> = Matrix<F>;

pub struct BaseRelation<F: ConvertibleField> {
    _marker: PhantomData<F>
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
                && x.commitment.ncols() == w.ncols() * pp.decomposition_length
                && x.inner_products.size() == w.ncols() * pp.decomposition_length;
        }
        is_well_defined
    }

    fn is_satisfied(pp: &Self::PublicParameters, x: &Self::Instance, w: &Self::Witness) -> bool {
        let norm_bound_sq = pp.norm_bound * pp.norm_bound;
        Self::is_well_defined(pp, x, Some(w))
            && &pp.commitment_mat * w == x.commitment
            && Into::<SymmetricMatrix::<F>>::into(w.transpose() * w) == x.inner_products
            && x.inner_products.diag().into_iter().all(|inner_product_ii| Into::<UnsignedRepresentative>::into(inner_product_ii).0 as f64 <= norm_bound_sq)
    }
}

pub fn norm_l2_squared_columnwise<F: Scalar + WithL2Norm + Copy>(mat: &Matrix<F>) -> Vec<u64> {
    mat.column_iter().map(|c| c.into_iter().map(|v| *v).collect::<Vec<F>>().l2_norm_squared()).collect()
}