#![allow(non_snake_case)]

use anyhow::Error;
use ark_ff::{Fp, MontBackend};
use ark_ff::fields::MontConfig;
use ark_relations::r1cs::ConstraintSystem;
use num_traits::{One, Zero};

use crate::labrador::binary_r1cs::util::*;
use crate::labrador::prover::prove_principal_relation;
use crate::labrador::setup::CommonReferenceString;
use crate::labrador::util::concat;
use crate::labrador::witness::Witness;
use crate::lattice_arithmetic::challenge_set::labrador_challenge_set::LabradorChallengeSet;
use crate::lattice_arithmetic::challenge_set::weighted_ternary::WeightedTernaryChallengeSet;
use crate::lattice_arithmetic::matrix::{Matrix, sample_uniform_mat, Vector};
use crate::lattice_arithmetic::poly_ring::PolyRing;
use crate::lattice_arithmetic::ring::{Ring, Zq};
use crate::lattice_arithmetic::traits::{FromRandomBytes, Modulus};
use crate::nimue::arthur::LatticeArthur;
use crate::relations::labrador::principal_relation::{ConstantQuadDotProdFunction, PrincipalRelation, QuadDotProdFunction};

#[derive(MontConfig)]
#[modulus = "2"]
#[generator = "1"]
pub struct F2Config;

pub type F2 = Fp<MontBackend<F2Config, 1>, 1>;
pub type Z2 = Zq<2>;

pub struct BinaryR1CSCRS<R: PolyRing> {
    pub A: Matrix<R>,
    pub num_constraints: usize,
    pub num_variables: usize,
    pub m: usize,
}


impl<R: PolyRing> BinaryR1CSCRS<R> {
    pub fn new(num_constraints: usize, num_variables: usize) -> Self {
        let k = num_constraints;
        let n = num_variables;
        let d = R::dimension();
        // TODO: pad instead
        assert_eq!(k % d, 0, "number of constraints {} must be multiple of d = {}", k, d);
        assert_eq!(n % d, 0, "number of variables {} must be multiple of d = {}", n, d);
        let m = 10usize; // TODO: choose such that MSIS_{m, 2n+6k} is hard with l_inf bound = 1

        let q = R::BaseRing::modulus() as usize;
        assert!(n + 3 * k < q, "n + 3k = {} must be less than q = {} for soundness", n + 3 * k, q);
        assert!(6 * k < q, "6k = {} must be less than q = {} for soundness", 6 * k, q);
        assert!(128 * (n + 3 * k) < 15 * q, "n + 3*k = {} must be less than 15q/128 = {} to be able to compose with Labrador-core", n + 3 * k, 15 * q / 128);

        Self {
            A: sample_uniform_mat(m.div_ceil(d), (3 * k + n).div_ceil(d)),
            num_constraints,
            num_variables,
            m,
        }
    }

    pub fn pr_crs(&self) -> CommonReferenceString<R> {
        let d = R::dimension();
        let r_pr: usize = 8;
        let n_pr = self.num_variables.div_ceil(d);
        let norm_bound = (R::BaseRing::modulus() as f64).sqrt();
        let k = 10usize; //TODO ensure MSIS is hard for this k
        let k1 = k; // TODO
        let k2 = k; // TODO

        let num_quad_constraints = self.m.div_ceil(d) + 3 * n_pr;
        let num_constant_quad_constraints = 4 + 1 + SECPARAM;

        let basis = R::BaseRing::from(3u128); // TODO

        CommonReferenceString::<R>::new(r_pr, n_pr, R::dimension(), norm_bound, k, k1, k2, num_quad_constraints, num_constant_quad_constraints, basis)
    }
}

pub struct R1CSInstance<R: Ring> {
    // TODO: use sparse matrices instead
    pub A: Matrix<R>,
    pub B: Matrix<R>,
    pub C: Matrix<R>,
}

pub type BinaryR1CSInstance = R1CSInstance<Z2>;

pub struct R1CSWitness<R: Ring> {
    pub w: Vector<R>,
}

pub type BinaryR1CSWitness = R1CSWitness<Z2>;

pub fn ark_prove_binary_r1cs<R: PolyRing>(cs: ConstraintSystem<F2>, crs: CommonReferenceString<R>, arthur: &mut LatticeArthur<R>) {
    let r1cs_mat = cs.to_matrices().unwrap();
    let w = Vector::<F2>::from_vec([cs.instance_assignment, cs.witness_assignment].concat());
    let k = cs.num_constraints;
    let n = cs.num_instance_variables + cs.num_witness_variables;

    // let A = to_sparse_matrix(&r1cs_mat.a);
    // let B = to_sparse_matrix(&r1cs_mat.b);
    // let C = to_sparse_matrix(&r1cs_mat.c);
}

pub fn prove_binary_r1cs<'a, R: PolyRing>(crs: &BinaryR1CSCRS<R>, arthur: &'a mut LatticeArthur<R>, instance: &BinaryR1CSInstance, witness: &BinaryR1CSWitness) -> Result<&'a [u8], Error>
    where LabradorChallengeSet<R>: FromRandomBytes<R>, WeightedTernaryChallengeSet<R>: FromRandomBytes<R>
{
    let (A, B, C) = (&instance.A, &instance.B, &instance.C);
    let w = &witness.w;
    let (k, n) = (A.nrows(), A.ncols());
    let d = R::dimension();
    assert_eq!(k, n, "the current implementation only support k = n"); // TODO: remove this restriction by splitting a,b,c or w into multiple vectors

    // TODO: add statement

    let a = A * w;
    let b = B * w;
    let c = C * w;

    let a_R = lift::<R>(&a);
    let b_R = lift::<R>(&b);
    let c_R = lift::<R>(&c);
    let w_R = lift::<R>(w);

    let v = concat(&[&a_R, &b_R, &c_R, &w_R]);
    let t = &crs.A * &v;

    arthur.absorb_vector(&t).unwrap();

    let alpha = arthur.squeeze_binary_matrix(SECPARAM, k).unwrap();
    let beta = arthur.squeeze_binary_matrix(SECPARAM, k).unwrap();
    let gamma = arthur.squeeze_binary_matrix(SECPARAM, k).unwrap();

    // delta_i is computed mod 2, i.e., over Z2
    let delta = &alpha * A + &beta * B + &gamma * C;

    // g_i is computed over Zq
    let (a_BR, b_BR, c_BR) = (a.map(embed::<R::BaseRing>), b.map(embed::<R::BaseRing>), c.map(embed::<R::BaseRing>));
    let (alpha_BR, beta_BR, gamma_BR) = (alpha.map(embed::<R::BaseRing>), beta.map(embed::<R::BaseRing>), gamma.map(embed::<R::BaseRing>));
    let (delta_BR, w_BR) = (delta.map(embed::<R::BaseRing>), w.map(embed::<R::BaseRing>));

    let g = &alpha_BR * &a_BR + &beta_BR * &b_BR + &gamma_BR * &c_BR - &delta_BR * &w_BR;
    arthur.absorb_vector(&g).unwrap();

    let a_tilde = R::sigma_vec(&a_R);
    let b_tilde = R::sigma_vec(&b_R);
    let c_tilde = R::sigma_vec(&c_R);
    let w_tilde = R::sigma_vec(&w_R);

    let transcript = BinaryR1CSTranscript { t, alpha, beta, gamma, g, delta };

    let instance_pr = reduce(&crs, &instance, &transcript);

    let witness_pr = Witness::<R> {
        s: vec![a_R, b_R, c_R, w_R, a_tilde, b_tilde, c_tilde, w_tilde] // see definition of indices above
    };

    arthur.ratchet()?;
    prove_principal_relation(arthur, &instance_pr, &witness_pr, &crs.pr_crs())
}

