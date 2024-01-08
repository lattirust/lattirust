#![allow(non_snake_case)]

use ark_std::iterable::Iterable;
use rayon::prelude::*;

use crate::labrador::setup::CommonReferenceString;
use crate::labrador::util::*;
pub use crate::labrador::witness::Witness;
use crate::lattice_arithmetic::balanced_decomposition::{decompose_balanced_polyring, decompose_balanced_vec};
use crate::lattice_arithmetic::challenge_set::labrador_challenge_set::LabradorChallengeSet;
use crate::lattice_arithmetic::challenge_set::weighted_ternary::WeightedTernaryChallengeSet;
use crate::lattice_arithmetic::matrix::{Matrix, norm_sq_vec, Vector};
use crate::lattice_arithmetic::poly_ring::PolyRing;
use crate::lattice_arithmetic::ring::Ring;
use crate::lattice_arithmetic::traits::{FromRandomBytes, WithLog2};
use crate::nimue::arthur::LatticeArthur;
use crate::relations::labrador::principal_relation::PrincipalRelation;

struct ProverState<R: PolyRing> {
    instance: PrincipalRelation<R>,
    witness: Witness<R>,
    crs: CommonReferenceString<R>,
    t: Option<Vec<Vector<R>>>,
    G: Option<Vec<Vec<R>>>,
    H: Option<Vec<Vec<R>>>,
    Pi: Option<Vec<Matrix<R>>>,
    psi: Option<Vec<Vector<R::BaseRing>>>,
    omega: Option<Vec<Vector<R::BaseRing>>>,
    alpha: Option<Vector<R>>,
    beta: Option<Vector<R>>,
    phi__: Option<Vec<Vec<Vector<R>>>>,
    c: Option<Vec<R>>,
}

// TODO: implement as Mul trait for Vector<R> once Vector is a struct a not a newtype anymore
fn mul_scalar_vector<R: Ring>(s: R, A: &Vector<R>) -> Vector<R> {
    A.map(|a_ij| s * a_ij)
}

// TODO: implement as Mul trait for Vector<R> once Vector is a struct a not a newtype anymore
fn mul_basescalar_vector<R: PolyRing>(s: R::BaseRing, A: &Vector<R>) -> Vector<R> {
    A.map(|a_ij| a_ij * s)
}

impl<R: PolyRing> ProverState<R> {
    fn new(instance: PrincipalRelation<R>, witness: Witness<R>, crs: CommonReferenceString<R>) -> Self {
        Self {
            instance,
            witness,
            crs,
            t: None,
            G: None,
            H: None,
            Pi: None,
            psi: None,
            omega: None,
            alpha: None,
            beta: None,
            phi__: None,
            c: None,
        }
    }
}


fn prove_1<R: PolyRing>(state: &mut ProverState<R>) -> Vector<R> {
    let (_, witness, crs) = (&state.instance, &state.witness, &state.crs);
    let t: Vec<Vector<R>> = witness.s.par_iter().map(
        |s_i| commit(&crs.A, s_i)
    ).collect();
    let t_decomp: Vec<Vec<Vector<R>>> = t.par_iter().map(|t_i| decompose_balanced_vec(t_i, crs.decomposition_basis, Some(crs.t1))).collect();

    let G = inner_products(&witness.s);

    let G_decomp: Vec<Vec<Vec<R>>> = G.par_iter().map(
        |G_i| G_i.par_iter().map(
            |G_ij| decompose_balanced_polyring(G_ij, crs.decomposition_basis, Some(crs.t2))
        ).collect()
    ).collect();

    let mut u1 = Vector::<R>::zeros(crs.k1);
    for i in 0..crs.r {
        assert_eq!(t_decomp[i].len(), crs.t1, "decomposition of t has the wrong number of elements");
        for k in 0..crs.t1 {
            u1 += &crs.B[i][k] * &t_decomp[i][k];
        }

        for j in 0..i + 1 {
            assert_eq!(G_decomp[i][j].len(), crs.t2, "decomposition of G has the wrong number of elements");
            for k in 0..crs.t2 {
                u1 += &crs.C[i][j][k] * G_decomp[i][j][k];
            }
        }
    }

    state.t.replace(t);
    state.G.replace(G);
    u1
}


fn prove_2<R: PolyRing>(state: &ProverState<R>) -> Vector<R::BaseRing> {
    let (_, witness, crs) = (&state.instance, &state.witness, &state.crs);
    let Pi = state.Pi.as_ref().expect("Pi not available");
    let mut p = Vector::<R::BaseRing>::zeros(256);
    for j in 0..256 {
        for i in 0..crs.r {
            let pi_ij = Pi[i].row(j).into_owned().transpose();
            p[j] += R::flattened(&pi_ij).dot(&R::flattened(&witness.s[i]));
        }
    }
    p
}


fn prove_3<R: PolyRing>(state: &mut ProverState<R>) -> Vec<R> {
    let (instance, witness, crs) = (&state.instance, &state.witness, &state.crs);
    let Pi = state.Pi.as_ref().expect("Pi not available");
    let psi = state.psi.as_ref().expect("psi not available");
    let omega = state.omega.as_ref().expect("omega not available");
    let num_aggregs = (128. / R::BaseRing::log2_q()).ceil() as usize;
    debug_assert_eq!(Pi.len(), crs.r);
    debug_assert_eq!(psi.len(), num_aggregs);
    debug_assert_eq!(omega.len(), num_aggregs);
    for Pi_i in Pi {
        debug_assert_eq!(Pi_i.nrows(), 256);
        debug_assert_eq!(Pi_i.ncols(), crs.n);
    }
    for psi_i in psi { debug_assert_eq!(psi_i.len(), instance.ct_quad_dot_prod_funcs.len()); }
    for omega_i in omega { debug_assert_eq!(omega_i.nrows(), 256); }


    let mut phi__ = vec![vec![Vector::<R>::zeros(crs.n); crs.r]; num_aggregs];
    let mut b = vec![R::zero(); num_aggregs];

    for k in 0..num_aggregs {
        for i in 0..crs.r {
            for j in 0..crs.r {
                // Compute a''_{ij}^{(k)}
                let mut a_ijk = R::zero();
                for l in 0..instance.ct_quad_dot_prod_funcs.len() {
                    a_ijk += instance.ct_quad_dot_prod_funcs[l].A[(i, j)] * psi[k][l];
                }

                b[k] += a_ijk * inner_prod(&witness.s[i], &witness.s[j]);
            }
            // Compute vec{phi}''_i^{(k)}
            for l in 0..instance.ct_quad_dot_prod_funcs.len() {
                phi__[k][i] += mul_basescalar_vector(psi[k][l], &instance.ct_quad_dot_prod_funcs[l].phi[i]);
            }
            for j in 0..256 {
                phi__[k][i] += mul_basescalar_vector(omega[k][j], &R::sigma_vec(&Pi[i].row(j).transpose()));
            }
            b[k] += inner_prod(&phi__[k][i], &witness.s[i]);
        }
    }
    state.phi__.replace(phi__);
    b
}

fn prove_4<R: PolyRing>(state: &mut ProverState<R>) -> Vector<R> {
    let (instance, witness, crs) = (&state.instance, &state.witness, &state.crs);
    let alpha = state.alpha.as_ref().expect("alpha not available");
    let beta = state.beta.as_ref().expect("beta not available");
    let phi__ = state.phi__.as_ref().expect("phi'' not available");

    // Compute phi (in parallel)
    let phi = (0..crs.r).into_par_iter().map(|i| {
        let mut phi_i = Vector::<R>::zeros(crs.n);
        for k in 0..instance.ct_quad_dot_prod_funcs.len() {
            phi_i += mul_scalar_vector(alpha[k], &instance.ct_quad_dot_prod_funcs[k].phi[i]);
        }
        for k in 0..instance.quad_dot_prod_funcs.len() {
            phi_i += mul_scalar_vector(beta[k], &phi__[k][i]);
        }
        phi_i
    }).collect::<Vec<_>>();

    // Compute H (in parallel)
    let mut H = inner_products2(&phi, &witness.s);
    let H_2 = inner_products2(&witness.s, &phi);
    for i in 0..crs.r {
        for j in 0..i + 1 {
            H[i][j] += H_2[i][j]; // TODO: divide by 2
        }
    }

    let mut u_2 = Vector::<R>::zeros(crs.k2);
    for i in 0..crs.r {
        for j in 0..i + 1 {
            let h_ij_decomp = decompose_balanced_polyring(&H[i][j], crs.decomposition_basis, Some(crs.t1));
            for k in 0..crs.t1 {
                u_2 += &crs.D[i][j][k] * h_ij_decomp[k];
            }
        }
    }

    state.H.replace(H);
    u_2
}


fn prove_5<R: PolyRing>(state: &ProverState<R>) -> (Vector<R>, Vec<Vector<R>>, Vec<Vec<R>>, Vec<Vec<R>>) {
    let (_, witness, crs) = (&state.instance, &state.witness, &state.crs);
    let c = state.c.as_ref().expect("c not available");
    let mut z = Vector::<R>::zeros(crs.n);
    debug_assert_eq!(c.len(), crs.r);
    for i in 0..crs.r {
        z += mul_scalar_vector(c[i], &witness.s[i]);
    }
    (z, state.t.clone().unwrap(), state.G.clone().unwrap(), state.H.clone().unwrap())
}

pub fn prove_principal_relation<R: PolyRing>(arthur: &mut LatticeArthur<R>, instance: PrincipalRelation<R>, witness: Witness<R>, crs: CommonReferenceString<R>) -> Result<&[u8], anyhow::Error>
    where LabradorChallengeSet<R>: FromRandomBytes<R>, WeightedTernaryChallengeSet<R>: FromRandomBytes<R>, i128: From<<R as PolyRing>::BaseRing>
{
    // Check dimensions and norms
    debug_assert_eq!(witness.s.len(), crs.r);
    for s_i in &witness.s {
        debug_assert_eq!(s_i.len(), crs.n);
    }
    let sum_norm_sq: u64 = witness.s.iter().map(|s_i| norm_sq_vec(s_i)).sum();
    debug_assert!(sum_norm_sq as f64 <= instance.norm_bound * instance.norm_bound);

    // Initialize Fiat-Shamir transcript
    // TODO: add public statement
    arthur.ratchet()?;

    // Prove
    let num_constraints = instance.quad_dot_prod_funcs.len();
    let num_ct_constraints = instance.ct_quad_dot_prod_funcs.len();
    let mut state = ProverState::new(instance, witness, crs.clone());

    let u_1 = prove_1(&mut state);
    assert_eq!(u_1.len(), crs.k1);
    arthur.absorb_vector(&u_1).expect("error absorbing prover message 1");

    let pis = arthur.squeeze_matrices::<R, WeightedTernaryChallengeSet<R>>(256, crs.n, crs.r).expect("error squeezing verifier message 1");
    state.Pi.replace(pis);

    let p = prove_2(&state);
    arthur.absorb_vector::<R::BaseRing>(&p).expect("error absorbing prover message 2");

    let psi = arthur.squeeze_vectors::<R::BaseRing, R::BaseRing>(num_ct_constraints, crs.num_aggregs).expect("error squeezing verifier message 2 (psi)");
    let omega = arthur.squeeze_vectors::<R::BaseRing, R::BaseRing>(256, crs.num_aggregs).expect("error squeezing verifier message 2 (omega)");
    state.psi.replace(psi);
    state.omega.replace(omega);

    let b__ = prove_3(&mut state);
    arthur.absorb_vec(&b__).expect("error absorbing prover message 3");

    let alpha = arthur.squeeze_vector::<R, R>(num_constraints).expect("error squeezing verifier message 3 (alpha)");
    let beta = arthur.squeeze_vector::<R, R>(crs.num_aggregs).expect("error squeezing verifier message 3 (beta)");
    state.alpha.replace(alpha);
    state.beta.replace(beta);

    let u_2 = prove_4(&mut state);
    arthur.absorb_vector(&u_2).expect("error absorbing prover message 4");

    let c = arthur.squeeze_vec::<R, LabradorChallengeSet<R>>(crs.r).expect("error squeezing verifier message 4");
    state.c.replace(c);

    let (z, t, G, H) = prove_5(&state);
    arthur.absorb_vector(&z).expect("error absorbing prover message 5 (z)");
    arthur.absorb_vectors(&t).expect("error absorbing prover message 5 (t)");
    arthur.absorb::<Vec<Vec<R>>>(&G).expect("error absorbing prover message 5 (G)");
    arthur.absorb::<Vec<Vec<R>>>(&H).expect("error absorbing prover message 5 (H)");

    return Ok(arthur.transcript());
    // TODO: absorb last prover message as well?
}