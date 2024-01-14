#![allow(non_snake_case)]

use nimue::InvalidTag;
use rayon::prelude::*;

use crate::labrador::common_reference_string::CommonReferenceString;
use crate::labrador::util::*;
use crate::lattice_arithmetic::balanced_decomposition::{decompose_balanced_polyring, decompose_balanced_vec};
use crate::lattice_arithmetic::challenge_set::labrador_challenge_set::LabradorChallengeSet;
use crate::lattice_arithmetic::challenge_set::weighted_ternary::WeightedTernaryChallengeSet;
use crate::lattice_arithmetic::matrix::{Matrix, norm_sq_ringelem, norm_sq_vec, norm_vec_basering, Vector};
use crate::lattice_arithmetic::poly_ring::PolyRing;
use crate::lattice_arithmetic::ring::Ring;
use crate::lattice_arithmetic::traits::FromRandomBytes;
use crate::nimue::merlin::LatticeMerlin;
use crate::relations::labrador::principal_relation::PrincipalRelation;

#[macro_export]
macro_rules! check {
    ($ cond: expr) => {
        {
            if !($cond) {
                return Err(InvalidTag::from("invalid proof"));
            }
        }
    };
    ( $ cond : expr , $ ( $ arg : tt ) + ) => {
        {
            if !($cond) {
                return Err(InvalidTag::from(format!("invalid proof: {}", format!($($arg)+))));
            }
        }
    };
}

#[macro_export]
macro_rules! check_eq {
    ( $ a : expr , $ b : expr ) => { check!($a == $b) };
    ( $ a : expr , $ b : expr , $ ( $ arg : tt ) + ) => { check!($a == $b, $($arg)+) };
}

pub struct VerifierState<R: PolyRing> {
    instance: PrincipalRelation<R>,
    crs: CommonReferenceString<R>,
    u1: Vector<R>,
    Pi: Vec<Matrix<R>>,
    psi: Vec<Vector<R::BaseRing>>,
    omega: Vec<Vector<R::BaseRing>>,
    b__: Vec<R>,
    alpha: Vector<R>,
    beta: Vector<R>,
    u2: Vector<R>,
    c: Vec<R>,
    z: Vector<R>,
    t: Vec<Vector<R>>,
    t_decomp: Vec<Vec<Vector<R>>>,
    G: Vec<Vec<R>>,
    G_decomp: Vec<Vec<Vec<R>>>,
    H: Vec<Vec<R>>,
    H_decomp: Vec<Vec<Vec<R>>>,
}

// Compute sum_{i,j in [r]} A_ij * c_i * c_j
fn linear_combination_symmetric_matrix<R: Ring>(A: &Vec<Vec<R>>, c: &Vec<R>) -> R {
    let n = A.len();
    debug_assert_eq!(c.len(), n);
    let mut lc = R::zero();
    for i in 0..n {
        debug_assert_eq!(A[i].len(), i + 1);
        for j in 0..i + 1 {
            lc += A[i][j] * c[i] * c[j];
        }
        for j in i + 1..n {
            // for j >= i+1, get A_ij from A_ji, since A is stored in lower triangular representation
            lc += A[j][i] * c[i] * c[j];
        }
    }
    lc
}

pub fn verify_dot_product_constraints<R: PolyRing>(state: &VerifierState<R>) -> Result<(), InvalidTag> {
    let (instance, crs) = (&state.instance, &state.crs);
    let (n, r, num_aggregs, K) = (crs.n, crs.r, crs.num_aggregs, instance.quad_dot_prod_funcs.len());
    let (z, t, c) = (&state.z, &state.t, &state.c);

    // Check Az = c1 * t1 + ... + c_r * t_r
    let Az = &crs.A * z;
    let t_lc = c.into_iter().zip(t).map(|(c_i, t_i)| mul_scalar_vector(*c_i, &t_i)).sum();
    check_eq!(Az, t_lc);

    // Check <z, z> = sum_{i,j in [r]} g_ij * c_i * c_j
    let g_lc = linear_combination_symmetric_matrix(&state.G, &c);
    check_eq!(inner_prod(&z, &z), g_lc);

    // Check sum_{i in [r]} <phi_i, z> * c_i = sum_{i,j in [r]} h_ij * c_i * c_j
    // Compute phi_i
    for psi_i in &state.psi { debug_assert_eq!(psi_i.len(), instance.ct_quad_dot_prod_funcs.len()); }
    for omega_i in &state.omega { debug_assert_eq!(omega_i.nrows(), 256); }

    let mut phi__ = vec![vec![Vector::<R>::zeros(crs.n); crs.r]; num_aggregs];

    for k in 0..num_aggregs {
        for i in 0..crs.r {
            for j in 0..crs.r {
                // Compute a''_{ij}^{(k)}
                let mut a_ijk = R::zero();
                for l in 0..K {
                    a_ijk += instance.ct_quad_dot_prod_funcs[l].A[(i, j)] * state.psi[k][l];
                }
            }
            // Compute vec{phi}''_i^{(k)}
            for l in 0..K {
                phi__[k][i] += mul_basescalar_vector(state.psi[k][l], &instance.ct_quad_dot_prod_funcs[l].phi[i]);
            }
            for j in 0..256 {
                phi__[k][i] += mul_basescalar_vector(state.omega[k][j], &R::sigma_vec(&state.Pi[i].row(j).transpose()));
            }
        }
    }
    let mut phi = vec![Vector::<R>::zeros(n); r];
    let mut phi_z_lc = R::zero();
    for i in 0..r {
        for k in 0..K {
            phi[i] += mul_scalar_vector(state.alpha[k], &instance.quad_dot_prod_funcs[k].phi[i]);
        }
        for k in 0..num_aggregs {
            phi[i] += mul_scalar_vector(state.beta[k], &phi__[i][k]);
        }
        phi_z_lc += inner_prod(&phi[i], &z) * c[i];
    }
    let h_lc = linear_combination_symmetric_matrix(&state.H, &c);
    check_eq!(phi_z_lc, h_lc);

    // Check sum_{i, j in [r]} a_ij * g_ij + sum_{i in [r]} h_ii = b
    let mut b = R::zero();
    for k in 0..K {
        b += state.alpha[k] * instance.quad_dot_prod_funcs[k].b;
    }
    for k in 0..num_aggregs {
        b += state.beta[k] * state.b__[k];
    }
    let mut a_g_h_lc = R::zero();
    let mut A = Matrix::<R>::zeros(r, r);
    for k in 0..K {
        A += mul_scalar_matrix(state.alpha[k], &instance.quad_dot_prod_funcs[k].A);
    }
    for k in 0..num_aggregs {
        A += mul_scalar_matrix(state.beta[k], &instance.ct_quad_dot_prod_funcs[k].A);
    }

    for i in 0..r {
        for j in 0..i + 1 {
            a_g_h_lc += A[(i, j)] * state.G[i][j];
        }
        for j in i + 1..r {
            a_g_h_lc += A[(i, j)] * state.G[j][i];
        }
        a_g_h_lc += state.H[i][i];
    }
    check_eq!(a_g_h_lc, b);

    // Check u1
    let mut t_g_lc = Vector::<R>::zeros(crs.k1);
    for i in 0..crs.r {
        for k in 0..crs.t1 {
            t_g_lc += &crs.B[i][k] * &state.t_decomp[i][k];
        }
        for j in 0..i + 1 { // TODO: or should we flip indices of G, since we're working with symmetric matrices?
            for k in 0..crs.t2 {
                t_g_lc += &crs.C[i][j][k] * state.G_decomp[i][j][k];
            }
        }
    }
    check_eq!(state.u1, t_g_lc);

    // Check u2
    let mut h_lc = Vector::<R>::zeros(crs.k2);
    for i in 0..r {
        for j in 0..i + 1 {
            for k in 0..crs.t1 {
                h_lc += &crs.D[i][j][k] * state.H_decomp[i][j][k];
            }
        }
    }
    check_eq!(state.u2, h_lc);

    Ok(())
}


pub fn verify_principal_relation<R: PolyRing>(merlin: &mut LatticeMerlin, instance: PrincipalRelation<R>, crs: CommonReferenceString<R>) -> Result<(), InvalidTag>
    where LabradorChallengeSet<R>: FromRandomBytes<R>, WeightedTernaryChallengeSet<R>: FromRandomBytes<R>, i128: From<<R as PolyRing>::BaseRing>, i64: From<<R as PolyRing>::BaseRing>
{
    // Init Fiat-Shamir transcript
    // TODO: add public statement
    merlin.ratchet()?;

    let (n, r) = (crs.n, crs.r);
    let num_constraints = instance.quad_dot_prod_funcs.len();
    let num_ct_constraints = instance.ct_quad_dot_prod_funcs.len();

    let u1 = merlin.next_vector::<R>(crs.k1).expect("error extracting prover message 1 from transcript");

    let Pi = merlin.challenge_matrices::<R, WeightedTernaryChallengeSet<R>>(256, n, r).expect("error extracting verifier message 1 from transcript");

    let p = merlin.next_vector::<R::BaseRing>(256).expect("error extracting prover message 2 from transcript");
    let norm_p = norm_vec_basering::<R>(&p);
    let p_norm_bound = 128f64.sqrt() * instance.norm_bound;
    check!(norm_p <=p_norm_bound, "||p||_2 = {} is not <= sqrt(128)*beta = {}", norm_p, p_norm_bound);

    let psi = merlin.challenge_vectors::<R::BaseRing, R::BaseRing>(num_ct_constraints, crs.num_aggregs).expect("error extracting verifier message 2 (psi) from transcript");
    let omega = merlin.challenge_vectors::<R::BaseRing, R::BaseRing>(256, crs.num_aggregs).expect("error extracting verifier message 2 (omega) from transcript");

    let b__ = merlin.next_vec::<R>(crs.num_aggregs).expect("error extracting prover message 3 from transcript");

    for k in 0..crs.num_aggregs {
        let mut rhs_k = omega[k].dot(&p);
        for l in 0..num_ct_constraints {
            rhs_k += psi[k][l] * instance.ct_quad_dot_prod_funcs[l].b.coeffs()[0];
        }
        check_eq!(b__[k].coeffs()[0], rhs_k, "constant coeff of b''^(k) = {:?} is not equal to rhs = {:?}", b__[k].coeffs()[0], rhs_k);
    }

    let alpha = merlin.challenge_vector::<R, R>(num_constraints).expect("error extracting verifier message 3 (alpha) from transcript");
    let beta = merlin.challenge_vector::<R, R>(crs.num_aggregs).expect("error extracting verifier message 3 (beta) from transcript");

    let u2 = merlin.next_vector::<R>(crs.k2).expect("error extracting prover message 4 from transcript");

    let c = merlin.challenge_vec::<R, LabradorChallengeSet<R>>(crs.r).expect("error extracting verifier message 4 from transcript");

    let z = merlin.next_vector::<R>(crs.n).expect("error extracting prover message 5 (z) from transcript");
    let t = merlin.next_vectors::<R>(crs.k, crs.t1).expect("error extracting prover message 5 (t) from transcript");
    let G = merlin.next_symmetric_matrix::<R>(crs.r).expect("error extracting prover message 5 (G) from transcript");
    let H = merlin.next_symmetric_matrix::<R>(crs.r).expect("error extracting prover message 5 (H) from transcript");

    // Bulk of the verification checks
    // TODO: compute a__ijk, etc.

    // We don't need to check if G and H are symmetric, since we only send the upper triangular part of the matrices

    let mut sum_norm_sq = 0u64;

    // Decompose z
    let z_decomp = decompose_balanced_vec(&z, crs.decomposition_basis, Some(2usize));
    assert_eq!(z_decomp.len(), 2);
    check!(l_inf_norm_vec(&z_decomp[0]) * 2 <= i64::from(crs.decomposition_basis) as u64);

    for z_i in z_decomp.iter() {
        sum_norm_sq += norm_sq_vec(z_i);
    }

    // Decompose t
    let t_decomp: Vec<Vec<Vector<R>>> = t.par_iter().map(|t_i| decompose_balanced_vec(t_i, crs.decomposition_basis, Some(crs.t1))).collect();
    for t_i in t_decomp.iter() {
        assert_eq!(t_i.len(), crs.t1);
        for k in 0..crs.t1 - 1 {
            check!(l_inf_norm_vec(&t_i[k]) * 2 <= i64::from(crs.decomposition_basis) as u64);
            sum_norm_sq += norm_sq_vec(&t_i[k]);
        }
    }

    // Decompose G
    let G_decomp: Vec<Vec<Vec<R>>> = G.par_iter().map(
        |G_i| G_i.par_iter().map(
            |G_ij| decompose_balanced_polyring(G_ij, crs.decomposition_basis, Some(crs.t2))
        ).collect()
    ).collect();
    for G_i in G_decomp.iter() {
        assert_eq!(G_i.len(), crs.r);
        for j in 0..crs.r {
            assert_eq!(G_i[j].len(), crs.t2);
            for k in 0..crs.t2 - 1 {
                check!(l_inf_norm(&G_i[j][k]) * 2 <= i64::from(crs.decomposition_basis) as u64);
                sum_norm_sq += norm_sq_ringelem(&G_i[j][k]);
            }
            sum_norm_sq += norm_sq_ringelem(&G_i[j][crs.t2 - 1]);
        }
    }

    // Decompose H
    let H_decomp: Vec<Vec<Vec<R>>> = H.par_iter().map(
        |H_i| H_i.par_iter().map(
            |H_ij| decompose_balanced_polyring(H_ij, crs.decomposition_basis, Some(crs.t2))
        ).collect()
    ).collect();
    for H_i in H_decomp.iter() {
        assert_eq!(H_i.len(), crs.r);
        for j in 0..crs.r {
            assert_eq!(H_i[j].len(), crs.t2);
            for k in 0..crs.t2 - 1 {
                check!(l_inf_norm(&H_i[j][k]) * 2 <= i64::from(crs.decomposition_basis) as u64);
                sum_norm_sq += norm_sq_ringelem(&H_i[j][k]);
            }
            sum_norm_sq += norm_sq_ringelem(&H_i[j][crs.t2 - 1]);
        }
    }

    // Consolidated norm check
    let beta_prime = 0u64;
    // TODO
    check!(sum_norm_sq <= beta_prime * beta_prime);


    let state = VerifierState {
        instance,
        crs,
        u1,
        Pi,
        psi,
        omega,
        b__,
        alpha,
        beta,
        z,
        t,
        t_decomp,
        G,
        G_decomp,
        H,
        H_decomp,
        c,
        u2,
    };

    // Dot Product Checks
    let end_recursion = true;
    if end_recursion {
        verify_dot_product_constraints(&state)?
    } else {
        // TODO: recurse
    }

    Ok(())
}

fn l_inf_norm_vec<R: PolyRing>(v: &Vector<R>) -> u64
    where i64: From<R::BaseRing>
{
    R::flattened_coeffs(v).into_iter().map(|x| i64::from(x).abs() as u64).max().unwrap()
}

fn l_inf_norm<R: PolyRing>(v: &R) -> u64
    where i64: From<R::BaseRing>
{
    v.coeffs().into_iter().map(|x| i64::from(x).abs() as u64).max().unwrap()
}