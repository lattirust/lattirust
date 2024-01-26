#![allow(non_snake_case)]

use nimue::InvalidTag;
use rayon::prelude::*;

use crate::labrador::common_reference_string::CommonReferenceString;
use crate::labrador::shared::{BaseTranscript, fold_instance, next_norm_bound_sq, recurse};
use crate::labrador::util::*;
use crate::lattice_arithmetic::balanced_decomposition::{decompose_balanced_polyring, decompose_balanced_vec};
use crate::lattice_arithmetic::challenge_set::labrador_challenge_set::LabradorChallengeSet;
use crate::lattice_arithmetic::challenge_set::weighted_ternary::WeightedTernaryChallengeSet;
use crate::lattice_arithmetic::matrix::{Matrix, norm_sq_ringelem, norm_sq_vec, norm_sq_vec_basering, norm_vec_basering, Vector};
use crate::lattice_arithmetic::poly_ring::PolyRing;
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

/// Verify the final dot product constraints and consolidated norm check, used in the last step of the recursion
pub fn verify_final<R: PolyRing>(transcript: &BaseTranscript<R>) -> Result<(), InvalidTag> {
    let (instance, crs) = (&transcript.instance, &transcript.crs);
    let (n, r, num_aggregs, K) = (crs.n, crs.r, crs.num_aggregs, instance.quad_dot_prod_funcs.len());
    let (z, t, c) = (transcript.z.as_ref().expect("z not available"), transcript.t.as_ref().expect("t not available"), transcript.c.as_ref().expect("c not available"));
    let G = transcript.G.as_ref().expect("G not available");
    let H = transcript.H.as_ref().expect("H not available");

    // Decompose z, t, G, H
    let mut sum_norm_sq = 0u64;
    let decomposition_basis_u64 = Into::<i64>::into(crs.decomposition_basis) as u64;

    let z_decomp = decompose_balanced_vec(&z, crs.decomposition_basis, Some(2usize));
    assert_eq!(z_decomp.len(), 2);
    check!(l_inf_norm_vec(&z_decomp[0]) * 2 <= decomposition_basis_u64);

    for z_i in z_decomp.iter() {
        sum_norm_sq += norm_sq_vec(z_i);
    }

    let t_decomp: Vec<Vec<Vector<R>>> = t.par_iter().map(|t_i| decompose_balanced_vec(t_i, crs.decomposition_basis, Some(crs.t1))).collect();
    for t_i in t_decomp.iter() {
        assert_eq!(t_i.len(), crs.t1);
        for k in 0..crs.t1 - 1 {
            check!(l_inf_norm_vec(&t_i[k]) * 2 <= decomposition_basis_u64);
            sum_norm_sq += norm_sq_vec(&t_i[k]);
        }
    }

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
                check!(l_inf_norm(&G_i[j][k]) * 2 <= decomposition_basis_u64);
                sum_norm_sq += norm_sq_ringelem(&G_i[j][k]);
            }
            sum_norm_sq += norm_sq_ringelem(&G_i[j][crs.t2 - 1]);
        }
    }

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
                check!(l_inf_norm(&H_i[j][k]) * 2 <= decomposition_basis_u64);
                sum_norm_sq += norm_sq_ringelem(&H_i[j][k]);
            }
            sum_norm_sq += norm_sq_ringelem(&H_i[j][crs.t2 - 1]);
        }
    }

    // Consolidated norm check
    let beta_prime_sq = next_norm_bound_sq(&transcript);
    check!(sum_norm_sq as f64 <= beta_prime_sq);

    // Check Az = c1 * t1 + ... + c_r * t_r
    let Az = &crs.A * z;
    let t_lc = c.into_iter().zip(t).map(|(c_i, t_i)| t_i * *c_i).sum();
    check_eq!(Az, t_lc);

    // Check <z, z> = sum_{i,j in [r]} g_ij * c_i * c_j
    let G = transcript.G.as_ref().expect("G not available");
    let g_lc = linear_combination_symmetric_matrix(G, &c);
    check_eq!(inner_prod(&z, &z), g_lc);

    // Check sum_{i in [r]} <phi_i, z> * c_i = sum_{i,j in [r]} h_ij * c_i * c_j
    // Compute phi_i
    let psi = transcript.psi.as_ref().expect("psi not available");
    let omega = transcript.omega.as_ref().expect("omega not available");
    let alpha = transcript.alpha.as_ref().expect("alpha not available");
    let beta = transcript.beta.as_ref().expect("beta not available");
    let phi = transcript.phi.as_ref().expect("phi not available");

    for psi_i in psi { debug_assert_eq!(psi_i.len(), instance.ct_quad_dot_prod_funcs.len()); }
    for omega_i in omega { debug_assert_eq!(omega_i.nrows(), 256); }

    let mut phi_z_lc = R::zero();
    for i in 0..r {
        phi_z_lc += inner_prod(&phi[i], &z) * c[i];
    }
    let h_lc = linear_combination_symmetric_matrix(&H, &c);
    check_eq!(phi_z_lc, h_lc);

    // Check sum_{i, j in [r]} a_ij * g_ij + sum_{i in [r]} h_ii = b
    let b__ = transcript.b__.as_ref().expect("b'' not available");
    let mut b = R::zero();
    for k in 0..K {
        b += alpha[k] * instance.quad_dot_prod_funcs[k].b;
    }
    for k in 0..num_aggregs {
        b += beta[k] * b__[k];
    }
    let mut a_g_h_lc = R::zero();
    let mut A = Matrix::<R>::zeros(r, r);
    for k in 0..K {
        A += &instance.quad_dot_prod_funcs[k].A * alpha[k];
    }
    for k in 0..num_aggregs {
        A += &instance.ct_quad_dot_prod_funcs[k].A * beta[k];
    }

    for i in 0..r {
        for j in 0..i + 1 {
            a_g_h_lc += A[(i, j)] * G[i][j];
        }
        for j in i + 1..r {
            a_g_h_lc += A[(i, j)] * G[j][i];
        }
        a_g_h_lc += H[i][i];
    }
    check_eq!(a_g_h_lc, b);

    // Check u1
    let u_1 = transcript.u_1.as_ref().expect("u_1 not available");
    let mut t_g_lc = Vector::<R>::zeros(crs.k1);
    for i in 0..crs.r {
        for k in 0..crs.t1 {
            t_g_lc += &crs.B[i][k] * &t_decomp[i][k];
        }
        for j in 0..i + 1 { // TODO: or should we flip indices of G, since we're working with symmetric matrices?
            for k in 0..crs.t2 {
                t_g_lc += &crs.C[i][j][k] * G_decomp[i][j][k];
            }
        }
    }
    check_eq!(*u_1, t_g_lc);

    // Check u2
    let u_2 = transcript.u_2.as_ref().expect("u_2 not available");
    let mut h_lc = Vector::<R>::zeros(crs.k2);
    for i in 0..r {
        for j in 0..i + 1 {
            for k in 0..crs.t1 {
                h_lc += &crs.D[i][j][k] * H_decomp[i][j][k];
            }
        }
    }
    check_eq!(*u_2, h_lc);

    Ok(())
}

/// Verify consistency for one instance of the core Labrador protocol, used in each step of the recursion
pub fn verify_core<'a, R: PolyRing>(crs: &'a CommonReferenceString<R>, instance: &'a PrincipalRelation<R>, merlin: &mut LatticeMerlin) -> Result<BaseTranscript<'a, R>, InvalidTag>
    where LabradorChallengeSet<R>: FromRandomBytes<R>, WeightedTernaryChallengeSet<R>: FromRandomBytes<R>
{
    let (n, r) = (crs.n, crs.r);
    let num_constraints = instance.quad_dot_prod_funcs.len();
    let num_ct_constraints = instance.ct_quad_dot_prod_funcs.len();

    let u_1 = merlin.next_vector::<R>(crs.k1).expect("error extracting prover message 1 from transcript");

    let Pi = merlin.challenge_matrices::<R, WeightedTernaryChallengeSet<R>>(256, n, r).expect("error extracting verifier message 1 from transcript");

    let p = merlin.next_vector::<R::BaseRing>(256).expect("error extracting prover message 2 from transcript");
    let norm_p_sq = norm_sq_vec_basering::<R>(&p);
    let p_norm_bound_sq = 128f64 * crs.norm_bound_squared;
    check!(norm_p_sq <= p_norm_bound_sq.floor() as u64 , "||p||_2^2 = {} is not <= 128*beta^2 = {}", norm_p_sq, p_norm_bound_sq);

    let psi = merlin.challenge_vectors::<R::BaseRing, R::BaseRing>(num_ct_constraints, crs.num_aggregs).expect("error extracting verifier message 2 (psi) from transcript");
    let omega = merlin.challenge_vectors::<R::BaseRing, R::BaseRing>(256, crs.num_aggregs).expect("error extracting verifier message 2 (omega) from transcript");

    let b__ = merlin.next_vec::<R>(crs.num_aggregs).expect("error extracting prover message 3 from transcript");

    for k in 0..crs.num_aggregs {
        let mut rhs_k = omega[k].dot(&p);
        for l in 0..num_ct_constraints {
            rhs_k += psi[k][l] * instance.ct_quad_dot_prod_funcs[l].b;
        }
        check_eq!(b__[k].coeffs()[0], rhs_k, "constant coeff of b''^(k) = {:?} is not equal to rhs = {:?}", b__[k].coeffs()[0], rhs_k);
    }

    let alpha = merlin.challenge_vector::<R, R>(num_constraints).expect("error extracting verifier message 3 (alpha) from transcript");
    let beta = merlin.challenge_vector::<R, R>(crs.num_aggregs).expect("error extracting verifier message 3 (beta) from transcript");

    let u_2 = merlin.next_vector::<R>(crs.k2).expect("error extracting prover message 4 from transcript");

    let c = merlin.challenge_vec::<R, LabradorChallengeSet<R>>(crs.r).expect("error extracting verifier message 4 from transcript");

    // Compute phi
    let mut phi__ = vec![vec![Vector::<R>::zeros(crs.n); crs.r]; crs.num_aggregs];
    let K = instance.quad_dot_prod_funcs.len();
    for k in 0..crs.num_aggregs {
        for i in 0..crs.r {
            for j in 0..crs.r {
                // Compute a''_{ij}^{(k)}
                let mut a_ijk = R::zero();
                for l in 0..K {
                    a_ijk += instance.ct_quad_dot_prod_funcs[l].A[(i, j)] * psi[k][l];
                }
            }
            // Compute vec{phi}''_i^{(k)}
            for l in 0..K {
                phi__[k][i] += mul_basescalar_vector(psi[k][l], &instance.ct_quad_dot_prod_funcs[l].phi[i]);
            }
            for j in 0..256 {
                phi__[k][i] += mul_basescalar_vector(omega[k][j], &R::sigma_vec(&Pi[i].row(j).transpose()));
            }
        }
    }
    let mut phi = vec![Vector::<R>::zeros(n); r];
    for i in 0..r {
        for k in 0..K {
            phi[i] += &instance.quad_dot_prod_funcs[k].phi[i] * alpha[k];
        }
        for k in 0..crs.num_aggregs {
            phi[i] += &phi__[k][i] * beta[k];
        }
    }

    Ok(BaseTranscript::new_core(crs, instance, u_1, Pi, psi, omega, b__, alpha, beta, u_2, c, phi))
}


pub fn verify_principal_relation<R: PolyRing>(merlin: &mut LatticeMerlin, mut instance: &PrincipalRelation<R>, mut crs: &CommonReferenceString<R>) -> Result<(), InvalidTag>
    where LabradorChallengeSet<R>: FromRandomBytes<R>, WeightedTernaryChallengeSet<R>: FromRandomBytes<R>
{
    // Init Fiat-Shamir transcript
    // TODO: add public statement
    // merlin.ratchet()?;

    let mut transcript: BaseTranscript<R>;
    let mut instance = instance.to_owned();
    let mut crs = crs.to_owned();
    loop {
        transcript = verify_core(&crs, &instance, merlin)?;
        let recurse = recurse(&transcript);
        if recurse {
            (instance, crs, _) = fold_instance(&transcript, false);
        } else {
            break;
        }
    }

    // Final checks
    let z = merlin.next_vector::<R>(crs.n).expect("error extracting prover message 5 (z) from transcript");
    let t = merlin.next_vectors::<R>(crs.k, crs.t1).expect("error extracting prover message 5 (t) from transcript");
    let G = merlin.next_symmetric_matrix::<R>(crs.r).expect("error extracting prover message 5 (G) from transcript");
    let H = merlin.next_symmetric_matrix::<R>(crs.r).expect("error extracting prover message 5 (H) from transcript");
    transcript.z.replace(z);
    transcript.t.replace(t);
    transcript.G.replace(G);
    transcript.H.replace(H);

    verify_final(&transcript)?;

    Ok(())
}