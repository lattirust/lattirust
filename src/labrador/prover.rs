#![allow(non_snake_case)]

use std::cmp::max;

use rayon::prelude::*;

use crate::labrador::setup::CommonReferenceString;
use crate::labrador::util::*;
pub use crate::labrador::witness::Witness;
use crate::lattice_arithmetic::balanced_decomposition::{decompose_balanced_polyring, decompose_balanced_vec};
use crate::lattice_arithmetic::challenge_set::labrador_challenge_set::LabradorChallengeSet;
use crate::lattice_arithmetic::challenge_set::weighted_ternary::WeightedTernaryChallengeSet;
use crate::lattice_arithmetic::matrix::{Matrix, norm_sq_vec, Vector};
use crate::lattice_arithmetic::poly_ring::PolyRing;
use crate::lattice_arithmetic::traits::{FromRandomBytes, WithLog2};
use crate::nimue::arthur::LatticeArthur;
use crate::relations::labrador::principal_relation::{PrincipalRelation, QuadDotProdFunction};

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
    phi__: Option<Vec<Vec<Vector<R>>>>,
    alpha: Option<Vector<R>>,
    beta: Option<Vector<R>>,
    phi: Option<Vec<Vector<R>>>,
    c: Option<Vec<R>>,
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
            phi__: None,
            alpha: None,
            beta: None,
            phi: None,
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
        for k in 0..instance.quad_dot_prod_funcs.len() {
            phi_i += &instance.ct_quad_dot_prod_funcs[k].phi[i] * alpha[k];
        }
        for k in 0..crs.num_aggregs {
            phi_i += &phi__[k][i] * beta[k];
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

    state.phi.replace(phi);
    state.H.replace(H);
    u_2
}


fn prove_5<R: PolyRing>(state: &ProverState<R>) -> (Vector<R>, Vec<Vector<R>>, Vec<Vec<R>>, Vec<Vec<R>>) {
    let (_, witness, crs) = (&state.instance, &state.witness, &state.crs);
    let c = state.c.as_ref().expect("c not available");
    let mut z = Vector::<R>::zeros(crs.n);
    debug_assert_eq!(c.len(), crs.r);
    for i in 0..crs.r {
        z += &witness.s[i] * c[i];
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
    arthur.absorb_lower_triangular_matrix(&G).expect("error absorbing prover message 5 (G)");
    arthur.absorb_lower_triangular_matrix(&H).expect("error absorbing prover message 5 (H)");

    // Generate instance for next iteration of the protocol
    // Set nu, mu such that 2 * n_next ≈ m_next,
    // i.e., set mu = 2 and nu ≈ 2 * n / m
    let m = crs.r * crs.t1 * crs.k + (crs.t1 + crs.t2) * (crs.r * (crs.r + 1)).div_ceil(2);
    let mu = 2;
    let nu = (2. * crs.n as f64 / m as f64).round() as usize;
    let n_next = max(crs.n.div_ceil(nu), m.div_ceil(mu));
    let m_next = (m as f64 / 2.).round() as usize;
    let r_next = 2 * nu + m;


    let z_decomp = decompose_balanced_vec(&z, crs.decomposition_basis, Some(2usize));
    assert_eq!(z_decomp.len(), 2);

    // TODO: concat decompositions instead
    let v = concat(&[flatten_vec_vector(&t), flatten_symmetric_matrix(&G), flatten_symmetric_matrix(&H)]);
    let z_0_split = split(&z_decomp[0], nu);
    let z_1_split = split(&z_decomp[1], nu);
    let v_split = split(&v, mu);

    let mut quad_dot_prod_funcs_next = Vec::<QuadDotProdFunction<R>>::with_capacity(crs.k + crs.k1 + crs.k2 + 3);
    // TODO

    //                           r * t1 * n
    // ┌─────────────┬─────────────┬─────────────┬──────────────┐
    // │    t_0 (0)  │    t_0 (1)  │     ...     │ t_{r-1} (t-1)│
    // └─────────────┴─────────────┴─────────────┴──────────────┘
    // ┌────────┬────────┬────────┐                          ┌───────────┐
    // │  s_2nu │s_2nu+1 │s_2nu+2 │       ...                │s_2nu+tau-1│
    // └────────┴────────┴────────┘                          └───────────┘
    //                        tau * n_next
    // tau = ceil(r * t1 * n / n_next)
    let tau = (crs.r * crs.t1 * crs.n).div_ceil(n_next);

    let offset_t_g = tau * n_next - crs.r * crs.t1 * crs.n;
    let gamma = ((crs.r * (crs.r + 1)).div_ceil(2) * crs.t2 + offset_t_g).div_ceil(n_next);

    let offset_g_h = gamma * n_next - ((crs.r * (crs.r + 1)).div_ceil(2) * crs.t2 + offset_t_g);
    let delta = ((crs.r * (crs.r + 1)).div_ceil(2) * crs.t1 + offset_g_h).div_ceil(n_next);

    let c = state.c.expect("c not available");

    let b_ring = R::from(crs.decomposition_basis);
    let mut b_pows = Vec::<R>::with_capacity(max(crs.t1, crs.t2));
    b_pows[0] = R::one();
    for i in 1..max(crs.t1, crs.t2) {
        b_pows.push(b_pows[i - 1] * b_pows[i - 1]);
    }

    // Constraints for Az = sum_{i in [r]} c_i t_i
    for l in 0..crs.r {
        let mut phis_next = vec![Vector::<R>::zeros(n_next); r_next];
        let A_l = crs.A.row(l).transpose();
        // <A_l, z_0>
        phis_next[0..nu].clone_from_slice(split(&A_l, nu).as_slice());
        // <A_l * b, z_1>
        phis_next[nu..2 * nu].clone_from_slice(split(&(&A_l * b_ring), nu).as_slice());

        // <(t_j^(k)_l)_{j, k}, (c_j * b^k)_{j, k})>
        let mut c_vec = vec![R::zero(); crs.r * crs.t1 * crs.n];
        for i in 0..crs.r {
            for k in 0..crs.t1 {
                c_vec[i * crs.r * crs.t1 * crs.n + k * crs.t1 * crs.n + l] = c[i] * b_pows[k];
            }
        }
        let c_vec_split = chunk_pad(&Vector::<R>::from_vec(c_vec), n_next);
        assert_eq!(c_vec_split.len(), tau);
        phis_next[2 * nu..2 * nu + tau].clone_from_slice(c_vec_split.as_slice());

        quad_dot_prod_funcs_next.push(QuadDotProdFunction::<R> {
            A: Matrix::<R>::zeros(r_next, r_next),
            phi: phis_next,
            b: R::zero(),
        });
    }

    // Constraints for <z, z> = sum_{i, j in [r]} g_ij c_i c_j
    let mut A = Matrix::<R>::zeros(r_next, r_next);
    let b_sq = R::from(crs.decomposition_basis * crs.decomposition_basis);
    for i in 0..r_next {
        A[(i, i)] = R::one();
        A[(i, i + nu)] = b_ring;
        A[(i + nu, i)] = A[(i, i + nu)];
        A[(i + nu, i + nu)] = b_sq;
    }
    let mut phis_next = vec![Vector::<R>::zeros(n_next); r_next];
    let two = R::from(2u128);

    let mut c_prods = Vec::<R>::with_capacity((crs.r * (crs.r + 1)).div_ceil(2) * crs.t2);
    for i in 0..crs.r {
        for j in 0..i + 1 {
            for k in 0..crs.t2 {
                if i == j {
                    c_prods.push(c[i] * c[i] * b_pows[k]);
                } else {
                    c_prods.push(c[i] * c[j] * b_pows[k] * two);
                }
            }
        }
    }
    c_prods = shift_right(&c_prods, offset_t_g);
    let c_prods_flat_split = chunk_pad(&Vector::<R>::from_vec(c_prods), n_next);
    assert_eq!(c_prods_flat_split.len(), gamma);
    phis_next[2 * nu + tau - 1..2 * nu + tau - 1 + gamma].clone_from_slice(c_prods_flat_split.as_slice());

    quad_dot_prod_funcs_next.push(QuadDotProdFunction::<R> {
        A,
        phi: phis_next,
        b: R::zero(),
    });

    // Constraints for sum_{i in [r]} <phi_i, z> c_i = sum_{i, j in [r]} h_ij c_i c_j
    let phi = state.phi.expect("phi not available");
    debug_assert_eq!(phi.len(), crs.n);
    let mut phis_next = vec![Vector::<R>::zeros(n_next); r_next];

    let mut phi_lc_0 = Vector::<R>::zeros(crs.n);
    let mut phi_lc_1 = Vector::<R>::zeros(crs.n);

    for i in 0..crs.r {
        phi_lc_0 += &phi[i] * c[i];
        phi_lc_1 += &phi[i] * c[i] * b_ring;
    }
    // <sum_{i in [r]} phi_i * c_i, z_0>
    phis_next[0..nu].clone_from_slice(split(&phi_lc_0, nu).as_slice());
    // <sum_{i in [r]} phi_i * c_i * b, z_1>
    phis_next[nu..2 * nu].clone_from_slice(split(&phi_lc_1, nu).as_slice());

    // <(c_i * c_j * b^k)_{i, j, k}, h>
    let mut c_prods = Vec::<R>::with_capacity((crs.r * (crs.r + 1)).div_ceil(2) * crs.t2);
    for i in 0..crs.r {
        for j in 0..i + 1 {
            for k in 0..crs.t1 {
                if i == j {
                    c_prods.push(c[i] * c[i] * b_pows[k]);
                } else {
                    c_prods.push(c[i] * c[j] * b_pows[k] * two);
                }
            }
        }
    }
    c_prods = shift_right(&c_prods, offset_g_h);
    let c_prods_flat_split = chunk_pad(&Vector::<R>::from_vec(c_prods), n_next);
    assert_eq!(c_prods_flat_split.len(), delta);
    phis_next[2 * nu + tau + gamma - 1..2 * nu + tau + gamma - 1 + delta].clone_from_slice(c_prods_flat_split.as_slice());

    quad_dot_prod_funcs_next.push(QuadDotProdFunction::<R> {
        A: Matrix::<R>::zeros(r_next, r_next),
        phi: phis_next,
        b: R::zero(),
    });

    // Constraints for u_1
    for l in 0..crs.k1 {
        let mut phis = vec![Vector::<R>::zeros(n_next); r_next];
        let mut B_flat = Vec::<R>::with_capacity(crs.r * crs.t1 * crs.n);
        for i in 0..crs.r {
            for k in 0..crs.t1 {
                B_flat.extend(crs.B[i][k].row(l).transpose().data.as_vec());

// // split B_ik's l-th row into vectors os size n_next, and set them to phi[2*nu..2*nu+tau]
// let B_ikl_split = chunk_pad(&crs.B[i][k].row(l).transpose(), n_next);
// let tau = (crs.r * crs.t1 * crs.n).div_ceil(n_next);
// assert_eq!(B_ikl_split.len(), tau);
// phis[2*nu..2*nu+tau].clone_from_slice(B_ikl_split.as_slice());
            }
        }
        let B_flat_split = chunk_pad(&Vector::<R>::from_vec(B_flat), n_next);
        phis[2 * nu..2 * nu + tau].clone_from_slice(B_flat_split.as_slice());

        let mut C_flat = Vec::<R>::with_capacity((crs.r * (crs.r + 1)).div_ceil(2) * crs.t2);
        for i in 0..crs.r {
            for j in 0..i + 1 {
                for k in 0..crs.t2 {
                    C_flat.push(crs.C[i][j][k][l]);
                }
            }
        }
        C_flat = shift_right(&C_flat, offset_t_g);
        let C_flat_split = chunk_pad(&Vector::<R>::from_vec(C_flat), n_next);
        assert_eq!(C_flat_split.len(), gamma);
        phis[2 * nu + tau - 1..2 * nu + tau - 1 + gamma].clone_from_slice(C_flat_split.as_slice());


        quad_dot_prod_funcs_next.push(QuadDotProdFunction::<R> {
            A: Matrix::<R>::zeros(r_next, r_next),
            phi: phis,
            b: u_1[l],
        });
    }

// Constraints for u_2
    for l in 0..crs.k2 {
        let mut phis = vec![Vector::<R>::zeros(n_next); r_next];

        let mut D_flat = Vec::<R>::with_capacity((crs.r * (crs.r + 1)).div_ceil(2) * crs.t1);
        for i in 0..crs.r {
            for j in 0..i + 1 {
                for k in 0..crs.t1 {
                    D_flat.push(crs.D[i][j][k][l]);
                }
            }
        }
        let D_flat = shift_right(&D_flat, offset_g_h);
        let D_flat_split = chunk_pad(&Vector::<R>::from_vec(D_flat), n_next);
        assert_eq!(D_flat_split.len(), delta);
        phis[2 * nu + tau + gamma - 1..2 * nu + tau + gamma - 1 + delta].clone_from_slice(D_flat_split.as_slice());

        quad_dot_prod_funcs_next.push(QuadDotProdFunction::<R> {
            A: Matrix::<R>::zeros(r_next, r_next),
            phi: phis,
            b: u_2[l],
        });
    }

    let b_f = i128::from(crs.decomposition_basis) as f64;
    let b_sq = b_f * b_f;
    let challenge_variance = LabradorChallengeSet::<R>::challenge_poly_sum_coeffs_variance();
    let gamma_sq = crs.beta * challenge_variance.sqrt();
    let gamma_1_sq = (b_sq * crs.t1 as f64) / 12. * (crs.r * crs.k * crs.d) as f64 + (b_sq * crs.t2 as f64) / 12. * ((crs.r * (crs.r + 1)).div_ceil(2) * crs.d) as f64;
    let gamma_2_sq = (b_sq * crs.t1 as f64) / 12. * ((crs.r * (crs.r + 1)).div_ceil(2) * crs.d) as f64;
    let beta_next: f64 = ((2. / b_sq) * gamma_sq + gamma_1_sq * gamma_2_sq).sqrt();

    let instance_next = PrincipalRelation::<R> {
        r: r_next,
        n: n_next,
        norm_bound: beta_next,
        quad_dot_prod_funcs: quad_dot_prod_funcs_next,
        ct_quad_dot_prod_funcs: vec![],
    };

    let witness_next = Witness::<R> {
        s: vec![z_0_split, z_1_split, v_split].concat(),
    };

    prove_principal_relation(arthur, instance_next, witness_next, crs)?;

    Ok(arthur.transcript())
}