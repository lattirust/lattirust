#![allow(non_snake_case)]

use std::cmp::max;

use crate::labrador::prover::Witness;
use crate::labrador::setup::CommonReferenceString;
use crate::labrador::util::{chunk_pad, concat, flatten_symmetric_matrix, flatten_vec_vector, shift_right};
use crate::lattice_arithmetic::balanced_decomposition::decompose_balanced_vec;
use crate::lattice_arithmetic::challenge_set::labrador_challenge_set::LabradorChallengeSet;
use crate::lattice_arithmetic::matrix::{Matrix, Vector};
use crate::lattice_arithmetic::poly_ring::PolyRing;
use crate::relations::labrador::principal_relation::{PrincipalRelation, QuadDotProdFunction};

/// A (subset of a) transcript of one execution of the core Labrador protocol, as needed to compute the next instance/crs/witness for recursion
pub struct BaseTranscript<'a, R: PolyRing> {
    pub instance: &'a PrincipalRelation<R>,
    pub crs: &'a CommonReferenceString<R>,
    pub(crate) u_1: Option<Vector<R>>,
    pub(crate) Pi: Option<Vec<Matrix<R>>>,
    pub(crate) p: Option<Vector<R::BaseRing>>,
    pub(crate) psi: Option<Vec<Vector<R::BaseRing>>>,
    pub(crate) omega: Option<Vec<Vector<R::BaseRing>>>,
    pub(crate) b__: Option<Vec<R>>,
    pub(crate) alpha: Option<Vector<R>>,
    pub(crate) beta: Option<Vector<R>>,
    pub(crate) u_2: Option<Vector<R>>,
    pub(crate) c: Option<Vec<R>>,
    pub(crate) z: Option<Vector<R>>,
    pub(crate) t: Option<Vec<Vector<R>>>,
    pub(crate) G: Option<Vec<Vec<R>>>,
    pub(crate) H: Option<Vec<Vec<R>>>,

    // Note: phi is not part of the transcript, but it is a linear combination of values in the transcript
    pub(crate) phi: Option<Vec<Vector<R>>>,

}

impl<'a, R: PolyRing> BaseTranscript<'a, R> {
    pub fn init(crs: &'a CommonReferenceString<R>, instance: &'a PrincipalRelation<R>) -> BaseTranscript<'a, R> {
        BaseTranscript {
            instance: &instance,
            crs: &crs,
            u_1: None,
            Pi: None,
            p: None,
            psi: None,
            omega: None,
            b__: None,
            alpha: None,
            beta: None,
            u_2: None,
            c: None,
            z: None,
            t: None,
            G: None,
            H: None,
            phi: None,
        }
    }

    pub fn new_core(crs: &'a CommonReferenceString<R>, instance: &'a PrincipalRelation<R>,
                    u_1: Vector<R>, Pi: Vec<Matrix<R>>, psi: Vec<Vector<R::BaseRing>>, omega: Vec<Vector<R::BaseRing>>, b__: Vec<R>, alpha: Vector<R>, beta: Vector<R>, u_2: Vector<R>, c: Vec<R>, phi: Vec<Vector<R>>) -> BaseTranscript<'a, R> {
        BaseTranscript {
            instance: &instance,
            crs: &crs,
            u_1: Some(u_1),
            Pi: Some(Pi),
            p: None,
            psi: Some(psi),
            omega: Some(omega),
            b__: Some(b__),
            alpha: Some(alpha),
            beta: Some(beta),
            u_2: Some(u_2),
            c: Some(c),
            z: None,
            t: None,
            G: None,
            H: None,
            phi: Some(phi),
        }
    }

    pub fn new(crs: &'a CommonReferenceString<R>, instance: &'a PrincipalRelation<R>,
               u_1: Vector<R>, Pi: Vec<Matrix<R>>, p: Vector<R::BaseRing>, psi: Vec<Vector<R::BaseRing>>, omega: Vec<Vector<R::BaseRing>>, b__: Vec<R>, alpha: Vector<R>, beta: Vector<R>, u_2: Vector<R>, c: Vec<R>, z: Vector<R>, t: Vec<Vector<R>>, G: Vec<Vec<R>>, H: Vec<Vec<R>>, phi: Vec<Vector<R>>) -> BaseTranscript<'a, R> {
        BaseTranscript {
            instance: &instance,
            crs: &crs,
            u_1: Some(u_1),
            Pi: Some(Pi),
            p: Some(p),
            psi: Some(psi),
            omega: Some(omega),
            b__: Some(b__),
            alpha: Some(alpha),
            beta: Some(beta),
            u_2: Some(u_2),
            c: Some(c),
            z: Some(z),
            t: Some(t),
            G: Some(G),
            H: Some(H),
            phi: Some(phi),
        }
    }
}

pub fn recurse<R: PolyRing>(transcript: &BaseTranscript<R>) -> bool {
    let crs = transcript.crs;
    let m = crs.r * crs.t1 * crs.k + (crs.t1 + crs.t2) * (crs.r * (crs.r + 1)).div_ceil(2);
    let mu = 2;
    let nu = (2. * crs.n as f64 / m as f64).round() as usize;

    // Recurse if 2*n >> m // TODO: select a recursion threshold with a bit more theoretical justification
    return 2 * 2 * crs.n > m;
}

pub fn next_norm_bound_sq<R: PolyRing>(transcript: &BaseTranscript<R>) -> f64 {
    let crs = transcript.crs;
    let b_f = Into::<i64>::into(crs.decomposition_basis) as f64;
    let b_sq = b_f * b_f;
    let challenge_variance = LabradorChallengeSet::<R>::challenge_poly_sum_coeffs_variance();
    let gamma_sq = crs.norm_bound_squared * challenge_variance;
    let gamma_1_sq = (b_sq * crs.t1 as f64) / 12. * (crs.r * crs.k * crs.d) as f64 + (b_sq * crs.t2 as f64) / 12. * ((crs.r * (crs.r + 1)).div_ceil(2) * crs.d) as f64;
    let gamma_2_sq = (b_sq * crs.t1 as f64) / 12. * ((crs.r * (crs.r + 1)).div_ceil(2) * crs.d) as f64;
    let beta_next_sq: f64 = (2. / b_sq) * gamma_sq + gamma_1_sq * gamma_2_sq;
    beta_next_sq
}

pub fn next_norm_bound<R: PolyRing>(transcript: &BaseTranscript<R>) -> f64 {
    next_norm_bound_sq(&transcript).sqrt()
}

pub fn fold_instance<'a, R: PolyRing>(transcript: &BaseTranscript<R>, compute_witness: bool) -> (PrincipalRelation<R>, CommonReferenceString<R>, Option<Witness<R>>) {
    assert!(recurse(transcript));
    let crs = transcript.crs;

    // Generate instance for next iteration of the protocol
    // Set nu, mu such that 2 * n_next ≈ m_next,
    // i.e., set mu = 2 and nu ≈ 2 * n / m
    let m = crs.r * crs.t1 * crs.k + (crs.t1 + crs.t2) * (crs.r * (crs.r + 1)).div_ceil(2);
    let mu = 2;
    let nu = (2. * crs.n as f64 / m as f64).round() as usize;
    let n_next = max(crs.n.div_ceil(nu), m.div_ceil(mu));
    let m_next = (m as f64 / 2.).round() as usize;
    let r_next = 2 * nu + m;


    let mut quad_dot_prod_funcs_next = Vec::<QuadDotProdFunction<R>>::with_capacity(crs.k + crs.k1 + crs.k2 + 3);

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

    let c = transcript.c.as_ref().expect("c not available");

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
        phis_next[0..nu].clone_from_slice(crate::labrador::util::split(&A_l, nu).as_slice());
        // <A_l * b, z_1>
        phis_next[nu..2 * nu].clone_from_slice(crate::labrador::util::split(&(&A_l * b_ring), nu).as_slice());

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

        quad_dot_prod_funcs_next.push(QuadDotProdFunction::<R>::new(Matrix::<R>::zeros(r_next, r_next), phis_next, R::zero()));
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

    quad_dot_prod_funcs_next.push(QuadDotProdFunction::<R>::new(A, phis_next, R::zero()));

    // Constraints for sum_{i in [r]} <phi_i, z> c_i = sum_{i, j in [r]} h_ij c_i c_j
    let phi = transcript.phi.as_ref().expect("phi not available");
    debug_assert_eq!(phi.len(), crs.n);
    let mut phis_next = vec![Vector::<R>::zeros(n_next); r_next];

    let mut phi_lc_0 = Vector::<R>::zeros(crs.n);
    let mut phi_lc_1 = Vector::<R>::zeros(crs.n);

    for i in 0..crs.r {
        phi_lc_0 += &phi[i] * c[i];
        phi_lc_1 += &phi[i] * c[i] * b_ring;
    }
    // <sum_{i in [r]} phi_i * c_i, z_0>
    phis_next[0..nu].clone_from_slice(crate::labrador::util::split(&phi_lc_0, nu).as_slice());
    // <sum_{i in [r]} phi_i * c_i * b, z_1>
    phis_next[nu..2 * nu].clone_from_slice(crate::labrador::util::split(&phi_lc_1, nu).as_slice());

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

    quad_dot_prod_funcs_next.push(QuadDotProdFunction::<R>::new(Matrix::<R>::zeros(r_next, r_next), phis_next, R::zero()));

    // TODO: add constraints for sum_{i,j in [r]} a_ij * g_ij + sum_{i in [r]} h_ii = b

    // Constraints for u_1
    let u_1 = transcript.u_1.as_ref().expect("u_1 not available");
    for l in 0..crs.k1 {
        let mut phis = vec![Vector::<R>::zeros(n_next); r_next];
        let mut B_flat = Vec::<R>::with_capacity(crs.r * crs.t1 * crs.n);
        for i in 0..crs.r {
            for k in 0..crs.t1 {
                B_flat.extend(crs.B[i][k].row(l).transpose().data.as_vec());
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


        quad_dot_prod_funcs_next.push(QuadDotProdFunction::<R>::new(Matrix::<R>::zeros(r_next, r_next), phis, u_1[l]));
    }

    // Constraints for u_2
    let u_2 = transcript.u_2.as_ref().expect("u_2 not available");
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

        quad_dot_prod_funcs_next.push(QuadDotProdFunction::<R>::new(Matrix::<R>::zeros(r_next, r_next), phis, u_2[l]));
    }

    let next_norm_bound_squared = next_norm_bound_sq(&transcript);

    let num_quad_constraints = quad_dot_prod_funcs_next.len();
    let instance_next = PrincipalRelation::<R> {
        quad_dot_prod_funcs: quad_dot_prod_funcs_next,
        ct_quad_dot_prod_funcs: vec![],
    };

    // TODO: check
    let crs_next = CommonReferenceString::<R>::new(r_next, n_next, crs.d, next_norm_bound_squared, crs.k, crs.k1, crs.k2, num_quad_constraints, 0, crs.decomposition_basis);

    let compute_witness = false;

    let witness = if compute_witness {
        let (z, t, G, H) = (transcript.z.as_ref().unwrap(), transcript.t.as_ref().unwrap(), transcript.G.as_ref().unwrap(), transcript.H.as_ref().unwrap());
        let z_decomp = decompose_balanced_vec(&z, crs.decomposition_basis, Some(2usize));
        assert_eq!(z_decomp.len(), 2);

        let v = concat(&[&flatten_vec_vector(&t), &flatten_symmetric_matrix(&G), &flatten_symmetric_matrix(&H)]);
        let z_0_split = crate::labrador::util::split(&z_decomp[0], nu);
        let z_1_split = crate::labrador::util::split(&z_decomp[1], nu);
        let v_split = crate::labrador::util::split(&v, mu);
        Some(Witness::<R> {
            s: vec![z_0_split, z_1_split, v_split].concat(),
        })
    } else {
        None
    };

    (instance_next, crs_next, witness)
}