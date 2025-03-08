#![allow(non_snake_case)]

use std::cmp::max;
use std::fmt::Debug;

use serde::Serialize;

use lattirust_arithmetic::balanced_decomposition::{
    decompose_balanced_vec_polyring, DecompositionFriendlySignedRepresentative,
};
use lattirust_arithmetic::linear_algebra::Matrix;
use lattirust_arithmetic::linear_algebra::Vector;
use lattirust_arithmetic::ring::PolyRing;
use lattirust_arithmetic::ring::representatives::WithSignedRepresentative;
use relations::principal_relation::{Index, Instance, QuadDotProdFunction, Witness};

use crate::shared::BaseTranscript;
use crate::util::{
    chunk_pad, concat, flatten_symmetric_matrix, flatten_vec_vector, mul_matrix_basescalar,
    shift_right,
};

/// Common reference string for one round of the LaBRADOR protocol
#[derive(Clone, Debug, Serialize)]
pub struct CommonReferenceString<R: PolyRing> {
    /// Number of witness vectors
    pub r: usize,
    /// Number of entries in a witness vector
    pub n: usize,
    /// Dimension of the polynomial ring Z_q[X]/(X^d + 1)
    pub d: usize,
    /// Square of the L2-norm bound on the concatenation of witness vectors
    pub norm_bound_squared: f64,
    /// Size of the first-level commitment
    pub k: usize,
    /// Size of the second-level commitment for elements decomposed in basis `b1`
    pub k1: usize,
    /// Size of the second-level commitment for elements decomposed in basis `b2`
    pub k2: usize,
    /// Length of decompositions in basis `b1`
    pub t1: usize,
    /// Length of decompositions in basis `b2`
    pub t2: usize,
    /// Number of aggregated constraints when reducing constant constraints, equal to  `security parameter / log(q)`
    pub num_aggregs: usize,
    /// Number of quadratic-linear constraints
    pub num_constraints: usize,
    /// Number of quadratic-linear constraints on constant coefficients
    pub num_constant_constraints: usize,
    /// First-level commitment matrix, of size k x n
    pub A: Matrix<R>,
    /// Second-level commitment matrices for first decomposition basis, (r x t1) instances, each of size k1 x k
    pub B: Vec<Vec<Matrix<R>>>,
    /// Second-level commitment matrices for second decomposition basis, (r x r x t2) instances, each of size k2 x 1
    pub C: Vec<Vec<Vec<Vector<R>>>>,
    /// Second-level commitment matrices for first decomposition basis, (r x r x t1) instances, each of size k2 x 1
    pub D: Vec<Vec<Vec<Vector<R>>>>,
    /// Decomposition basis for z-vectors, roughly equal to `b1` and `b2`
    pub b: u128,
    /// Decomposition basis for first-level commitments (t-vectors), roughly equal to `b` and `b2`
    pub b1: u128,
    /// Decomposition basis for inner product terms (g), roughly equal to `b1` and `b2`
    pub b2: u128,
    /// A reference to the CRS for the next recursive round, or `None` if this is the CRS for the last round
    pub next_crs: Option<Box<CommonReferenceString<R>>>,
}

pub fn fold_instance<R: PolyRing>(
    transcript: &BaseTranscript<R>,
    compute_witness: bool,
) -> (Index<R>, Instance<R>, Option<Witness<R>>)
where
    <R as PolyRing>::BaseRing: WithSignedRepresentative,
    <<R as PolyRing>::BaseRing as WithSignedRepresentative>::SignedRepresentative:
        DecompositionFriendlySignedRepresentative,

{
    // assert!(transcript.crs.recurse());
    let crs = transcript.crs;

    // Generate instance for next iteration of the protocol
    // Set nu, mu such that 2 * n_next ≈ m_next,
    // i.e., set mu = 2 and nu ≈ 2 * n / m
    let m = crs.r * crs.t1 * crs.k + (crs.t1 + crs.t2) * (crs.r * (crs.r + 1)).div_ceil(2);
    let mu = 2;
    let nu = (2. * crs.n as f64 / m as f64).round() as usize;
    let n_next = max(crs.n.div_ceil(nu), m.div_ceil(mu));
    let _m_next = (m as f64 / 2.).round() as usize;
    let r_next = 2 * nu + m;

    let mut quad_dot_prod_funcs_next =
        Vec::<QuadDotProdFunction<R>>::with_capacity(crs.k + crs.k1 + crs.k2 + 3);

    // The new witness is the concatenation of the flattened vectors z^(0), z^(1), (t_i^(j))_{i in [r], j in [t1]}, (g_ij)_{i >= j in [r]}, (h_ij)_{i >= j in [r]}

    // │<---------------------- r * t1 * n -------------------->│
    // ┌─────────────┬─────────────┬─────────────┬──────────────┐        ┌─────────────┬─────────────┬─────────────┬──────────────┐
    // │    t_0^(0)  │    t_0^(1)  │     ...     │t_{r-1}^(t1-1)│<-o_tg->│     g_0,0   │     g_0,1   │     ...     │     g_r,r    │
    // └─────────────┴─────────────┴─────────────┴──────────────┘        └─────────────┴─────────────┴─────────────┴──────────────┘
    // ┌────────┬────────┬────────┬                          ┬───────────┬─────────┬
    // │  s_2nu │s_2nu+1 │s_2nu+2 │             ...          │s_2nu+tau-1│s_2nu+tau│   ...
    // └────────┴────────┴────────┘                          └───────────┴─────────┴
    // │<------------------------ tau * n_next ------------------------->│
    // tau = ceil(r * t1 * n / n_next)
    let tau = (crs.r * crs.t1 * crs.n).div_ceil(n_next); // Number of n_next-sized s-vectors over R_q needed to represent all t_i^(j)'s

    let offset_t_g = tau * n_next - crs.r * crs.t1 * crs.n; // Offset to align the t_i^(j)'s with the g_ij's (o_tg in the diagram above)
    let gamma = ((crs.r * (crs.r + 1)).div_ceil(2) * crs.t2 + offset_t_g).div_ceil(n_next); // Number of n_next-sized s-vectors over R_q needed to represent all g_ij's

    let offset_g_h = gamma * n_next - ((crs.r * (crs.r + 1)).div_ceil(2) * crs.t2 + offset_t_g); // Offset to align the g_ij's with the h_ij's
    let eta = ((crs.r * (crs.r + 1)).div_ceil(2) * crs.t1 + offset_g_h).div_ceil(n_next); // Number of n_next-sized s-vectors over R_q needed to represent all h_ij's

    let c = transcript.c.as_ref().expect("c not available");

    let b_ring = R::try_from(crs.b).unwrap();
    let mut b_pows = Vec::<R>::with_capacity(max(crs.t1, crs.t2));
    b_pows[0] = R::one();
    for i in 1..max(crs.t1, crs.t2) {
        b_pows.push(b_pows[i - 1] * b_pows[i - 1]);
    }
    let two = R::try_from(2u128).unwrap();

    // Constraints for Az = sum_{i in [r]} c_i t_i
    {
        for l in 0..crs.k {
            let mut phis_next = vec![Vector::<R>::zeros(n_next); r_next];
            let A_l = crs.A.row(l).transpose(); // <A_l, z_0>
            phis_next[0..nu].clone_from_slice(crate::util::split(&A_l, nu).as_slice()); // <A_l * b, z_1>
            phis_next[nu..2 * nu]
                .clone_from_slice(crate::util::split(&(&A_l * b_ring), nu).as_slice());

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

            quad_dot_prod_funcs_next.push(QuadDotProdFunction::<R>::new(
                Matrix::<R>::zeros(r_next, r_next),
                phis_next,
                R::zero(),
            ));
        }
    }

    // Constraint for <z, z> = sum_{i, j in [r]} g_ij c_i c_j
    {
        let mut A = Matrix::<R>::zeros(r_next, r_next);
        let b_sq = R::try_from(crs.b * crs.b).unwrap();
        for i in 0..r_next {
            A[(i, i)] = R::one();
            A[(i, i + nu)] = b_ring;
            A[(i + nu, i)] = A[(i, i + nu)];
            A[(i + nu, i + nu)] = b_sq;
        }
        let mut phis_next = vec![Vector::<R>::zeros(n_next); r_next];

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
        phis_next[2 * nu + tau - 1..2 * nu + tau - 1 + gamma]
            .clone_from_slice(c_prods_flat_split.as_slice());

        quad_dot_prod_funcs_next.push(QuadDotProdFunction::<R>::new(A, phis_next, R::zero()));
    }

    // Constraint for sum_{i in [r]} <phi_i, z> c_i = sum_{i, j in [r]} h_ij c_i c_j
    {
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
        phis_next[0..nu].clone_from_slice(crate::util::split(&phi_lc_0, nu).as_slice());
        // <sum_{i in [r]} phi_i * c_i * b, z_1>
        phis_next[nu..2 * nu].clone_from_slice(crate::util::split(&phi_lc_1, nu).as_slice());

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
        assert_eq!(c_prods_flat_split.len(), eta);
        phis_next[2 * nu + tau + gamma - 1..2 * nu + tau + gamma - 1 + eta]
            .clone_from_slice(c_prods_flat_split.as_slice());

        quad_dot_prod_funcs_next.push(QuadDotProdFunction::<R>::new(
            Matrix::<R>::zeros(r_next, r_next),
            phis_next,
            R::zero(),
        ));
    }

    // Constraint for sum_{i,j in [r]} a_ij * g_ij + sum_{i in [r]} h_ii = b
    {
        let b__ = transcript.b__.as_ref().expect("b'' not available");
        let alpha = transcript.alpha.as_ref().expect("alpha not available");
        let beta = transcript.beta.as_ref().expect("beta not available");
        let psi = transcript.psi.as_ref().expect("psi not available");

        // Compute b = sum_{k in [K]} alpha_k * b^(k) + sum_{k in [K']} beta_k * b''^(k)
        let mut b = R::zero();
        for k in 0..crs.num_constraints {
            b += alpha[k] * transcript.instance.quad_dot_prod_funcs[k].b;
        }
        for k in 0..crs.num_aggregs {
            b += beta[k] * b__[k];
        }
        // Compute a_ij = sum_{k in [K]} alpha_k * a_ij^(k) + sum_{k in [K']} beta_k * a_ij''^(k), where a_ij''^(k) = sum_{l in [L]} psi_l^(k) * a_ij'^(l)
        let mut A = Matrix::<R>::zeros(crs.r, crs.r);
        for k in 0..crs.num_constraints {
            if let Some(ref A_k) = transcript.instance.quad_dot_prod_funcs[k].A {
                A += A_k * alpha[k];
            }
        }
        for k in 0..crs.num_aggregs {
            let mut a_k__ = Matrix::<R>::zeros(crs.r, crs.r);
            for l in 0..crs.num_constant_constraints {
                if let Some(ref A_l_) = transcript.instance.ct_quad_dot_prod_funcs[l].A {
                    a_k__ += mul_matrix_basescalar(A_l_, psi[k][l]);
                }
            }
            A += a_k__ * beta[k];
        }
        let mut phis_next = vec![Vector::<R>::zeros(n_next); r_next];
        // Set phis for a_ij * g_ij
        let mut A_vec = Vec::<R>::with_capacity((crs.r * (crs.r + 1)).div_ceil(2));
        for i in 0..crs.r {
            for j in 0..i {
                A_vec.push(A[(i, j)] + A[(j, i)]);
            }
            A_vec.push(A[(i, i)]);
        }
        phis_next[2 * nu + tau..2 * nu + tau + gamma]
            .clone_from_slice(chunk_pad(&Vector::<R>::from_vec(A_vec), n_next).as_slice());

        // Set phis for 1 * h_ii, i.e.
        let mut indic = Vector::<R>::zeros(crs.r * (crs.r + 1).div_ceil(2));
        for i in 0..crs.r {
            indic[((i + 1) * (i + 2)) / 2] = R::one();
        }
        phis_next[2 * nu + tau + gamma..2 * nu + tau + gamma + eta]
            .clone_from_slice(chunk_pad(&indic, n_next).as_slice());

        quad_dot_prod_funcs_next.push(QuadDotProdFunction::<R>::new(
            Matrix::<R>::zeros(r_next, r_next),
            phis_next,
            b,
        ));
    }

    // Constraints for u_1
    let u_1 = transcript.u_1.as_ref().expect("u_1 not available");
    for l in 0..crs.k1 {
        let mut phis = vec![Vector::<R>::zeros(n_next); r_next];
        let mut B_flat = Vec::<R>::with_capacity(crs.r * crs.t1 * crs.n);
        for i in 0..crs.r {
            for k in 0..crs.t1 {
                B_flat.extend(crs.B[i][k].row(l).transpose().iter());
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

        quad_dot_prod_funcs_next.push(QuadDotProdFunction::<R>::new(
            Matrix::<R>::zeros(r_next, r_next),
            phis,
            u_1[l],
        ));
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
        assert_eq!(D_flat_split.len(), eta);
        phis[2 * nu + tau + gamma - 1..2 * nu + tau + gamma - 1 + eta]
            .clone_from_slice(D_flat_split.as_slice());

        quad_dot_prod_funcs_next.push(QuadDotProdFunction::<R>::new(
            Matrix::<R>::zeros(r_next, r_next),
            phis,
            u_2[l],
        ));
    }

    debug_assert_eq!(
        crs.next_crs.as_ref().unwrap().num_constraints,
        quad_dot_prod_funcs_next.len()
    );
    debug_assert_eq!(crs.next_crs.as_ref().unwrap().num_constant_constraints, 0);

    let instance_next = Instance::<R> {
        quad_dot_prod_funcs: quad_dot_prod_funcs_next,
        ct_quad_dot_prod_funcs: vec![],
    };

    let witness = if compute_witness {
        let (z, t, G, H) = (
            transcript.z.as_ref().unwrap(),
            transcript.t.as_ref().unwrap(),
            transcript.G.as_ref().unwrap(),
            transcript.H.as_ref().unwrap(),
        );
        let z_decomp = decompose_balanced_vec_polyring(z, crs.b, Some(2usize));
        assert_eq!(z_decomp.len(), 2);

        let v = concat(&[
            flatten_vec_vector(t).as_slice(),
            flatten_symmetric_matrix(G).as_slice(),
            flatten_symmetric_matrix(H).as_slice(),
        ]);
        let z_0_split = crate::util::split(&z_decomp[0], nu);
        let z_1_split = crate::util::split(&z_decomp[1], nu);
        let v_split = crate::util::split(&v, mu);
        Some(Witness::<R> {
            s: [z_0_split, z_1_split, v_split].concat(),
        })
    } else {
        None
    };

    let index_next = *crs.clone().next_crs.unwrap();
    (index_next, instance_next, witness)
}
