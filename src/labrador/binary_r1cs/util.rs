#![allow(non_snake_case)]

use ark_ff::Field;
use num_traits::{One, Zero};

use crate::labrador::binary_r1cs::prover::{BinaryR1CSCRS, BinaryR1CSInstance, Z2};
use crate::lattice_arithmetic::matrix::{Matrix, Vector};
use crate::lattice_arithmetic::poly_ring::PolyRing;
use crate::lattice_arithmetic::ring::Ring;
use crate::relations::labrador::principal_relation::{ConstantQuadDotProdFunction, PrincipalRelation, QuadDotProdFunction};

pub const SECPARAM: usize = 128;

// fn inner_prod<F: Field>(row: &[(F, usize)], witness: &[F]) -> F {
//     let mut acc = F::zero();
//     for &(ref coeff, i) in row {
//         acc += &(if coeff.is_one() { witness[i] } else { witness[i] * coeff });
//     }
//     acc
// }

// pub fn mul_sparse_mat_vec<R: Field>(matrix: &Vec<Vec<(R, usize)>>, vec: &Vec<R>) -> Vec<R> {
//     ark_std::cfg_iter!(matrix)
//         .map(|row| inner_prod(row, vec))
//         .collect::<Vec<_>>()
// }

pub fn Z2_to_R_vec<R: Ring>(vec: &Vec<Z2>) -> Vec<R> {
    vec.iter().map(|x| if x.is_zero() { R::zero() } else { R::one() }).collect()
}

/// Reinterprets a vector of k = k' * d binary coefficients as k' vectors of d binary coefficients, represented as a vector of k' elements of the polynomial ring R with dimension d.
pub fn lift<R: PolyRing>(vec: &Vector<Z2>) -> Vector<R> {
    let d = R::dimension();
    assert_eq!(vec.len() % d, 0, "vector length {} must be multiple of dimension {}", vec.len(), d);
    let coeffs = vec.as_slice().chunks(d).map(|chunk| R::from(Z2_to_R_vec::<R::BaseRing>(&chunk.to_vec()))).collect();
    Vector::<R>::from_vec(coeffs)
}

/// Upcast an element in Z2 to an element in R
pub fn embed<R: Ring>(x: Z2) -> R {
    if x.is_zero() { R::zero() } else { R::one() }
}

pub fn basis_vector<R: PolyRing>(i: usize, n: usize) -> Vector<R> {
    assert!(i < n, "i = {} must be less than n = {}", i, n);
    let mut coeffs = vec![R::zero(); n];
    coeffs[i] = R::one();
    Vector::<R>::from_vec(coeffs)
}

// pub fn to_sparse_matrix(mat: &ark_relations::r1cs::Matrix<F2>) -> nalgebra_sparse::CsrMatrix<F2> {
//     let mut coo_mat = nalgebra_sparse::CooMatrix::<F2>::new(mat.len(), mat[0].len());
//     for (i, row) in mat.iter().enumerate() {
//         for (coeff, j) in row.into_iter() {
//             if !coeff.is_zero() {
//                 coo_mat.push(i, *j, *coeff);
//             }
//         }
//     }
//     nalgebra_sparse::CsrMatrix::<F2>::from(&coo_mat)
// }

pub struct BinaryR1CSTranscript<R: PolyRing> {
    pub t: Vector<R>,
    pub alpha: Matrix<Z2>,
    pub beta: Matrix<Z2>,
    pub gamma: Matrix<Z2>,
    pub g: Vector<R::BaseRing>,
    pub delta: Matrix<Z2>, // Not technically part of the transcript, but computed by prover and verifier
}

/// Express the constraint <alpha_i, a> = 0 as a constraint on the polynomial <alphaR_i, a_R> = 0, where alphaR_i is the element of R such that the constant term of alphaR_i * a_R (as polynomial multiplication over R) is equal to <alpha_i, a>
fn embed_Zqlinear_Rqlinear<R: PolyRing>(alpha_i: &Vector<Z2>, k: usize, n_pr: usize) -> Vector::<R> {
    let mut phi_a_idx = Vec::<R>::with_capacity(n_pr);
    let d = R::dimension();
    let k_ = k.div_ceil(d);
    debug_assert_eq!(k, d * k_);

    for j in 0..k_ {
        // Embed alpha_i as an element alphaR_i of R such that the constant term of alphaR_i * a_R (as polynomial multiplication over R) is equal to <alpha_i, a>
        let mut coeffs = vec![R::BaseRing::zero(); d];
        coeffs[0] = embed::<R::BaseRing>(alpha_i[j * d].clone());
        for l in 1..d {
            coeffs[d - 1 - l] = -embed::<R::BaseRing>(alpha_i[j * d + l].clone());
        }
        phi_a_idx.push(R::from(coeffs));
    }
    Vector::<R>::from(phi_a_idx)
}

pub fn reduce<R: PolyRing>(crs: &BinaryR1CSCRS<R>, instance: &BinaryR1CSInstance, transcript: &BinaryR1CSTranscript<R>) -> PrincipalRelation<R> {
    let (A, B, C) = (&instance.A, &instance.B, &instance.C);
    let (k, n) = (A.nrows(), A.ncols());
    let d = R::dimension();
    assert_eq!(k, n, "the current implementation only support k = n"); // TODO: remove this restriction by splitting a,b,c or w into multiple vectors

    let r_pr: usize = 8;
    let n_pr = n.div_ceil(d);
    let [a_idx, b_idx, c_idx, w_idx, a_tilde_idx, b_tilde_idx, c_tilde_idx, w_tilde_idx] = [0, 1, 2, 3, 4, 5, 6, 7];
    let (t, alpha, beta, gamma, g, delta) = (&transcript.t, &transcript.alpha, &transcript.beta, &transcript.gamma, &transcript.g, &transcript.delta);


    // F_1 = {A_i * (a || b || c || w) = t_i}_{i in [m/d]}
    let mut quad_dot_prod_funcs = Vec::<QuadDotProdFunction::<R>>::with_capacity(t.len() + 3 * n_pr);
    for i in 0..t.len() {
        let mut phi = crs.A.row(i).transpose().as_slice().chunks(n_pr).map(|v| Vector::<R>::from_row_slice(v)).collect::<Vec<Vector<R>>>();
        phi.append(&mut vec![Vector::<R>::zeros(n_pr); r_pr / 2]); // pad with zeros for "tilde witnesses"
        quad_dot_prod_funcs.push(QuadDotProdFunction::<R>::new(Matrix::<R>::zeros(r_pr, r_pr), phi, t[i]));
    }

    // TODO: the paper claims that the following constraints should be expressd as a "constant" constraint, but I don't see how this is possible. Added as standard constraints instead
    // ã = sigma_{-1}(a) <=>
    // ã_0 = a_0 and -ã_{n-i} = a_i for i in [n] <=>
    // <e_0, a> - <e_0, ã> = 0 and <e_i, a> + <e_{n-i}, ã> = 0, where e_i denote the i-th standard basis vector
    for i in 0..n_pr {
        for (idx, tilde_idx) in [(a_idx, a_tilde_idx), (b_idx, b_tilde_idx), (c_idx, c_tilde_idx)] {
            let mut phis = vec![Vector::<R>::zeros(n_pr); r_pr];
            phis[idx] = basis_vector(i, n_pr);
            phis[tilde_idx] = if i == 0 {
                -basis_vector(0, n_pr)
            } else {
                basis_vector(n_pr - i, n_pr)
            };
            quad_dot_prod_funcs.push(QuadDotProdFunction::<R>::new(Matrix::<R>::zeros(r_pr, r_pr), phis, R::zero()));
        }
    }

    // F_2
    let mut ct_quad_dot_prod_funcs = Vec::<ConstantQuadDotProdFunction::<R>>::with_capacity(5 + SECPARAM);
    // <a, ã - 1> = 0 <=>
    // <a, ã> + <ã, a> - <2, a> = 0
    for (idx, tilde_idx) in [(a_idx, a_tilde_idx), (b_idx, b_tilde_idx), (c_idx, c_tilde_idx), (w_idx, w_tilde_idx)] {
        let mut A = Matrix::<R>::zeros(r_pr, r_pr);
        A[(idx, tilde_idx)] = R::one();
        A[(tilde_idx, idx)] = R::one();
        let mut phi = vec![Vector::<R>::zeros(n_pr); r_pr];
        phi[idx] = Vector::<R>::from_element(n_pr, -R::from(2u128));
        ct_quad_dot_prod_funcs.push(ConstantQuadDotProdFunction::<R>::new(A, phi, R::BaseRing::zero()));
    }

    // <a + b - 2c, ã + ~b - 2~c - 1> = 0 <=>
    // <a, ã> + <a, ~b> -2*<a, ~c> + <-1, a> +
    // <b, ã> + <b, ~b> -2*<b, ~c> + <-1, b> +
    // -2*<c, ã> -2*<c, ~b> +4*<c, ~c> + <2, c> = 0
    // => double everything to make sure A is symmetric
    let mut A = Matrix::<R>::zeros(r_pr, r_pr);
    let min_two = -R::from(2u128);
    let vals = [
        (a_idx, a_tilde_idx, R::one()),
        (a_idx, b_tilde_idx, R::one()),
        (a_idx, c_tilde_idx, min_two),
        (b_idx, a_tilde_idx, R::one()),
        (b_idx, b_tilde_idx, R::one()),
        (b_idx, c_tilde_idx, min_two),
        (c_idx, a_tilde_idx, min_two),
        (c_idx, b_tilde_idx, min_two),
        (c_idx, c_tilde_idx, R::from(4u128))
    ];
    for (i, j, v) in vals {
        debug_assert!(A[(i, j)].is_zero());
        debug_assert!(A[(j, i)].is_zero());
        A[(i, j)] = v;
        A[(j, i)] = v;
    }
    let mut phi = vec![Vector::<R>::zeros(n_pr); r_pr];
    phi[a_idx] = Vector::<R>::from_element(n_pr, min_two);
    phi[b_idx] = Vector::<R>::from_element(n_pr, min_two);
    phi[c_idx] = Vector::<R>::from_element(n_pr, R::from(4u128));
    ct_quad_dot_prod_funcs.push(ConstantQuadDotProdFunction::<R>::new(A, phi, R::BaseRing::zero()));

    for i in 0..SECPARAM {
        // Constrain <alpha_i, a_i> + <beta_i, b_i> + <gamma_i, c_i> - <delta_i, w_i> = g_i (over the constant coefficients)
        let mut phi = vec![Vector::<R>::zeros(n_pr); r_pr];
        phi[a_idx] = embed_Zqlinear_Rqlinear(&alpha.row(i).transpose(), k, n_pr);
        phi[b_idx] = embed_Zqlinear_Rqlinear(&beta.row(i).transpose(), k, n_pr);
        phi[c_idx] = embed_Zqlinear_Rqlinear(&gamma.row(i).transpose(), k, n_pr);
        phi[w_idx] = -embed_Zqlinear_Rqlinear(&delta.row(i).transpose(), k, n_pr);

        ct_quad_dot_prod_funcs.push(ConstantQuadDotProdFunction::<R>::new(Matrix::<R>::zeros(r_pr, r_pr), phi, g[i]));
    }

    PrincipalRelation::<R> {
        quad_dot_prod_funcs,
        ct_quad_dot_prod_funcs,
    }
}