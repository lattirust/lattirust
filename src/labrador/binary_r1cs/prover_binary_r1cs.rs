#![allow(non_snake_case)]

use anyhow::Error;
use ark_ff::{Field, Fp, MontBackend};
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

const SECPARAM: usize = 128;

#[derive(MontConfig)]
#[modulus = "2"]
#[generator = "1"]
pub struct F2Config;

pub type F2 = Fp<MontBackend<F2Config, 1>, 1>;
pub type Z2 = Zq<2>;

pub struct BinaryR1CSCRS<R: PolyRing> {
    pub A: Matrix<R>,
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
        }
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

/// Upcast an element in Z2 to an element in Zq
pub fn embed<R: PolyRing>(x: Z2) -> R::BaseRing {
    if x.is_zero() { R::BaseRing::zero() } else { R::BaseRing::one() }
}

pub fn basis_vector<R: PolyRing>(i: usize, n: usize) -> Vector<R> {
    assert!(i < n, "i = {} must be less than n = {}", i, n);
    let mut coeffs = vec![R::zero(); n];
    coeffs[i] = R::one();
    Vector::<R>::from_vec(coeffs)
}

pub fn prove_binary_r1cs<'a, R: PolyRing>(crs: BinaryR1CSCRS<R>, arthur: &'a mut LatticeArthur<R>, instance: &BinaryR1CSInstance, witness: &BinaryR1CSWitness) -> Result<&'a [u8], Error>
    where LabradorChallengeSet<R>: FromRandomBytes<R>, WeightedTernaryChallengeSet<R>: FromRandomBytes<R>
{
    let (A, B, C) = (&instance.A, &instance.B, &instance.C);
    let w = &witness.w;
    let (k, n) = (A.nrows(), A.ncols());
    let d = R::dimension();
    assert_eq!(k, n, "the current implementation only support k = n"); // TODO: remove this restriction by splitting a,b,c or w into multiple vectors

    // TODO: add statement
    arthur.ratchet()?;

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
    let (a_BR, b_BR, c_BR) = (a.map(embed::<R>), b.map(embed::<R>), c.map(embed::<R>));
    let (alpha_BR, beta_BR, gamma_BR) = (alpha.map(embed::<R>), beta.map(embed::<R>), gamma.map(embed::<R>));
    let (delta_BR, w_BR) = (delta.map(embed::<R>), w.map(embed::<R>));

    let g = &alpha_BR * &a_BR + &beta_BR * &b_BR + &gamma_BR * &c_BR - &delta_BR * &w_BR;
    arthur.absorb_vector(&g).unwrap();

    let a_tilde = R::sigma_vec(&a_R);
    let b_tilde = R::sigma_vec(&b_R);
    let c_tilde = R::sigma_vec(&c_R);
    let w_tilde = R::sigma_vec(&w_R);
    let [a_idx, b_idx, c_idx, w_idx, a_tilde_idx, b_tilde_idx, c_tilde_idx, w_tilde_idx] = [0, 1, 2, 3, 4, 5, 6, 7];

    let r_pr: usize = 8;
    let n_pr = n.div_ceil(d);

    // F_1 = {A_i * (a || b || c || w) = t_i}_{i in [m/d]}
    let mut quad_dot_prod_funcs = Vec::<QuadDotProdFunction::<R>>::with_capacity(t.len());
    for i in 0..t.len() {
        let phi = crs.A.row(i).transpose().as_slice().chunks(n_pr).map(|v| Vector::<R>::from_row_slice(v)).collect::<Vec<Vector<R>>>();
        quad_dot_prod_funcs.push(QuadDotProdFunction::<R> {
            A: Matrix::<R>::zeros(r_pr, r_pr),
            phi,
            b: t[i],
        });
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
            quad_dot_prod_funcs.push(QuadDotProdFunction::<R> {
                A: Matrix::<R>::zeros(r_pr, r_pr),
                phi: phis,
                b: R::zero(),
            });
        }
    }

    // F_2
    let mut ct_quad_dot_prod_funcs = Vec::<ConstantQuadDotProdFunction::<R>>::with_capacity(3 * k); // TODO
    // <a, ã - 1> = 0 <=>
    // <a, ã> + <ã, a> - <2, a> = 0
    for (idx, tilde_idx) in [(a_idx, a_tilde_idx), (b_idx, b_tilde_idx), (c_idx, c_tilde_idx), (w_idx, w_tilde_idx)] {
        let mut A = Matrix::<R>::zeros(r_pr, r_pr);
        A[(idx, tilde_idx)] = R::one();
        A[(tilde_idx, idx)] = R::one();
        let mut phi = vec![Vector::<R>::zeros(n_pr); r_pr];
        phi[idx] = Vector::<R>::from_element(n_pr, R::from(2u128));
        ct_quad_dot_prod_funcs.push(ConstantQuadDotProdFunction::<R> {
            A,
            phi,
            b: R::BaseRing::zero(),
        })
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
    ct_quad_dot_prod_funcs.push(ConstantQuadDotProdFunction::<R> {
        A,
        phi,
        b: R::BaseRing::zero(),
    });

    for i in 0..SECPARAM {
        let mut phi = vec![Vector::<R>::zeros(n_pr); r_pr];

        let k_ = k.div_ceil(d);
        debug_assert_eq!(k, d * k_);
        let mut phi_a_idx = Vec::<R>::with_capacity(n_pr);
        for j in 0..k_ {
            // Embed alpha_i as an element alphaR_i of R such that the constant term of alphaR_i * a_R (as polynomial multiplication over R) is equal to <alpha_i, a>
            let mut coeffs = vec![R::BaseRing::zero(); d];
            coeffs[0] = embed::<R>(alpha[(i, j * d)].clone());
            for l in 1..d {
                coeffs[d - 1 - l] = -embed::<R>(alpha[(i, j * d + l)].clone());
            }
            phi_a_idx.push(R::from(coeffs));
        }
        phi[a_idx] = Vector::<R>::from_vec(phi_a_idx);

        ct_quad_dot_prod_funcs.push(ConstantQuadDotProdFunction::<R> {
            A: Matrix::<R>::zeros(r_pr, r_pr),
            phi,
            b: g[i],
        })
    }


    let norm_bound = (R::BaseRing::modulus() as f64).sqrt();
    let instance_pr = PrincipalRelation::<R> {
        r: r_pr,
        n: n_pr,
        norm_bound,
        quad_dot_prod_funcs,
        ct_quad_dot_prod_funcs,
    };

    let witness_pr = Witness::<R> {
        s: vec![a_R, b_R, c_R, w_R, a_tilde, b_tilde, c_tilde, w_tilde] // see definition of indices above
    };

    let basis = R::BaseRing::from(2u128);
    let crs_pr = CommonReferenceString::<R>::new(r_pr, n_pr, R::dimension(), norm_bound, 0, 0, 0, basis);

    prove_principal_relation(arthur, &instance_pr, &witness_pr, &crs_pr)
}

#[cfg(test)]
mod test {
    use nimue::hash::Keccak;
    use crate::labrador::iopattern::LabradorIOPattern;
    use crate::lattice_arithmetic::pow2_cyclotomic_poly_ring::Pow2CyclotomicPolyRing;
    use crate::nimue::iopattern::LatticeIOPattern;
    use super::*;


    const Q: u64 = 4294967291;
    // 2^32-5, prime
    const D: usize = 64;

    type BR = Zq<Q>;
    type R = Pow2CyclotomicPolyRing<BR, D>;

    #[test]
    fn test_prove_binary_r1cs() {
        let A = Matrix::<Z2>::identity(D, D);
        let B = Matrix::<Z2>::identity(D, D);
        let C = Matrix::<Z2>::identity(D, D);

        let crs = BinaryR1CSCRS::<R>::new(D, D);
        let instance = BinaryR1CSInstance { A, B, C };
        let witness = BinaryR1CSWitness { w: Vector::<Z2>::from_element(D, Z2::one()) };

        let io = LatticeIOPattern::<R, Keccak>::new("labrador_binaryr1cs")
            .ratchet()
            .labrador_binaryr1cs_io(&instance, &crs);
        let mut arthur = io.to_arthur();
        prove_binary_r1cs(crs, &mut arthur, &instance, &witness).unwrap();
    }
}