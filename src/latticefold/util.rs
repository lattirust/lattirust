use std::rc::Rc;
use std::vec;

use ark_ff::{Field, Fp4Config, MontConfig};
use ark_ff::fields::Fp2Config;
use ark_ff::fields::models::quadratic_extension::QuadExtConfig;
use ark_linear_sumcheck::ml_sumcheck::data_structures::ListOfProductsOfPolynomials;
use ark_poly::{DenseMultilinearExtension, MultilinearExtension};
use ark_std::iterable::Iterable;
use num_traits::Zero;

use crate::lattice_arithmetic::matrix::{Matrix, Vector};
use crate::lattice_arithmetic::ntt::NTT;
use crate::lattice_arithmetic::poly_ring::PolyRing;
use crate::lattice_arithmetic::pow2_cyclotomic_poly_ring_ntt::Pow2CyclotomicPolyRingNTT;
use crate::lattice_arithmetic::ring::Fq;
use crate::relations::{ccs, cm_ccs};
use crate::sumcheck::boolean_hypercube::BooleanHypercube;

const TAU: usize = 4;
const Q: u64 = ((1u128 << 64) - 59u128) as u64;
const D: usize = 64;  //TODO: check?

pub type F = Fq<Q>;
pub type R = Pow2CyclotomicPolyRingNTT<Q, D>;

// TODO: implement extension fields F_q^tau and partial NTT

pub struct Crs<R: PolyRing> {
    pub m: usize,
    pub b: usize,
    pub B: usize,
    pub k: usize,
    // Invariant: b^k = B
    _marker: std::marker::PhantomData<R>,
}


pub fn transpose<T>(v: &Vec<Vec<T>>) -> Vec<Vec<T>> where T: Clone {
    let rows = v.len();
    let cols = v[0].len();
    debug_assert!(v.iter().all(|row| row.len() == cols));
    (0..cols).map(|col| {
        (0..rows)
            .map(|row| v[row][col].clone())
            .collect()
    }).collect()
}

// v in (R_q = Z_q[X]/(X^d+1))^m -> (f_i)_{i in [d]}, where f_i: {0,1}^log(m) -> R_q is the MLE of the i-th coefficients of v in NTT representation
pub fn mle_vec<R: PolyRing, const Q: u64, const N: usize>(v: &Vector<R>) -> Vec<DenseMultilinearExtension<Fq<Q>>>
    where R: NTT<Q, N>
{
    let m = v.len();
    assert!(m.is_power_of_two());
    let log_m = m.next_power_of_two().ilog2() as usize;
    // v_coeffs[i][j] = v[j].ntt_coeffs()[i] holds the i-th NTT coefficient of v[j]
    let v_coeffs = transpose(&v.iter().map(|v_i| v_i.ntt_coeffs()).collect::<Vec<_>>());
    v_coeffs.into_iter().map(|v_i| DenseMultilinearExtension::from_evaluations_vec(log_m, v_i)).collect()
}

/// M in (R_q = Z_q[X]/(X^d+1))^{m x n} -> (f_i)_{i in [d]}, where f_i: {0,1}^{log(m) + log(n)} -> R_q is the MLE of the i-th coefficients of M in NTT representation
/// For R_i, C_j denoting row and column indices, respectively, f_i will use variables [R_1, ..., R_logm, C_1, ..., C_logn] if column_major is True, and [C_1, ..., C_logn, R_1, ..., R_logm] otherwise
pub fn mle_mat<R: PolyRing, const Q: u64, const N: usize>(mat: &Matrix<R>, column_major: bool) -> Vec<DenseMultilinearExtension<Fq<Q>>>
    where R: NTT<Q, N>
{
    let (m, n) = mat.shape();
    assert!(m.is_power_of_two() && n.is_power_of_two());
    if column_major {
        mle_vec(&Vector::<R>::from_row_slice(mat.as_slice()))
    } else {
        mle_vec(&Vector::<R>::from_row_slice(mat.transpose().as_slice()))
    }
}

/// On input b in F^n, return eq in F[X_1,...,X_n], eq(b, X_1,...,X_n) = prod_{i in [n]} (1 - b_i) * (1 - X_i) +  b_i * X_i
pub fn eq<F: Field>(b: &Vector<F>) -> DenseMultilinearExtension<F> {
    // TODO: can we make this more efficient?
    // TODO: is it better to return a vec of all sub-products for use later on?
    let n = b.len();
    let evals = BooleanHypercube::new(n as u8).map(|x: Vec<F>|
        (0..n).map(|i| F::one() - x[i] + b[i]).product()
    ).collect();
    DenseMultilinearExtension::from_evaluations_vec(n, evals)
}

pub fn scale_mle(mle: DenseMultilinearExtension<F>, scalar: F) -> DenseMultilinearExtension<F> {
    DenseMultilinearExtension::from_evaluations_vec(mle.num_vars, mle.evaluations.into_iter().map(|e| e * scalar).collect())
}

/// Return the MLE corresponding to the multiplication of f and g. The caller must ensure that f and g are defined over the same variables.
fn mul_mle(f: &DenseMultilinearExtension<F>, g: &DenseMultilinearExtension<F>) -> DenseMultilinearExtension<F> {
    assert_eq!(f.num_vars, g.num_vars);
    DenseMultilinearExtension::from_evaluations_vec(f.num_vars, f.evaluations.iter().zip(&g.evaluations).map(|(f_i, g_i)| f_i * g_i).collect())
}

/// Return sum_{b in {0,1}^n_c} mle[M_j](X, b) * mle[z_ccs](b)
fn linearization_sumcheck_poly_inner(mat: &Matrix<R>, z: &Vector<R>) -> Vec<DenseMultilinearExtension<F>> {
    let log_nc = mat.ncols().next_power_of_two().ilog2();
    let mle_mat = mle_mat::<R, Q, D>(mat, false); // Ensure we can partially evaluate on column indices without relabelling
    let mle_z = mle_vec::<R, Q, D>(&z);
    let mut mle_res = vec![DenseMultilinearExtension::zero(); D];
    for k in 0..D {
        for b in BooleanHypercube::new(log_nc as u8) {
            mle_res[k] += scale_mle(mle_mat[k].fix_variables(&b), mle_z[k].evaluate(&b).unwrap());
        }
    }
    mle_res
}

pub fn linearization_sumcheck_poly(crs: &cm_ccs::Crs<R>, x: &cm_ccs::Instance<R>, w: &cm_ccs::Witness<R>, beta: &Vector<F>) -> Vec<ListOfProductsOfPolynomials<F>> {
    let log_m = crs.crs_cm.m.next_power_of_two().ilog2() as usize;
    let mut list = vec![ListOfProductsOfPolynomials::<F>::new(log_m); D];

    let eq = eq(&beta);
    let z_ccs = ccs::concat(&x.x_ccs, &w.w_ccs);

    // In each iteration of the outer loop, terms[k] = [ eq(X) ] ++ [ sum_b mle[M_j](X, b) * mle[z](b) ]_{j in S_i}
    let mut terms = vec![Vec::with_capacity(crs.crs_ccs.deg + 1); D];

    for i in 0..crs.crs_ccs.n_matrices {
        terms.clear();
        for k in 0..D {
            terms[k].push(eq.clone());
        }
        for j in crs.crs_ccs.multisets[i].iter() {
            let mut ps = linearization_sumcheck_poly_inner(&crs.crs_ccs.matrices[*j], &z_ccs);

            for (k, p_k) in ps.into_iter().enumerate() {
                terms[k].push(p_k);
            }
        }
        for k in 0..D {
            let c_i = Fq::<Q>::from(crs.crs_ccs.scalars[i].ntt_coeffs()[k]);
            list[k].add_product(
                terms[k].iter().map(|p| Rc::new(p.clone())),
                c_i);
        }
    }

    list
}

#[cfg(test)]
mod tests {
    use ark_poly::MultilinearExtension;
    use ark_std::UniformRand;

    use crate::lattice_arithmetic::matrix::{sample_uniform_mat, Vector};
    use crate::lattice_arithmetic::pow2_cyclotomic_poly_ring_ntt::Pow2CyclotomicPolyRingNTT;
    use crate::lattice_arithmetic::ring::Fq;
    use crate::sumcheck::boolean_hypercube::index;

    use super::*;

    const Q: u64 = 2u64.pow(16) + 1;

    type F = Fq::<Q>;

    const M: usize = 1 << 1;
    const N: usize = 1 << 2;
    const D: usize = 64;

    type R = Pow2CyclotomicPolyRingNTT<Q, D>;

    #[test]
    pub fn test_mle_vec() {
        let rng = &mut ark_std::test_rng();
        let v = Vector::from_vec((0..M).map(|x| R::rand(rng)).collect());
        let mle = mle_vec::<R, Q, D>(&v);
        assert_eq!(mle.len(), D);
        for k in 0..D {
            for i in 0..M {
                assert_eq!(mle[k].evaluate(&index(i as u128, M as u128)).unwrap(), v[i].ntt_coeffs()[k]);
            }
        }
    }

    #[test]
    pub fn test_mle_mat() {
        let rng = &mut ark_std::test_rng();
        let mat = sample_uniform_mat(M, N, rng);
        let mle = mle_mat::<R, Q, D>(&mat, true);
        assert_eq!(mle.len(), D);
        for k in 0..D {
            for i in 0..M {
                let i_bin = index(i as u128, M as u128);
                let mle_ki = mle[k].fix_variables(&i_bin);
                for j in 0..N {
                    let j_bin = index(j as u128, N as u128);
                    let mut ij_bin = i_bin.clone();
                    ij_bin.extend_from_slice(&j_bin);
                    assert_eq!(mle_ki.evaluate(&j_bin).unwrap(), mat[(i, j)].ntt_coeffs()[k]);
                    assert_eq!(mle[k].evaluate(&ij_bin).unwrap(), mat[(i, j)].ntt_coeffs()[k]);
                }
            }
        }
    }

    #[test]
    pub fn test_sumcheck() {
        let rng = &mut ark_std::test_rng();
        let v = Vector::from_vec((0..M).map(|x| R::rand(rng)).collect());
        let mle = mle_vec::<R, Q, D>(&v);
        let log_m = M.next_power_of_two().ilog2() as usize;
    }
}