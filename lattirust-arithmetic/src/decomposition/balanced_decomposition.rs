use std::iter::Sum;
use std::ops::Mul;

use nalgebra::Scalar;
use num_bigint::BigUint;
use num_traits::{One, Signed, ToPrimitive, Zero};
use rayon::prelude::*;
use rounded_div::RoundedDiv;

use crate::decomposition::pad_zeros;
use crate::linear_algebra::{
    ClosedAddAssign, ClosedMulAssign, Matrix, RowVector, SymmetricMatrix, Vector,
};
use crate::ring::representatives::WithSignedRepresentative;
use crate::ring::{PolyRing, Ring};

use super::{
    pad_and_transpose, DecompositionFriendlySignedRepresentative,
};

/// Returns the maximum number of terms in the balanced decomposition in basis `b` of any `x` with $|\textt{x}| \leq \textt{max}$.
pub fn balanced_decomposition_max_length(b: u128, max: BigUint) -> usize {
    if max.is_zero() {
        0
    } else {
        max.to_f64().unwrap().log(b as f64).floor() as usize + 1 + 1 // +1 because we are using balanced decomposition
    }
}

/// Returns the balanced decomposition of a slice as a Vec of Vecs.
///
/// # Arguments
/// * `v`: input element
/// * `b`: basis for the decomposition, must be even
/// * `padding_size`: indicates whether the output should be padded with zeros to a specified length `k` if `padding_size` is `Some(k)`, or if it should be padded to the largest decomposition length required for `v` if `padding_size` is `None`
///
/// # Output
/// Returns `d`, the decomposition in basis `b` as a Vec of size `decomp_size`, i.e.,
/// $\texttt{v}\[i\] = \sum_{j \in \[k\]} \texttt{b}^j \texttt{d}\[j\]$ and $|\texttt{d}\[j\]| \leq \left\lfloor\frac{\texttt{b}}{2}\right\rfloor$.
pub fn decompose_balanced<R>(v: &R, basis: u128, padding_size: Option<usize>) -> Vec<R>
where
    R: Ring + WithSignedRepresentative,
    R::SignedRepresentative: DecompositionFriendlySignedRepresentative,
{
    assert!(
        !basis.is_zero() && !basis.is_one(),
        "cannot decompose in basis 0 or 1"
    );
    // TODO: not sure if this really necessary, but having b be even allow for more efficient divisions/remainders
    assert_eq!(basis % 2, 0, "decomposition basis must be even");

    let basis_half_floor = basis.div_euclid(2);
    let b = R::SignedRepresentative::try_from(basis as i128).unwrap();
    let b_half_floor = R::SignedRepresentative::try_from(basis_half_floor as i128).unwrap();

    let mut decomp_bal_signed = Vec::<R>::new();
    let mut curr = Into::<R::SignedRepresentative>::into(*v);
    loop {
        let rem = curr.clone() % b.clone(); // rem = curr % b is in [-(b-1), (b-1)]

        // Ensure digit is in [-b/2, b/2]
        if rem.abs() <= b_half_floor {
            decomp_bal_signed.push(R::try_from(rem).unwrap());
            curr = curr.clone() / b.clone(); // Rust integer division rounds towards zero
        } else {
            // The next element in the decomposition is sign(rem) * (|rem| - b)
            let digit = if rem.is_negative() {
                rem.clone() + b.clone()
            } else {
                rem.clone() - b.clone()
            };
            decomp_bal_signed.push(R::try_from(digit).unwrap());
            let carry = rem.rounded_div(b.clone()); // Round toward nearest integer, not towards 0
            curr = curr / b.clone() + carry;
        }

        if curr.is_zero() {
            break;
        }
    }

    pad_zeros(&mut decomp_bal_signed, padding_size);

    decomp_bal_signed
}

/// Returns the balanced decomposition of a slice as a Vec of Vecs.
///
/// # Arguments
/// * `v`: input slice, of length `l`
/// * `b`: basis for the decomposition, must be even
/// * `padding_size`: indicates whether the output should be padded with zeros to a specified length `k` if `padding_size` is `Some(k)`, or if it should be padded to the largest decomposition length required for `v` if `padding_size` is `None`
///
/// # Output
/// Returns `d` the decomposition in basis `b` as a Vec of size `decomp_size`, with each item being a Vec of length `l`, i.e.,
/// for all $i \in \[l\]: \texttt{v}\[i\] = \sum_{j \in \[k\]} \texttt{b}^j \texttt{d}\[i\]\[j\]$ and $|\texttt{d}\[i\]\[j\]| \leq \left\lfloor\frac{\texttt{b}}{2}\right\rfloor$.
pub fn decompose_balanced_vec<R>(v: &[R], b: u128, padding_size: Option<usize>) -> Vec<Vec<R>>
where
    R: Ring + WithSignedRepresentative,
    R::SignedRepresentative: DecompositionFriendlySignedRepresentative,
{
    let decomp: Vec<Vec<R>> = v
        .par_iter()
        .map(|v_i| decompose_balanced(v_i, b, padding_size))
        .collect(); // v.len() x decomp_size
    pad_and_transpose(decomp, padding_size) // decomp_size x v.len()
}

/// Returns the balanced decomposition of a [`PolyRing`] element as a Vec of [`PolyRing`] elements.
///
/// # Arguments
/// * `v`: `PolyRing` element to be decomposed
/// * `b`: basis for the decomposition, must be even
/// * `padding_size`: indicates whether the output should be padded with zeros to a specified length `k` if `padding_size` is `Some(k)`, or if it should be padded to the largest decomposition length required for `v` if `padding_size` is `None`
///
/// # Output
/// Returns `d` the decomposition in basis `b` as a Vec of size `decomp_size`, i.e.,
/// for all $\texttt{v} = \sum_{j \in \[k\]} \texttt{b}^j \texttt{d}\[j\]$ and $|\texttt{d}\[j\]| \leq \left\lfloor\frac{\texttt{b}}{2}\right\rfloor$.
pub fn decompose_balanced_polyring<PR: PolyRing>(
    v: &PR,
    b: u128,
    padding_size: Option<usize>,
) -> Vec<PR>
where
    PR::BaseRing: WithSignedRepresentative,
    <PR::BaseRing as WithSignedRepresentative>::SignedRepresentative:
        DecompositionFriendlySignedRepresentative,
{
    decompose_balanced_vec::<PR::BaseRing>(v.coefficients().as_slice(), b, padding_size)
        .into_par_iter()
        .map(|v_i| PR::from(v_i))
        .collect()
}

/// Returns the balanced decomposition of a slice of [`PolyRing`] elements as a Vec of [`Vector`] of [`PolyRing`] elements.
///
/// # Arguments
/// * `v`: input slice, of length `l`
/// * `b`: basis for the decomposition, must be even
/// * `padding_size`: indicates whether the output should be padded with zeros to a specified length `k` if `padding_size` is `Some(k)`, or if it should be padded to the largest decomposition length required for `v` if `padding_size` is `None`
///
/// # Output
/// Returns `d` the decomposition in basis `b` as a Vec of size `decomp_size`, with each item being a Vec of length `l`, i.e.,
/// for all $i \in \[l\]: \texttt{v}\[i\] = \sum_{j \in \[k\]} \texttt{b}^j \texttt{d}\[i\]\[j\]$ and $|\texttt{d}\[i\]\[j\]| \leq \left\lfloor\frac{\texttt{b}}{2}\right\rfloor$.
pub fn decompose_balanced_vec_polyring<PR: PolyRing>(
    v: &[PR],
    b: u128,
    padding_size: Option<usize>,
) -> Vec<Vector<PR>>
where
    PR::BaseRing: WithSignedRepresentative,
    <PR::BaseRing as WithSignedRepresentative>::SignedRepresentative:
        DecompositionFriendlySignedRepresentative,
{
    let decomp: Vec<Vec<PR>> = v
        .par_iter()
        .map(|ring_elem| decompose_balanced_polyring(ring_elem, b, padding_size))
        .collect(); // v.len() x decomp_size
    pad_and_transpose(decomp, padding_size)
        .into_par_iter()
        .map(Vector::from)
        .collect() // decomp_size x v.len()
}

// In: v: m x n; output: k x m x n where k is the decomposition length
pub fn decompose_vec_vector_dimfirst<PR: PolyRing>(
    v: &Vec<Vector<PR>>,
    b: u128,
    padding_size: Option<usize>,
) -> Vec<Vec<Vector<PR>>>
where
    PR::BaseRing: WithSignedRepresentative,
    <PR::BaseRing as WithSignedRepresentative>::SignedRepresentative:
        DecompositionFriendlySignedRepresentative,
{
    // // v: m x n
    // let decomp: Vec<Vec<Vector<PR>>> = v
    //     .as_slice()
    //     .par_iter()
    //     .map(|v_i| decompose_balanced_vec_polyring(v_i.as_slice(), b, padding_size))
    //     .collect(); // m x decomp_size x n

    // let b_ = PR::try_from(b).unwrap();
    // let m = v.len();
    // let n = v[0].len();
    // for i in 0..v.len() {
    //    let mut res_i = Vector::<PR>::zeros(n);
    //    for j in 0..decomp[i].len() {
    //         res_i += decomp[i][j].clone() * b_.pow([j as u64]);
    //     }
    //     assert_eq!(res_i, v[i]);
    // }
    // pad_and_transpose_vec_vec_vector(decomp, padding_size)
    //     .into_par_iter()
    //     .map(|v_i| v_i.into_par_iter().map(Vector::from).collect())
    //     .collect() // decomp_size x m x n

    let decomp: Vec<Vec<Vec<PR>>> = v.iter().map(
        |v_i| v_i.iter().map(
            |v_ij| decompose_balanced_polyring(v_ij, b, padding_size)
        ).collect()
    ).collect(); // m x n x k

    let m = v.len();
    let n = v[0].len();
    let decomp_len = decomp[0][0].len();
    let mut res = vec![vec![Vector::zeros(n); m]; decomp_len];
    for i in 0..m {
        for j in 0..n {
            let recomp = recompose(&decomp[i][j], PR::try_from(b).unwrap());
            assert_eq!(recomp, v[i][j]);
            for k in 0..decomp_len {
                res[k][i][j] = decomp[i][j][k];
            }
        }
    }

    let mut recomp = vec![Vector::<PR>::zeros(n); m];
    let b_ring = PR::try_from(b).unwrap();
    for k in 0..decomp_len {
        for i in 0..m {
            for j in 0..n {
                recomp[i][j] += res[k][i][j] * b_ring.pow([k as u64]);
            }
        }
    }
    assert_eq!(recomp, *v);

    res
}

// // In: v: m x n; output: k x m x n where k is the decomposition length
// pub fn decompose_matrix_dimfirst<PR: PolyRing>(
//     mat: &Matrix<PR>,
//     b: u128,
//     padding_size: Option<usize>,
// ) -> Vec<Matrix<PR>>
// where
//     PR::BaseRing: WithSignedRepresentative,
//     <PR::BaseRing as WithSignedRepresentative>::SignedRepresentative:
//         DecompositionFriendlySignedRepresentative,
// {
//     let decomp = mat.map(|v_ij| decompose_balanced_polyring(&v_ij, b, padding_size));
   
//     // Given a matrix of size m x n x k, output a vec of size k of m x n matrices 
//     let decomp_vec: Vec<Vec<PR>> = decomp
//         .as_slice()
//         .par_iter()
//         .map(|v_i| v_i.as_slice().to_vec())
//         .collect(); // m x n x decomp_size
//     pad_and_transpose_vec_vec_vector(decomp_vec, padding_size)
// }

// // v: k x m x n; output: m x n
// pub fn recompose_vec_vector_dimfirst<PR: PolyRing>(
//     v: &Vec<Vec<Vector<PR>>>,
//     basis: PR,
// ) -> Vec<Vector<PR>>
// where
//     PR::BaseRing: WithSignedRepresentative,
//     <PR::BaseRing as WithSignedRepresentative>::SignedRepresentative:
//         DecompositionFriendlySignedRepresentative,
// {
//     let vec_vector_vector: Vec<Matrix<PR>> = v.iter().map(|v_i| Vector::<Vector<PR>>::from(v_i.clone())).collect();
//     let recomp_vector_vector: Vector<Vector<PR>> = recompose(&vec_vector_vector, basis);
//     recomp_vector_vector.as_slice().iter().collect()
// }

/// Returns the balanced gadget decomposition of a [`Matrix`] of dimensions `m × n` as a matrix of dimensions `m × (k * n)`.
///
/// # Arguments
/// * `mat`: input matrix of dimensions `m × n`
/// * `b`: basis for the decomposition, must be even
/// * `padding_size`: indicates whether the decomposition length is the specified `k` if `padding_size` is `Some(k)`, or if `k` is the largest decomposition length required for `mat` if `padding_size` is `None`
///
/// # Output
/// Returns `d` the decomposition in basis `b` as a Matrix of dimensions `m × (k * n)`, i.e.,
/// $\texttt{mat} = \texttt{d} \times G_\texttt{n}$ where $G_\texttt{n} = I_\texttt{n} \otimes (1, \texttt{b}, \ldots, \texttt{b}^k) \in R^{\texttt{k}\texttt{n} \times \texttt{n}}$ and $|\texttt{d}\[i\]\[j\]| \leq \left\lfloor\frac{\texttt{b}}{2}\right\rfloor$.
pub fn decompose_matrix<R>(
    mat: &Matrix<R>,
    decomposition_basis: u128,
    decomposition_length: usize,
) -> Matrix<R>
where
    R: Ring + WithSignedRepresentative,
    R::SignedRepresentative: DecompositionFriendlySignedRepresentative,
{
    Matrix::<R>::from_rows(
        mat.par_row_iter()
            .map(|s_i| {
                RowVector::<R>::from(
                    s_i.map(|s_ij| {
                        decompose_balanced(&s_ij, decomposition_basis, Some(decomposition_length))
                    })
                    .as_slice()
                    .concat(),
                )
            })
            .collect::<Vec<_>>()
            .as_slice(),
    )
}

pub fn recompose<A, B>(vec: &Vec<A>, basis: B) -> A
where
    A: Mul<B, Output = A> + Copy + Sum + Send + Sync,
    B: Ring,
{
    vec.par_iter()
        .enumerate()
        .map(|(i, v_i)| *v_i * basis.pow([i as u64]))
        .sum()
}

/// Given a `m x n*k` matrix `mat` decomposed in basis b and a slice \[1, b, ..., b^(k-1)] `powers_of_basis`, returns the `m x n` recomposed matrix.
pub fn recompose_matrix<F: Scalar + ClosedAddAssign + ClosedMulAssign + Zero + Send + Sync>(
    mat: &Matrix<F>,
    powers_of_basis: &[F],
) -> Matrix<F> {
    let nk = mat.ncols();
    let k = powers_of_basis.len();
    debug_assert!(nk % k == 0, "matrix `mat` to be recomposed should have dimensions m x nk, where k is the length of `powers_of_basis`, but its number of columns is not a multiple of k");
    let n = nk / k;
    let pows = Vector::<F>::from_slice(powers_of_basis).transpose();
    let rows = mat
        .par_row_iter()
        .map(|r_i| {
            RowVector::from(
                (0..n)
                    .into_par_iter()
                    .map(|j| r_i.columns_range(j * k..j * k + k).dot(&pows))
                    .collect::<Vec<_>>(),
            )
        })
        .collect::<Vec<_>>();
    Matrix::<F>::from_rows(rows.as_slice())
}

/// Given a `n*d x n*d` symmetric matrix `mat` and a slice `\[1, b, ..., b^(d-1)\]` `powers_of_basis`, returns the `n x n` symmetric matrix corresponding to $G^T \textt{mat} G$, where $G = I_n \otimes (1, b, ..., b^(\textt{d}-1))$ is the gadget matrix of dimensions `n*d x n`.
pub fn recompose_left_right_symmetric_matrix<F: Clone + Sum + Send + Sync>(
    mat: &SymmetricMatrix<F>,
    powers_of_basis: &[F],
) -> SymmetricMatrix<F>
where
    for<'a> &'a F: Mul<&'a F, Output = F>,
{
    let (nd, d) = (mat.size(), powers_of_basis.len());
    assert_eq!(nd % d, 0);

    let n = nd / d;
    (0..n)
        .into_par_iter()
        .map(|i| {
            (0..=i)
                .into_par_iter()
                .map(|j| {
                    (0..nd)
                        .filter(|k| k / d == i)
                        .flat_map(|k| {
                            (0..nd).filter(|l| l / d == j).map(move |l| {
                                &mat[(k, l)] * &(&powers_of_basis[k % d] * &powers_of_basis[l % d])
                            })
                        })
                        .sum()
                })
                .collect()
        })
        .collect::<Vec<_>>()
        .into()
}

#[cfg(test)]
mod tests {
    use ark_std::test_rng;

    use crate::ring;
    use crate::ring::ntt::ntt_prime;
    use crate::ring::pow2_cyclotomic_poly_ring_ntt::Pow2CyclotomicPolyRingNTT;
    use crate::ring::util::powers_of_basis;
    use crate::ring::Zq1;
    use crate::traits::WithLinfNorm;

    use super::*;

    const N: usize = 128;
    const Q: u64 = ntt_prime::<N>(12);
    const VEC_LENGTH: usize = 32;
    const BASIS_TEST_RANGE: [u128; 5] = [2, 4, 8, 16, 32];

    type R = Zq1<Q>;
    type PolyR = Pow2CyclotomicPolyRingNTT<R, N>;

    #[test]
    fn test_decompose_balanced() {
        let vs: Vec<R> = (0..Q).map(|v| R::try_from(v).unwrap()).collect();
        for b in BASIS_TEST_RANGE {
            for v in &vs {
                let decomp = decompose_balanced(v, b, None);

                // Check that the decomposition length is correct
                let max = v.linf_norm();
                let max_decomp_length = balanced_decomposition_max_length(b, max);
                assert!(decomp.len() <= max_decomp_length);

                // Check that all entries are smaller than b/2
                for v_i in &decomp {
                    assert!(v_i.linf_norm() <= BigUint::from(b.div_floor(2)));
                }

                // Check that the decomposition is correct
                assert_eq!(*v, recompose(&decomp, R::try_from(b).unwrap()));
            }
        }
    }

    fn get_test_vec() -> [R; N] {
        core::array::from_fn(|i| R::try_from((i as u64 * Q) / (N as u64)).unwrap())
    }

    #[test]
    fn test_decompose_balanced_vec() {
        let v = get_test_vec();
        for b in BASIS_TEST_RANGE {
            let b_half = b / 2;
            let decomp = decompose_balanced_vec(&v, b, None);

            // Check that the decomposition length is correct
            let max = v.linf_norm();
            let max_decomp_length = balanced_decomposition_max_length(b, max);
            assert!(decomp.len() <= max_decomp_length);

            // Check that all entries are smaller than b/2 in absolute value
            for d_i in &decomp {
                for d_ij in d_i {
                    assert!(d_ij.linf_norm() <= b_half.into());
                }
            }

            for i in 0..decomp.len() {
                // Check that the decomposition is correct
                let decomp_i = decomp.iter().map(|d_j| d_j[i]).collect::<Vec<_>>();
                assert_eq!(v[i], recompose(&decomp_i, R::try_from(b).unwrap()));
            }
        }
    }

    #[test]
    fn test_decompose_balanced_polyring() {
        let v = PolyR::from_coefficient_array(get_test_vec());
        for b in BASIS_TEST_RANGE {
            let b_half = b / 2;
            let decomp: Vec<PolyR> = decompose_balanced_polyring(&v, b, None);

            for d_i in &decomp {
                for d_ij in d_i.coefficients() {
                    assert!(d_ij.linf_norm() <= b_half.into());
                }
            }

            assert_eq!(v, recompose(&decomp, R::try_from(b).unwrap()));
        }
    }

    #[test]
    fn test_decompose_balanced_vec_polyring() {
        let v = Vector::<PolyR>::from_fn(VEC_LENGTH, |i, _| {
            PolyR::try_from_coefficients(
                get_test_vec()
                    .into_iter()
                    .map(|v| v + R::try_from(i as u64).unwrap())
                    .collect::<Vec<_>>()
                    .as_slice(),
            ).unwrap()
        });
        for b in BASIS_TEST_RANGE {
            let b_half = b / 2;
            let decomp = decompose_balanced_vec_polyring::<PolyR>(&v.as_slice(), b, None);

            for v_i in &decomp {
                for v_ij in v_i.as_slice() {
                    for v_ijk in v_ij.coefficients() {
                        assert!(v_ijk.linf_norm() <= b_half.into());
                    }
                }
            }

            let mut recomposed = Vector::<PolyR>::zeros(v.len());
            for (i, v_i) in decomp.iter().enumerate() {
                recomposed +=
                    v_i * PolyR::from_scalar(Ring::pow(&R::try_from(b).unwrap(), [i as u64]));
            }
            assert_eq!(v, recomposed);
        }
    }

    #[test]
    pub fn test_decompose_balanced_matrix() {
        const NROWS: usize = 101;
        const NCOLS: usize = 42;

        let rng = &mut test_rng();
        let mat = Matrix::<R>::rand(NROWS, NCOLS, rng);

        for b in BASIS_TEST_RANGE {
            let b_half = b.div_floor(2);

            let decomp_size = balanced_decomposition_max_length(b, ((Q - 1) as u128).into());

            let decomp = decompose_matrix(&mat, b, decomp_size);
            assert_eq!(decomp.nrows(), mat.nrows());
            assert_eq!(decomp.ncols(), mat.ncols() * decomp_size);
            for d_ij in decomp.iter() {
                assert!(d_ij.as_signed_representative().0.abs() <= b_half.into());
            }

            let pows_b = ring::util::powers_of_basis(R::try_from(b).unwrap(), decomp_size);
            let recomp = recompose_matrix(&decomp, &pows_b);

            assert_eq!(mat, recomp);

            // Check that recompose is equivalent to right-multiplication by the gadget matrix
            let gadget_vector = Vector::<R>::from_slice(&pows_b);
            let gadget_matrix = Matrix::<R>::identity(NCOLS, NCOLS).kronecker(&gadget_vector);
            let mat_g = &decomp * &gadget_matrix;
            assert_eq!(recomp, mat_g);

            // Check that A * M == recompose(A * decompose(M))
            let mat_l = Matrix::<R>::rand(64, NROWS, rng);
            let lhs = &mat_l * &mat;
            let rhs = recompose_matrix(&(&mat_l * &decomp), &pows_b);
            assert_eq!(lhs, rhs);
        }
    }

    #[test]
    fn test_left_right_decompose_symmetric_matrix() {
        const SIZE: usize = 13;
        let rng = &mut test_rng();

        for b in BASIS_TEST_RANGE {
            let decomp_size = balanced_decomposition_max_length(b, (Q - 1).into());

            let pows_b = powers_of_basis(R::try_from(b).unwrap(), decomp_size);
            let mat = SymmetricMatrix::<R>::rand(SIZE * decomp_size, rng);

            let expected: SymmetricMatrix<R> = recompose_matrix(
                &recompose_matrix(&mat.clone().into(), pows_b.as_slice()).transpose(),
                pows_b.as_slice(),
            )
            .into();

            let actual = recompose_left_right_symmetric_matrix(&mat, pows_b.as_slice());
            assert_eq!(expected, actual);
        }
    }
}
