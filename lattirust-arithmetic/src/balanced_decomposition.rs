use std::iter::Sum;

use ark_ff::Field;
use nalgebra::{ClosedAdd, ClosedMul, Scalar};
use num_traits::{One, Zero};
use rounded_div::RoundedDiv;

use crate::linear_algebra::Matrix;
use crate::linear_algebra::{RowVector, Vector};
use crate::ring::{ConvertibleRing, PolyRing, SignedRepresentative};

/// Returns the maximum number of terms in the balanced decomposition in basis `b` of any `x` with $|\textt{x}| \leq \textt{max}$.
pub fn balanced_decomposition_max_length(b: u128, max: u128) -> usize {
    if max == 0 {
        0
    } else {
        (max as f64).log(b as f64).floor() as usize + 1 + 1 // +1 because we are using balanced decomposition
    }
}

/// Given a vector of vectors `v`, pads each row to the same length $l = \max_i \texttt{v}\[i\]\texttt{.len()}$ and transposes the result. The output is a Vec of Vec of dimensionts `l` times `v.len()`.
pub fn pad_and_transpose<F: Copy + Zero>(mut v: Vec<Vec<F>>) -> Vec<Vec<F>> {
    let rows = v.len();
    let cols = v.iter().map(|d_i| d_i.len()).max().unwrap();
    // Pad each row to the same length `cols
    for row in &mut v {
        row.resize(cols, F::zero());
    }

    // Reshape as cols x rows
    (0..cols)
        .map(|col| (0..rows).map(|row| v[row][col]).collect::<Vec<F>>())
        .collect()
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
pub fn decompose_balanced<R: ConvertibleRing>(
    v: &R,
    b: u128,
    padding_size: Option<usize>,
) -> Vec<R> {
    assert!(
        !b.is_zero() && !b.is_one(),
        "cannot decompose in basis 0 or 1"
    );
    // TODO: not sure if this really necessary, but having b be even allow for more efficient divisions/remainders
    assert_eq!(b % 2, 0, "decomposition basis must be even");

    let b_half_floor = b.div_euclid(2);
    let b = b as i128;
    let mut decomp_bal_signed = Vec::<i128>::new();
    let mut curr = Into::<SignedRepresentative>::into(*v).0;
    loop {
        let rem = curr % b; // rem = curr % b is in [-(b-1), (b-1)]

        // Ensure digit is in [-b/2, b/2]
        if rem.abs() as u128 <= b_half_floor {
            decomp_bal_signed.push(rem);
            curr = curr / b; // Rust integer division rounds towards zero
        } else {
            // The next element in the decomposition is sign(rem) * (|rem| - b)
            if rem < 0 {
                decomp_bal_signed.push(rem + b);
            } else {
                decomp_bal_signed.push(rem - b);
            }
            let carry = rem.rounded_div(b); // Round toward nearest integer, not towards 0
            curr = (curr / b) + carry;
        }

        if curr.is_zero() {
            break;
        }
    }

    if let Some(padding_size) = padding_size {
        assert!(
            decomp_bal_signed.len() <= padding_size,
            "padding_size = {} must be at least decomp_bal.len() = {}",
            padding_size,
            decomp_bal_signed.len(),
        );
        decomp_bal_signed.resize(padding_size, 0);
    } else {
        // Strip trailing zeros
        while decomp_bal_signed.last() == Some(&0) {
            decomp_bal_signed.pop();
        }
    }

    decomp_bal_signed
        .into_iter()
        .map(|x| Into::<R>::into(SignedRepresentative(x)))
        .collect::<Vec<R>>()
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
pub fn decompose_balanced_vec<F: ConvertibleRing>(
    v: &[F],
    b: u128,
    padding_size: Option<usize>,
) -> Vec<Vec<F>> {
    let decomp: Vec<Vec<F>> = v
        .iter()
        .map(|v_i| decompose_balanced(v_i, b, padding_size))
        .collect(); // v.len() x decomp_size
    pad_and_transpose(decomp) // decomp_size x v.len()
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
pub fn decompose_balanced_polyring<R: PolyRing>(
    v: &R,
    b: u128,
    padding_size: Option<usize>,
) -> Vec<R> {
    decompose_balanced_vec::<R::BaseRing>(v.coeffs().as_slice(), b, padding_size)
        .into_iter()
        .map(|v_i| R::from(v_i))
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
pub fn decompose_balanced_vec_polyring<R: PolyRing>(
    v: &Vector<R>,
    b: u128,
    padding_size: Option<usize>,
) -> Vec<Vector<R>> {
    let decomp: Vec<Vec<R>> = v
        .as_slice()
        .iter()
        .map(|ring_elem| decompose_balanced_polyring(ring_elem, b, padding_size))
        .collect(); // v.len() x decomp_size
    pad_and_transpose(decomp)
        .into_iter()
        .map(|v_i| Vector::from(v_i))
        .collect() // decomp_size x v.len()
}

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
pub fn decompose_matrix<F: ConvertibleRing>(
    mat: &Matrix<F>,
    decomposition_basis: u128,
    decomposition_length: usize,
) -> Matrix<F> {
    Matrix::<F>::from_rows(
        mat.row_iter()
            .map(|s_i| {
                RowVector::<F>::from(
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

pub fn recompose<A, B>(v: &Vec<A>, b: B) -> A
where
    A: std::ops::Mul<B, Output = A> + Copy + Sum,
    B: Field,
{
    v.iter()
        .enumerate()
        .map(|(i, v_i)| *v_i * b.pow([i as u64]))
        .sum()
}

/// Given a m x n*k matrix `mat` decomposed in basis b and a slice \[1, b, ..., b^(k-1)] `powers_of_basis`, returns the m x n recomposed matrix.
pub fn recompose_matrix<F: Scalar + ClosedAdd + ClosedMul + Zero>(
    mat: &Matrix<F>,
    powers_of_basis: &[F],
) -> Matrix<F> {
    let (m, nk) = (mat.nrows(), mat.ncols());
    let k = powers_of_basis.len();
    debug_assert!(nk % k == 0, "matrix `mat` to be recomposed should have dimensions m x nk, where k is the length of `powers_of_basis`, but its number of columns is not a multiple of k");
    let n = nk / k;
    let pows = Vector::<F>::from_slice(powers_of_basis).transpose();
    Matrix::<F>::from_fn(m, n, |i, j| {
        mat.row(i).columns_range(j * k..j * k + k).dot(&pows)
    })
}

#[cfg(test)]
mod tests {
    use ark_std::test_rng;

    use crate::ntt::ntt_modulus;
    use crate::ring::pow2_cyclotomic_poly_ring_ntt::Pow2CyclotomicPolyRingNTT;
    use crate::ring::Zq;

    use super::*;

    const N: usize = 128;
    const Q: u64 = ntt_modulus::<N>(16);
    const VEC_LENGTH: usize = 32;
    const BASIS_TEST_RANGE: [u128; 5] = [2, 4, 8, 16, 32];

    type R = Zq<Q>;
    type PolyR = Pow2CyclotomicPolyRingNTT<Q, N>;

    #[test]
    fn test_decompose_balanced() {
        let vs: Vec<R> = (0..Q).map(|v| R::from(v)).collect();
        for b in BASIS_TEST_RANGE {
            let b_half = R::from(b / 2);
            for v in &vs {
                let decomp = decompose_balanced(v, b, None);

                // Check that the decomposition length is correct
                let max = SignedRepresentative::from(*v).0.abs() as u128;
                let max_decomp_length = balanced_decomposition_max_length(b, max);
                assert!(decomp.len() <= max_decomp_length);

                // Check that all entries are smaller than b/2
                for v_i in &decomp {
                    assert!(*v_i <= b_half || *v_i >= -b_half);
                }

                // Check that the decomposition is correct
                assert_eq!(*v, recompose(&decomp, R::from(b)));
            }
        }
    }

    fn get_test_vec() -> Vec<R> {
        (0..(N as u64) * Q)
            .step_by((Q / (N as u64)) as usize)
            .map(|x| R::from(x))
            .collect()
    }

    #[test]
    fn test_decompose_balanced_vec() {
        let v = get_test_vec();
        for b in BASIS_TEST_RANGE {
            let b_half = b / 2;
            let decomp = decompose_balanced_vec(&v, b, None);

            // Check that the decomposition length is correct
            let max = v
                .iter()
                .map(|x| SignedRepresentative::from(*x).0.abs() as u128)
                .max()
                .unwrap();
            let max_decomp_length = balanced_decomposition_max_length(b, max);
            assert!(decomp.len() <= max_decomp_length);

            // Check that all entries are smaller than b/2 in absolute value
            for d_i in &decomp {
                for d_ij in d_i {
                    let s_ij = SignedRepresentative::from(*d_ij).0;
                    assert!(s_ij.abs() as u128 <= b_half);
                }
            }

            for i in 0..decomp.len() {
                // Check that the decomposition is correct
                let decomp_i = decomp.iter().map(|d_j| d_j[i]).collect::<Vec<_>>();
                assert_eq!(v[i], recompose(&decomp_i, R::from(b)));
            }
        }
    }

    #[test]
    fn test_decompose_balanced_polyring() {
        let v = PolyR::from(get_test_vec());
        for b in BASIS_TEST_RANGE {
            let b_half = b / 2;
            let decomp = decompose_balanced_polyring(&v, b, None);

            for d_i in &decomp {
                for d_ij in d_i.coeffs() {
                    let s_ij = SignedRepresentative::from(d_ij).0;
                    assert!(s_ij.abs() as u128 <= b_half);
                }
            }

            assert_eq!(v, recompose(&decomp, R::from(b)));
        }
    }

    #[test]
    fn test_decompose_balanced_vec_polyring() {
        let v = Vector::<PolyR>::from_fn(VEC_LENGTH, |i, _| {
            PolyR::from(
                (0..(N as u64) * Q)
                    .step_by((Q / (N as u64)) as usize)
                    .map(|x| R::from(x + i as u64))
                    .collect::<Vec<_>>(),
            )
        });
        for b in BASIS_TEST_RANGE {
            let b_half = R::from(b / 2);
            let decomp = decompose_balanced_vec_polyring::<PolyR>(&v, b, None);

            for v_i in &decomp {
                for v_ij in v_i.as_slice() {
                    for v_ijk in v_ij.coeffs() {
                        assert!(v_ijk <= b_half || v_ijk >= -b_half);
                    }
                }
            }

            let mut recomposed = Vector::<PolyR>::zeros(v.len());
            for (i, v_i) in decomp.iter().enumerate() {
                recomposed += v_i * PolyR::from_scalar(R::from(b).pow(&[i as u64]));
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

            let decomp_size = balanced_decomposition_max_length(b, (Q - 1) as u128);
            let decomp = decompose_matrix(&mat, b, decomp_size);
            assert_eq!(decomp.nrows(), mat.nrows());
            assert_eq!(decomp.ncols(), mat.ncols() * decomp_size);
            for d_ij in decomp.iter() {
                assert!(Into::<SignedRepresentative>::into(*d_ij).0.abs() <= b_half as i128);
            }

            let pows_b = (0..decomp_size)
                .map(|i| R::from(b.pow(i as u32)))
                .collect::<Vec<_>>();
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
}
