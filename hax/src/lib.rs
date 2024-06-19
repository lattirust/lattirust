use std::iter::Sum;
use std::ops::{Add, Mul};

use num_traits::{One, Pow};
use num_traits::Zero;
use rayon::prelude::*;
use rounded_div::RoundedDiv;

/// Returns the floor of the logarithm of `x` in base `base`.
pub fn floor_log(mut x: u128, base: u128) -> u128 {
    assert!(base > 1);
    assert!(x >= 1);
    let mut floor_log = 0;
    for _ in 0.. {
        if x >= base {
            x /= base;
            floor_log += 1;
        } else {
            // break;
        }
    }
    // while x >= base {
    //     x /= base;
    //     floor_log += 1;
    // }
    floor_log
}

/// Returns the maximum number of terms in the balanced decomposition in basis `b` of any `x` with $|\textt{x}| \leq \textt{max}$.
pub fn balanced_decomposition_max_length(b: u128, max: u128) -> usize {
    if max == 0 {
        0
    } else {
        floor_log(max, b) as usize + 1 + 1 // +1 because we are using balanced decomposition
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
pub fn decompose_balanced(
    v: &i128,
    b: u128,
    padding_size: Option<usize>,
    length: usize,
) -> Vec<i128> {
    assert!(
        !b.is_zero() && !b.is_one(),
        "cannot decompose in basis 0 or 1"
    );
    // TODO: not sure if this really necessary, but having b be even allow for more efficient divisions/remainders
    assert_eq!(b % 2, 0, "decomposition basis must be even");

    let b_half_floor = b.div_euclid(2);
    let b = b as i128;
    let mut decomp_bal_signed = Vec::<i128>::new();
    let mut curr = *v;
    for _ in 0..length {
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
    }

    if let Some(padding_size) = padding_size {
        assert!(
            decomp_bal_signed.len() <= padding_size,
            "padding_size = {} must be at least decomp_bal.len() = {}",
            padding_size,
            decomp_bal_signed.len(),
        );
        decomp_bal_signed.resize(padding_size, 0);
    }
    // else {
    //     // Strip trailing zeros
    //     while decomp_bal_signed.last() == Some(&0) {
    //         decomp_bal_signed.pop();
    //     }
    // }

    decomp_bal_signed
}

/*
pub const fn pad<F: Copy + Zero>(v: &mut Vec<Vec<F>>) {
    let cols = v.iter().map(|d_i| d_i.len()).max().unwrap();

    // Pad each row to the same length `cols
    for mut row in &v {
        row.resize(cols, F::zero());
    }
}

pub const fn transpose<F: Copy + Zero>(v: Vec<Vec<F>>) -> Vec<Vec<F>> {
    let rows = v.len();
    let cols = v[0].len();
    debug_assert!(v.iter().all(|v_i| v_i.len() == cols));

    // Reshape as cols x rows
    (0..cols)
        .map(|col| (0..rows).map(|row| v[row][col]).collect::<Vec<F>>())
        .collect()
}

*/

/// Given a vector of vectors `v`, pads each row to the same length $l = \max_i \texttt{v}\[i\]\texttt{.len()}$ and transposes the result. The output is a Vec of Vec of dimensionts `l` times `v.len()`.
/*pub fn pad_and_transpose<F: Copy + Zero>(mut v: Vec<Vec<F>>) -> Vec<Vec<F>> {
    pad()
}*/

pub fn recompose<A, B>(v: &Vec<A>, b: B) -> A
where
    A: std::ops::Mul<B, Output = A> + Copy + Sum + Send + Sync,
    B: Clone + Mul<B, Output = B> + Add<B, Output = B> + Pow<[u64; 1], Output = B> + Send + Sync,
{
    v.par_iter()
        .enumerate()
        .map(|(i, v_i)| *v_i * b.clone().pow([i as u64]))
        .sum()
}

// #[cfg(test)]
// mod tests {
//     use ark_std::test_rng;
//
//     use super::*;
//
//     const N: usize = 128;
//     const Q: u64 = (1<<16)-1;
//     const VEC_LENGTH: usize = 32;
//     const BASIS_TEST_RANGE: [u128; 5] = [2, 4, 8, 16, 32];
//
//     #[test]
//     fn test_decompose_balanced() {
//         let vs: Vec<R> = (0..Q).map(|v| R::from(v)).collect();
//         for b in BASIS_TEST_RANGE {
//             let b_half = R::from(b / 2);
//             for v in &vs {
//                 let decomp = decompose_balanced(v, b, None);
//
//                 // Check that the decomposition length is correct
//                 let max = SignedRepresentative::from(*v).0.abs() as u128;
//                 let max_decomp_length = balanced_decomposition_max_length(b, max);
//                 assert!(decomp.len() <= max_decomp_length);
//
//                 // Check that all entries are smaller than b/2
//                 for v_i in &decomp {
//                     assert!(*v_i <= b_half || *v_i >= -b_half);
//                 }
//
//                 // Check that the decomposition is correct
//                 assert_eq!(*v, recompose(&decomp, R::from(b)));
//             }
//         }
//     }
//
//     fn get_test_vec() -> Vec<R> {
//         (0..(N as u64) * Q)
//             .step_by((Q / (N as u64)) as usize)
//             .map(|x| R::from(x))
//             .collect()
//     }
//
//     #[test]
//     fn test_decompose_balanced_vec() {
//         let v = get_test_vec();
//         for b in BASIS_TEST_RANGE {
//             let b_half = b / 2;
//             let decomp = decompose_balanced_vec(&v, b, None);
//
//             // Check that the decomposition length is correct
//             let max = v
//                 .iter()
//                 .map(|x| SignedRepresentative::from(*x).0.abs() as u128)
//                 .max()
//                 .unwrap();
//             let max_decomp_length = balanced_decomposition_max_length(b, max);
//             assert!(decomp.len() <= max_decomp_length);
//
//             // Check that all entries are smaller than b/2 in absolute value
//             for d_i in &decomp {
//                 for d_ij in d_i {
//                     let s_ij = SignedRepresentative::from(*d_ij).0;
//                     assert!(s_ij.abs() as u128 <= b_half);
//                 }
//             }
//
//             for i in 0..decomp.len() {
//                 // Check that the decomposition is correct
//                 let decomp_i = decomp.iter().map(|d_j| d_j[i]).collect::<Vec<_>>();
//                 assert_eq!(v[i], recompose(&decomp_i, R::from(b)));
//             }
//         }
//     }
//
//     #[test]
//     fn test_decompose_balanced_polyring() {
//         let v = PolyR::from(get_test_vec());
//         for b in BASIS_TEST_RANGE {
//             let b_half = b / 2;
//             let decomp = decompose_balanced_polyring(&v, b, None);
//
//             for d_i in &decomp {
//                 for d_ij in d_i.coeffs() {
//                     let s_ij = SignedRepresentative::from(d_ij).0;
//                     assert!(s_ij.abs() as u128 <= b_half);
//                 }
//             }
//
//             assert_eq!(v, recompose(&decomp, R::from(b)));
//         }
//     }
//
//     #[test]
//     fn test_decompose_balanced_vec_polyring() {
//         let v = Vector::<PolyR>::from_fn(VEC_LENGTH, |i, _| {
//             PolyR::from(
//                 (0..(N as u64) * Q)
//                     .step_by((Q / (N as u64)) as usize)
//                     .map(|x| R::from(x + i as u64))
//                     .collect::<Vec<_>>(),
//             )
//         });
//         for b in BASIS_TEST_RANGE {
//             let b_half = R::from(b / 2);
//             let decomp = decompose_balanced_vec_polyring::<PolyR>(&v, b, None);
//
//             for v_i in &decomp {
//                 for v_ij in v_i.as_slice() {
//                     for v_ijk in v_ij.coeffs() {
//                         assert!(v_ijk <= b_half || v_ijk >= -b_half);
//                     }
//                 }
//             }
//
//             let mut recomposed = Vector::<PolyR>::zeros(v.len());
//             for (i, v_i) in decomp.iter().enumerate() {
//                 recomposed += v_i * PolyR::from_scalar(ring::Ring::pow(&R::from(b), &[i as u64]));
//             }
//             assert_eq!(v, recomposed);
//         }
//     }
//
//     #[test]
//     pub fn test_decompose_balanced_matrix() {
//         const NROWS: usize = 101;
//         const NCOLS: usize = 42;
//
//         let rng = &mut test_rng();
//         let mat = Matrix::<R>::rand(NROWS, NCOLS, rng);
//
//         for b in BASIS_TEST_RANGE {
//             let b_half = b.div_floor(2);
//
//             let decomp_size = balanced_decomposition_max_length(b, (Q - 1) as u128);
//             let decomp = decompose_matrix(&mat, b, decomp_size);
//             assert_eq!(decomp.nrows(), mat.nrows());
//             assert_eq!(decomp.ncols(), mat.ncols() * decomp_size);
//             for d_ij in decomp.iter() {
//                 assert!(Into::<SignedRepresentative>::into(*d_ij).0.abs() <= b_half as i128);
//             }
//
//             let pows_b = (0..decomp_size)
//                 .map(|i| R::from(b.pow(i as u32)))
//                 .collect::<Vec<_>>();
//             let recomp = recompose_matrix(&decomp, &pows_b);
//
//             assert_eq!(mat, recomp);
//
//             // Check that recompose is equivalent to right-multiplication by the gadget matrix
//             let gadget_vector = Vector::<R>::from_slice(&pows_b);
//             let gadget_matrix = Matrix::<R>::identity(NCOLS, NCOLS).kronecker(&gadget_vector);
//             let mat_g = &decomp * &gadget_matrix;
//             assert_eq!(recomp, mat_g);
//
//             // Check that A * M == recompose(A * decompose(M))
//             let mat_l = Matrix::<R>::rand(64, NROWS, rng);
//             let lhs = &mat_l * &mat;
//             let rhs = recompose_matrix(&(&mat_l * &decomp), &pows_b);
//             assert_eq!(lhs, rhs);
//         }
//     }
//
//     #[test]
//     fn test_left_right_decompose_symmetric_matrix() {
//         const SIZE: usize = 13;
//         let rng = &mut test_rng();
//
//         for b in BASIS_TEST_RANGE {
//             let decomp_size = balanced_decomposition_max_length(b, (Q - 1) as u128);
//
//             let pows_b = (0..decomp_size)
//                 .map(|i| R::from(b.pow(i as u32)))
//                 .collect::<Vec<_>>();
//             let mat = SymmetricMatrix::<R>::rand(SIZE * decomp_size, rng);
//
//             let expected: SymmetricMatrix<R> = recompose_matrix(
//                 &recompose_matrix(&mat.clone().into(), &pows_b.as_slice()).transpose(),
//                 &pows_b.as_slice(),
//             )
//             .into();
//
//             let actual = recompose_left_right_symmetric_matrix(&mat, &pows_b.as_slice());
//             assert_eq!(expected, actual);
//         }
//     }
// }
