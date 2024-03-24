use ark_ff::Field;
use ark_std::iterable::Iterable;
use crypto_bigint::Zero;
use nalgebra::Scalar;
use crate::lattice_arithmetic::matrix::{Matrix, Vector};
use crate::lattice_arithmetic::poly_ring::{ConvertibleField, SignedRepresentative};
use crate::lattice_arithmetic::traits::FromRandomBytes;

pub struct TernaryChallengeSet<R> {
    _marker: std::marker::PhantomData<R>,
}

const SEC_PARAM: usize = 128;

impl<F: ConvertibleField> FromRandomBytes<F> for TernaryChallengeSet<F> {
    fn byte_size() -> usize {
        SEC_PARAM.next_power_of_two().ilog2() as usize + 2
    }

    fn try_from_random_bytes(bytes: &[u8]) -> Option<F> {
        assert_eq!(bytes.len(), Self::byte_size());
        let v = (bytes.last()? % 3) as i128;
        let s = SignedRepresentative(v - 1);
        Some(Into::<F>::into(s))
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Trit {
    MinusOne,
    Zero,
    One,
}

impl FromRandomBytes<Trit> for TernaryChallengeSet<Trit> {
    fn byte_size() -> usize {
        SEC_PARAM.next_power_of_two().ilog2() as usize + 2
    }

    fn try_from_random_bytes(bytes: &[u8]) -> Option<Trit> {
        assert_eq!(bytes.len(), Self::byte_size());
        let v = (bytes.last()? % 3) as i8;
        let res = match v - 1 {
            -1 => Trit::MinusOne,
            0 => Trit::Zero,
            1 => Trit::One,
            _ => unreachable!(),
        };
        Some(res)
    }
}

pub fn mul_f_trit<F: Field>(a: &Matrix<F>, b: &Matrix<Trit>) -> Matrix<F> {
    let mut c = Matrix::<F>::zeros(a.nrows(), b.ncols());
    for (i, a_i) in a.row_iter().enumerate() {
        for (j, b_j) in b.column_iter().enumerate() {
            for (a_ik, b_jk) in a_i.iter().zip(b_j.iter()) {
                match b_jk {
                    Trit::MinusOne => c[(i, j)] -= a_ik,
                    Trit::One => c[(i, j)] += a_ik,
                    Trit::Zero => {}
                }
            }
        }
    }
    c
}

pub fn mul_f_trit_sym<F: Field>(a: &Vec<Vec<F>>, b: &Matrix<Trit>) -> Matrix<F> {
    let mut c = Matrix::<F>::zeros(a.len(), b.ncols());
    for (i, a_i) in a.iter().enumerate() {
        for (j, b_j) in b.column_iter().enumerate() {
            for (k, b_jk) in b.iter().enumerate() {
                // Use the fact that a is symmetric
                let a_ik = if (k <= i+1) {a[i][k]} else {a[k][i]};
                match b_jk {
                    Trit::MinusOne => c[(i, j)] -= a_ik,
                    Trit::One => c[(i, j)] += a_ik,
                    Trit::Zero => {}
                }
            }
        }
    }
    c
}

// pub fn mul_F_Trit_vec<F: Field>(a: Matrix<F>, b: &Vector<Trit>) -> Matrix<F> {
//     let mut c = Vector::<F>::zeros(a.nrows());
//     for (i, a_i) in a.row_iter().enumerate() {
//         for (a_ik, b_jk) in a_i.iter().zip(b.iter()) {
//             match b_jk {
//                 Trit::MinusOne => c[i] -= a_ik,
//                 Trit::One => c[i] += a_ik,
//                 Trit::Zero => {}
//             }
//         }
//     }
//     c
// }

pub fn mul_trit_f<F: Field>(a: Matrix<Trit>, b: &Matrix<F>) -> Matrix<F> {
    let mut c = Matrix::<F>::zeros(a.nrows(), b.ncols());
    for (i, a_i) in a.row_iter().enumerate() {
        for (j, b_j) in b.column_iter().enumerate() {
            for (a_ik, b_jk) in a_i.iter().zip(b_j.iter()) {
                match a_ik {
                    Trit::MinusOne => c[(i, j)] -= b_jk,
                    Trit::One => c[(i, j)] += b_jk,
                    Trit::Zero => {}
                }
            }
        }
    }
    c
}

