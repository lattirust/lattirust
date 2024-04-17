use std::ops::{AddAssign, Neg, SubAssign};

use nalgebra::Scalar;
use num_traits::{One, Zero};

use crate::linear_algebra::Matrix;
use crate::ring::{ConvertibleRing, SignedRepresentative};
use crate::linear_algebra::SymmetricMatrix;
use crate::traits::FromRandomBytes;

pub struct TernaryChallengeSet<R> {
    _marker: std::marker::PhantomData<R>,
}

const SEC_PARAM: usize = 128;

impl<F: ConvertibleRing> FromRandomBytes<F> for TernaryChallengeSet<F> {
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

pub fn field_to_trit<F: Zero + One + Neg<Output = F> + PartialEq>(f: F) -> Option<Trit> {
    if f.is_zero() {
        Some(Trit::Zero)
    } else if f.is_one() {
        Some(Trit::One)
    } else if (-f).is_one() {
        Some(Trit::MinusOne)
    } else {
        None
    }
}

pub fn trit_to_field<F: Zero + One + Neg<Output = F>>(trit: Trit) -> F {
    match trit {
        Trit::MinusOne => -F::one(),
        Trit::Zero => F::zero(),
        Trit::One => F::one(),
    }
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

pub fn mul_f_trit<F: Scalar + Zero + SubAssign + AddAssign>(
    a: &Matrix<F>,
    b: &Matrix<Trit>,
) -> Matrix<F> {
    let mut c = Matrix::<F>::zeros(a.nrows(), b.ncols());
    for (i, a_i) in a.row_iter().enumerate() {
        for (j, b_j) in b.column_iter().enumerate() {
            for (a_ik, b_jk) in a_i.iter().zip(b_j.iter()) {
                match b_jk {
                    Trit::MinusOne => c[(i, j)] -= a_ik.clone(),
                    Trit::One => c[(i, j)] += a_ik.clone(),
                    Trit::Zero => {}
                }
            }
        }
    }
    c
}

/// Returns the symmetric matrix equal to c.transpose() * a * c
pub fn mul_trit_transpose_sym_trit<F: Clone + Zero + AddAssign + SubAssign>(
    a: &SymmetricMatrix<F>,
    c: &Matrix<Trit>,
) -> SymmetricMatrix<F> {
    assert_eq!(a.size(), c.nrows());
    let mut res = SymmetricMatrix::<F>::zero(c.ncols());
    for (l, c_l) in c.row_iter().enumerate() {
        for (k, c_k) in c.row_iter().enumerate() {
            for (i, c_li) in c_l
                .into_iter()
                .enumerate()
                .filter(|(_, c_li)| **c_li != Trit::Zero)
            {
                // We only need to set entries for j <= i (otherwise we'll have update off-diagonal entries twice), hence the .take(i+1)
                for (j, c_kj) in c_k
                    .into_iter()
                    .enumerate()
                    .take(i + 1)
                    .filter(|(_, c_kj)| **c_kj != Trit::Zero)
                {
                    let positive = c_li == c_kj;
                    if positive {
                        res.at_mut(i, j).add_assign(a.at(l, k).clone());
                    } else {
                        res.at_mut(i, j).sub_assign(a.at(l, k).clone());
                    }
                }
            }
        }
    }
    res
}

#[cfg(test)]
mod tests {
    use ark_std::test_rng;

    use crate::ring::Zq;

    use super::*;

    const Q: u64 = 65537;

    type F = Zq<Q>;

    const M: usize = 64;
    const N: usize = 128;
    const K: usize = 256;

    #[test]
    fn test_mul_f_trit() {
        let rng = &mut test_rng();
        let mat: Matrix<F> = Matrix::<F>::rand(M, N, rng);
        let trits_field: Matrix<F> = Matrix::<F>::rand_ternary(N, K, rng);
        let trits = trits_field.map(|f| field_to_trit(f).unwrap());

        let mat_trits = mul_f_trit(&mat, &trits);
        let mat_trits_field = mat * trits_field;

        println!("is:     {:?}", mat_trits);
        println!("should: {:?}", mat_trits_field);
        assert_eq!(mat_trits_field, mat_trits);
    }

    #[test]
    fn test_mul_trit_transpose_sym_trit() {
        let rng = &mut test_rng();
        let mat = SymmetricMatrix::<F>::rand(N, rng);
        let trits_field: Matrix<F> = Matrix::<F>::rand_ternary(N, M, rng);
        let trits = trits_field.map(|f| field_to_trit(f).unwrap());

        let mat_trits = mul_trit_transpose_sym_trit(&mat, &trits);
        let mat_trits_field: SymmetricMatrix<F> =
            (trits_field.transpose() * Into::<Matrix<F>>::into(mat) * trits_field).into();

        assert_eq!(mat_trits_field, mat_trits);
    }
}
