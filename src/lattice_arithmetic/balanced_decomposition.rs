use std::iter::Sum;

use ark_ff::Field;
use num_traits::{One, Pow, Zero};
use rounded_div::RoundedDiv;

use crate::lattice_arithmetic::matrix::Vector;
use crate::lattice_arithmetic::poly_ring::{ConvertibleField, PolyRing, SignedRepresentative};
use crate::lattice_arithmetic::ring::Ring;
use crate::lattice_arithmetic::traits::IntegerDiv;

/// Returns the decomposition of `v' in basis `b', where centered representatives are used, i.e.,
/// v = \sum_i b^i v_i, with ||v_i||_\infty <= b/2.
pub fn decompose_balanced<R: ConvertibleField>(v: &R, b: u128, padding_size: Option<usize>) -> Vec<R>
{
    assert!(!b.is_zero() && !b.is_one(), "cannot decompose in basis 0 or 1");
    // TODO: not sure if this really necessary, but having even b's allow for more efficient divisions/remainders
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

    let mut decomp_bal = decomp_bal_signed.into_iter().map(|x| Into::<R>::into(SignedRepresentative(x))).collect::<Vec<R>>();

    if let Some(padding_size) = padding_size {
        assert!(decomp_bal.len() <= padding_size, "decomp_bal.len() = {} > padding_size = {}", decomp_bal.len(), padding_size);
        decomp_bal.resize(padding_size, R::zero());
    }
    decomp_bal
}

pub fn decompose_balanced_polyring<R: PolyRing>(v: &R, b: u128, padding_size: Option<usize>) -> Vec<R>
{
    let mut decomp: Vec<Vec<R::BaseRing>> = v.coeffs().iter().map(|ring_elem| decompose_balanced(ring_elem, b, padding_size)).collect(); // len(v) x decomp_size
    let rows = v.coeffs().len();
    let cols = decomp.iter().map(|d_i| d_i.len()).max().unwrap();
    // Pad each row to the same length `cols'
    for row in &mut decomp {
        row.resize(cols, R::BaseRing::zero());
    }

    (0..cols).map(|col| {
        R::from(
            (0..rows)
                .map(|row| decomp[row][col])
                .collect::<Vec<R::BaseRing>>()
        )
    }).collect()
}


pub fn decompose_balanced_vec<R: PolyRing>(v: &Vector<R>, b: u128, padding_size: Option<usize>) -> Vec<Vector<R>>
{
    let mut decomp: Vec<Vec<R>> = v.as_slice().iter().map(|ring_elem| decompose_balanced_polyring(ring_elem, b, padding_size)).collect(); // len(v) x decomp_size
    let rows = v.len();
    let cols = decomp.iter().map(|d_i| d_i.len()).max().unwrap();
    // Pad each row to the same length `cols'
    for row in &mut decomp {
        row.resize(cols, R::zero());
    }

    (0..cols).map(|col| {
        Vector::<R>::from((0..rows)
            .map(|row| decomp[row][col])
            .collect::<Vec<_>>())
    }).collect()
}

pub fn recompose<A, B>(v: &Vec<A>, b: B) -> A
    where A: std::ops::Mul<B, Output=A> + Copy + Sum, B: Field
{
    v.iter().enumerate().map(|(i, v_i)| *v_i * b.pow([i as u64])).sum()
}

#[cfg(test)]
mod tests {
    use crate::lattice_arithmetic::ntt::ntt_modulus;
    use crate::lattice_arithmetic::poly_ring::PolyRing;
    use crate::lattice_arithmetic::pow2_cyclotomic_poly_ring_ntt::Pow2CyclotomicPolyRingNTT;
    use crate::lattice_arithmetic::ring::Fq;

    use super::*;

    const N: usize = 128;
    const Q: u64 = ntt_modulus::<N>(16);
    const VEC_LENGTH: usize = 32;
    const BASIS_TEST_RANGE: [u128; 5] = [2, 4, 8, 16, 32];

    type R = Fq<Q>;
    type PolyR = Pow2CyclotomicPolyRingNTT<Q, N>;

    #[test]
    fn test_decompose_balanced() {
        let vs: Vec<R> = (0..Q).map(|v| R::from(v)).collect();
        for b in BASIS_TEST_RANGE {
            let b_half = R::from(b / 2);
            for v in &vs {
                let decomp = decompose_balanced(v, b, None);
                // assert ||v_i||_\infty <= b/2
                for v_i in &decomp {
                    assert!(*v_i <= b_half || *v_i >= -b_half);
                }

                // assert v = \sum_i b^i v_i
                assert_eq!(*v, recompose(&decomp, R::from(b)));
            }
        }
    }

    #[test]
    fn test_decompose_balanced_polyring() {
        let v = PolyR::from((0..(N as u64) * Q).step_by((Q / (N as u64)) as usize).map(|x| R::from(x)).collect::<Vec<_>>());
        for b in BASIS_TEST_RANGE {
            let b_half = R::from(b / 2);
            let decomp = decompose_balanced_polyring(&v, b, None);

            for v_i in &decomp {
                for v_ij in v_i.coeffs() {
                    assert!(v_ij <= b_half || v_ij >= -b_half);
                }
            }

            assert_eq!(v, recompose(&decomp, R::from(b)));
        }
    }

    #[test]
    fn test_decompose_balanced_vec() {
        let v =
            Vector::<PolyR>::from_fn(VEC_LENGTH,
                                     |i, _|
                                         PolyR::from((0..(N as u64) * Q).step_by((Q / (N as u64)) as usize).map(|x| R::from(x + i as u64)).collect::<Vec<_>>()),
            );
        for b in BASIS_TEST_RANGE {
            let b_half = R::from(b / 2);
            let decomp = decompose_balanced_vec::<PolyR>(&v, b, None);

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
}