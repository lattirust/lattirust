use std::io::Write;

use num_traits::Zero;

use crate::lattice_arithmetic::matrix::Vector;
use crate::lattice_arithmetic::poly_ring::PolyRing;
use crate::lattice_arithmetic::ring::Ring;
use crate::lattice_arithmetic::traits::IntegerDiv;

/// Returns the decomposition of `v' in basis `b', where centered representatives are used, i.e.,
/// v = \sum_i b^i v_i, with ||v_i||_\infty <= b/2.
pub fn decompose_balanced<R: Ring + IntegerDiv>(v: &R, b: R, padding_size: Option<usize>) -> Vec<R> {
    assert!(b > R::one(), "cannot decompose in basis 0 or 1");
    let two = R::from(2u128);

    let rem_2 = b - (two * b.integer_div(two));
    assert!(rem_2.is_one(), "decomposition basis must be odd");

    let b_half_floor = b.integer_div(R::from(2u128));
    let mut decomp_bal = Vec::<R>::new();
    let mut curr = v.clone();
    loop {
        let rem = curr - (b * curr.integer_div(b)); // rem = curr % b


        // Ensure digit is in [-b/2, b/2)
        if rem <= b_half_floor {
            decomp_bal.push(rem);
            curr = curr.integer_div(b);
        } else {
            decomp_bal.push(rem - b);
            let carry = rem.div_round(b);
            curr = curr.integer_div(b) + carry;
        }

        if curr.is_zero() {
            break;
        }
    }

    if let Some(padding_size) = padding_size {
        assert!(decomp_bal.len() <= padding_size, "decomp_bal.len() = {} > padding_size = {}", decomp_bal.len(), padding_size);
        decomp_bal.resize(padding_size, R::zero());
    }
    decomp_bal
}

pub fn decompose_balanced_polyring<R: PolyRing>(v: &R, b: R::BaseRing, padding_size: Option<usize>) -> Vec<R> {
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


pub fn decompose_balanced_vec<R: PolyRing>(v: &Vector<R>, b: R::BaseRing, padding_size: Option<usize>) -> Vec<Vector<R>> {
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

pub fn recompose<A: Ring, B: Ring>(v: Vec<A>, b: B) -> A
    where A: std::ops::Mul<B, Output=A> {
    v.iter().enumerate().map(|(i, v_i)| *v_i * b.pow(&[i as u64])).sum()
}

mod tests {
    use crate::lattice_arithmetic::balanced_decomposition::{decompose_balanced, decompose_balanced_polyring, decompose_balanced_vec, recompose};
    use crate::lattice_arithmetic::matrix::Vector;
    use crate::lattice_arithmetic::poly_ring::PolyRing;
    use crate::lattice_arithmetic::pow2_cyclotomic_poly_ring_ntt::Pow2CyclotomicPolyRingNTT;
    use crate::lattice_arithmetic::ring::{Ring, Zq};
    use crate::lattice_arithmetic::traits::{IntegerDiv, Modulus};

    const Q: u64 = 8191;
    const N: usize = 128;
    const VEC_LENGTH: usize = 32;

    type R = Zq<Q>;
    type PolyR = Pow2CyclotomicPolyRingNTT<R, N>;

    #[test]
    fn test_decompose_balanced() {
        let vs: Vec<R> = (0..Q).map(|v| R::from(v)).collect();
        let bs: Vec<R> = (2..Q / 100).filter(|x| *x % 2 == 1).map(|b| R::from(b)).collect();
        for b in bs {
            let b_half = b.integer_div(R::from(2u128));
            for v in &vs {
                let decomp = decompose_balanced(v, b, None);
                // assert ||v_i||_\infty <= b/2
                for v_i in &decomp {
                    assert!(*v_i <= b_half || *v_i >= -b_half);
                }

                // assert v = \sum_i b^i v_i
                assert_eq!(*v, recompose(decomp, b));
            }
        }
    }

    // fn pp(v: &PolyR) -> String {
    //     v.coeffs().iter().map(|v_i| format!("{:?}", u64::from(*v_i))).collect::<Vec<_>>().join(", ")
    // }

    #[test]
    fn test_decompose_balanced_polyring() {
        let v = PolyR::from((0..(N as u64) * Q).step_by((Q / (N as u64)) as usize).map(|x| R::from(x)).collect::<Vec<_>>());
        let bs: Vec<R> = (2..Q / 100).filter(|x| *x % 2 == 1).map(|b| R::from(b)).collect();
        for b in bs {
            let b_half = b.integer_div(R::from(2u128));
            let decomp = decompose_balanced_polyring(&v, b, None);

            for v_i in &decomp {
                for v_ij in v_i.coeffs() {
                    assert!(v_ij <= b_half || v_ij >= -b_half);
                }
            }

            assert_eq!(v, recompose(decomp, b));
        }
    }

    #[test]
    fn test_decompose_balanced_vec() {
        let v =
            Vector::<PolyR>::from_fn(VEC_LENGTH,
                                     |i, _|
                                         PolyR::from((0..(N as u64) * Q).step_by((Q / (N as u64)) as usize).map(|x| R::from(x + i as u64)).collect::<Vec<_>>()),
            );
        let bs: Vec<R> = (2..Q / 100).filter(|x| *x % 2 == 1).map(|b| R::from(b)).collect();
        for b in bs {
            let b_half = b.integer_div(R::from(2u128));
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
                recomposed += v_i * PolyR::from_value(b.pow(&[i as u64]));
            }
            assert_eq!(v, recomposed);
        }
    }
}