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

use super::{pad_and_transpose, DecompositionFriendlySignedRepresentative};

//include lower bound, exclude upper bound
pub fn balanced_decomposition_max_length(b: u128, max: BigUint) -> usize {
    if max.is_zero() {
        0
    } else {
        max.to_f64().unwrap().log(b as f64).floor() as usize + 1 
    }
}

// standard base-b expansion (repeated Euclidean division) with digit balanced around zero + OpenFHE twist of dropping first digit
pub fn approximate_decompose_balanced<R>(v: &R, basis: u128, padding_size: Option<usize>) -> Vec<R>
where
    R: Ring + WithSignedRepresentative,
    R::SignedRepresentative: DecompositionFriendlySignedRepresentative,
{
    assert!(
        !basis.is_zero() && !basis.is_one() && !(basis == 2),
        "cannot decompose in basis 0,1 or 2"
    );
   assert!(basis % 2 == 0, "decomposition basis must be even");

    let basis_half_floor = basis.div_euclid(2);
    let b  = R::SignedRepresentative::try_from(basis as i128).unwrap();
    let b_half_floor = R::SignedRepresentative::try_from(basis_half_floor as i128).unwrap();

    let mut curr = Into::<R::SignedRepresentative>::into(*v);

    // approximate gadget decomposition is used; the first digit is ignored
    let mut rem = curr.clone() % b.clone();
    if rem.is_negative() { rem = rem + b.clone(); }
    let r0 = if rem > b_half_floor.clone() { rem - b.clone() } else { rem };
    curr = (curr - r0) / b.clone(); 

    let mut decomp_bal_signed = Vec::<R>::new();
    while !curr.is_zero() {
        let mut rem = curr.clone() % b.clone();
        if rem.is_negative() { rem = rem + b.clone(); }
        let digit = if rem > b_half_floor.clone() { rem - b.clone() } else { rem };

        decomp_bal_signed.push(R::try_from(digit.clone()).unwrap());
        curr = (curr - digit) / b.clone();
    }

    pad_zeros(&mut decomp_bal_signed, padding_size);
    decomp_bal_signed
}

pub fn approximate_decompose_balanced_vec<R>(v: &[R], b: u128, padding_size: Option<usize>) -> Vec<Vec<R>>
where
    R: Ring + WithSignedRepresentative,
    R::SignedRepresentative: DecompositionFriendlySignedRepresentative,
{
    let decomp: Vec<Vec<R>> = v
        .par_iter()
        .map(|v_i| approximate_decompose_balanced(v_i, b, padding_size))
        .collect(); // v.len() x decomp_size
    pad_and_transpose(decomp, padding_size) // decomp_size x v.len()
}

pub fn approximate_decompose_balanced_polyring<PR: PolyRing>(
    v: &PR,
    b: u128,
    padding_size: Option<usize>,
) -> Vec<PR>
where
    PR::BaseRing: WithSignedRepresentative,
    <PR::BaseRing as WithSignedRepresentative>::SignedRepresentative:
        DecompositionFriendlySignedRepresentative,
{
    approximate_decompose_balanced_vec::<PR::BaseRing>(v.coefficients().as_slice(), b, padding_size)
        .into_par_iter()
        .map(|v_i| PR::from(v_i))
        .collect()
}

    // Recompose sum_{i=0}^{vec.len()-1} vec[i] * basis^(i + start_pow)
    pub fn approximate_recompose<A, B>(vec: &[A], basis: B, start_pow: u64) -> A
    where
        A: Mul<B, Output = A> + Copy + Sum + Send + Sync,
        B: Ring,
    {
        // Precompute powers: basis^{start_pow}, basis^{start_pow+1}, ...
        let mut pow = B::one();
        for _ in 0..start_pow {
            pow = pow * basis;
        }
        let mut powers: Vec<B> = Vec::with_capacity(vec.len());
        if vec.len() > 0 {
            powers.push(pow);
            for _ in 1..vec.len() {
                pow = pow * basis;
                powers.push(pow);
            }
        }

        vec.par_iter()
            .zip(powers.par_iter())
            .map(|(d, p)| (*d) * (*p))
            .sum()
    }

    pub fn approximate_decompose_balanced_vec_polyring<PR: PolyRing>(
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
        .map(|ring_elem| approximate_decompose_balanced_polyring(ring_elem, b, padding_size))
        .collect(); // v.len() x decomp_size
    pad_and_transpose(decomp, padding_size)
        .into_par_iter()
        .map(Vector::from)
        .collect() // decomp_size x v.len()
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
    const BASIS_TEST_RANGE: [u128; 4] = [4, 8, 16, 32];

    type R = Zq1<Q>;
    type PolyR = Pow2CyclotomicPolyRingNTT<R, N>;

    fn first_digit_r0<R>(v: &R, b: u128) -> R
    where
        R: Ring + WithSignedRepresentative,
        R::SignedRepresentative: DecompositionFriendlySignedRepresentative,
    {
        let bs  = R::SignedRepresentative::try_from(b as i128).unwrap();
        let bh  = R::SignedRepresentative::try_from((b/2) as i128).unwrap();
        let mut rem = Into::<R::SignedRepresentative>::into(*v) % bs.clone();
        if rem.is_negative() { rem = rem + bs.clone(); }
        let r0s = if rem > bh { rem - bs } else { rem };
        R::try_from(r0s).unwrap()
    }


    #[test]
    fn test_decompose_balanced() {
        let vs: Vec<R> = (0..Q).map(|v| R::try_from(v).unwrap()).collect();
        for b in BASIS_TEST_RANGE {
            for v in &vs {
                let decomp = approximate_decompose_balanced(v, b, None);

                // Check that the decomposition length is correct
                let max = v.linf_norm();
                let max_decomp_length = balanced_decomposition_max_length(b, max);
                assert!(decomp.len() <= max_decomp_length);

                // Check that all entries are smaller than b/2
                for v_i in &decomp {
                    assert!(v_i.linf_norm() <= BigUint::from(b.div_floor(2)));
                }
                // Check that the decomposition is 
                let r0 = first_digit_r0(v, b);
                assert_eq!(*v, r0 + approximate_recompose(&decomp, R::try_from(b).unwrap(), 1));
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
            let decomp = approximate_decompose_balanced_vec(&v, b, None);

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
                let r0_i = first_digit_r0(&v[i], b);
                let decomp_i = decomp.iter().map(|d_j| d_j[i]).collect::<Vec<_>>();
                assert_eq!(v[i],r0_i + approximate_recompose(&decomp_i, R::try_from(b).unwrap(), 1));
            }
        }
    }

    #[test]
    fn test_decompose_balanced_polyring() {
        let v = PolyR::from_coefficient_array(get_test_vec());
        for b in BASIS_TEST_RANGE {
            let b_half = b / 2;
            let decomp: Vec<PolyR> = approximate_decompose_balanced_polyring(&v, b, None);

            for d_i in &decomp {
                for d_ij in d_i.coefficients() {
                    assert!(d_ij.linf_norm() <= b_half.into());
                }
            }
            let r0_coeffs: Vec<R> = v.coefficients().iter().map(|c| first_digit_r0(c, b)).collect();
            let r0_poly = PolyR::try_from_coefficients(r0_coeffs.as_slice()).unwrap();

            assert_eq!(v,r0_poly + approximate_recompose(&decomp, R::try_from(b).unwrap(), 1));
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
            )
            .unwrap()
        });
        for b in BASIS_TEST_RANGE {
            let b_half = b / 2;
            let decomp = approximate_decompose_balanced_vec_polyring::<PolyR>(&v.as_slice(), b, None);

            for v_i in &decomp {
                for v_ij in v_i.as_slice() {
                    for v_ijk in v_ij.coefficients() {
                        assert!(v_ijk.linf_norm() <= b_half.into());
                    }
                }
            }

            let r0_vec = Vector::<PolyR>::from_fn(VEC_LENGTH, |i, _| {
                let coeffs: Vec<R> = v[i].coefficients().iter().map(|c| first_digit_r0(c, b)).collect();
                PolyR::try_from_coefficients(coeffs.as_slice()).unwrap()
            });

            let mut recomposed = r0_vec.clone();
            for (i, v_i) in decomp.iter().enumerate() {
                recomposed +=
                    v_i * PolyR::from_scalar(Ring::pow(&R::try_from(b).unwrap(), i as u64 + 1));
            }
            assert_eq!(v, recomposed);
        }
    }

}
