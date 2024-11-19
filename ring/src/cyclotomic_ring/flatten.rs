//!
//! This module introduces a function `flatten_to_coeffs` on vectors of cyclotomic ring elements
//! to cheaply cast them into vectors of corresponding base field coefficients and its "inverse" `promote_from_coeffs`.
//!
use crate::PolyRing;
use ark_std::vec::*;

use super::{CyclotomicConfig, CyclotomicPolyRingGeneral, CyclotomicPolyRingNTTGeneral};

pub trait Flatten: PolyRing {
    fn flatten_to_coeffs(vec: Vec<Self>) -> Vec<Self::BaseRing> {
        let dimension = Self::dimension();

        let (ptr, len, cap) = vec.into_raw_parts();

        unsafe { Vec::from_raw_parts(ptr as *mut Self::BaseRing, len * dimension, cap * dimension) }
    }

    fn promote_from_coeffs(mut vec: Vec<Self::BaseRing>) -> Option<Vec<Self>> {
        let dimension = Self::dimension();

        if vec.len() % dimension != 0 {
            return None;
        }

        if vec.capacity() % dimension != 0 {
            vec.shrink_to_fit();
        }

        let (ptr, len, cap) = vec.into_raw_parts();

        Some(unsafe { Vec::from_raw_parts(ptr as *mut Self, len / dimension, cap / dimension) })
    }
}

impl<C: CyclotomicConfig<N>, const N: usize, const D: usize> Flatten
    for CyclotomicPolyRingNTTGeneral<C, N, D>
{
}

impl<C: CyclotomicConfig<N>, const N: usize, const D: usize> Flatten
    for CyclotomicPolyRingGeneral<C, N, D>
{
}

#[cfg(test)]
mod tests {
    use ark_ff::{One, UniformRand};
    use ark_std::vec::*;
    use rand::SeedableRng;
    use rand_chacha::ChaCha8Rng;

    use crate::cyclotomic_ring::{
        flatten::Flatten,
        models::goldilocks::{Fq3, RqNTT, RqPoly},
    };

    #[test]
    fn test_flatten_ntt() {
        let vec: Vec<RqNTT> = vec![RqNTT::one(), RqNTT::from(3u32), RqNTT::from(42u32)];
        let flattened = RqNTT::flatten_to_coeffs(vec);

        assert_eq!(
            flattened,
            vec![
                Fq3::one(),
                Fq3::one(),
                Fq3::one(),
                Fq3::one(),
                Fq3::one(),
                Fq3::one(),
                Fq3::one(),
                Fq3::one(),
                Fq3::from(3u32),
                Fq3::from(3u32),
                Fq3::from(3u32),
                Fq3::from(3u32),
                Fq3::from(3u32),
                Fq3::from(3u32),
                Fq3::from(3u32),
                Fq3::from(3u32),
                Fq3::from(42u32),
                Fq3::from(42u32),
                Fq3::from(42u32),
                Fq3::from(42u32),
                Fq3::from(42u32),
                Fq3::from(42u32),
                Fq3::from(42u32),
                Fq3::from(42u32),
            ]
        )
    }

    #[test]
    fn test_promote_ntt() {
        let vec: Vec<Fq3> = vec![
            Fq3::one(),
            Fq3::one(),
            Fq3::one(),
            Fq3::one(),
            Fq3::one(),
            Fq3::one(),
            Fq3::one(),
            Fq3::one(),
            Fq3::from(3u32),
            Fq3::from(3u32),
            Fq3::from(3u32),
            Fq3::from(3u32),
            Fq3::from(3u32),
            Fq3::from(3u32),
            Fq3::from(3u32),
            Fq3::from(3u32),
            Fq3::from(42u32),
            Fq3::from(42u32),
            Fq3::from(42u32),
            Fq3::from(42u32),
            Fq3::from(42u32),
            Fq3::from(42u32),
            Fq3::from(42u32),
            Fq3::from(42u32),
        ];
        let promoted = RqNTT::promote_from_coeffs(vec).unwrap();

        assert_eq!(
            promoted,
            vec![RqNTT::one(), RqNTT::from(3u32), RqNTT::from(42u32)]
        )
    }

    #[test]
    fn test_flatten_promote_coeff() {
        let mut rng = ChaCha8Rng::seed_from_u64(0);

        let orig: Vec<RqPoly> = (0..3).map(|_| RqPoly::rand(&mut rng)).collect();
        let flattened = RqPoly::flatten_to_coeffs(orig.clone());
        let promoted = RqPoly::promote_from_coeffs(flattened).unwrap();

        assert_eq!(orig, promoted);
    }
}
