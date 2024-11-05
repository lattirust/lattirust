//!
//! This module introduces a function `flatten_to_coeffs` on vectors of cyclotomic ring elements
//! to cheaply cast them into vectors of corresponding base field coefficients and its "inverse" `promote_from_coeffs`.
//!
use crate::PolyRing;

use super::{CyclotomicConfig, CyclotomicPolyRingNTTGeneral};

pub trait Flatten: PolyRing {
    fn flatten_to_coeffs(vec: Vec<Self>) -> Vec<Self::BaseRing>;
    fn promote_from_coeffs(vec: Vec<Self::BaseRing>) -> Option<Vec<Self>>;
}

impl<C: CyclotomicConfig<N>, const N: usize, const D: usize> Flatten
    for CyclotomicPolyRingNTTGeneral<C, N, D>
{
    fn flatten_to_coeffs(vec: Vec<Self>) -> Vec<C::BaseCRTField> {
        let (ptr, len, cap) = vec.into_raw_parts();

        unsafe { Vec::from_raw_parts(ptr as *mut C::BaseCRTField, len * D, cap * D) }
    }

    fn promote_from_coeffs(
        mut vec: Vec<C::BaseCRTField>,
    ) -> Option<Vec<CyclotomicPolyRingNTTGeneral<C, N, D>>> {
        if vec.len() % D != 0 {
            return None;
        }

        if vec.capacity() % D != 0 {
            vec.shrink_to_fit();
        }

        let (ptr, len, cap) = vec.into_raw_parts();

        Some(unsafe {
            Vec::from_raw_parts(
                ptr as *mut CyclotomicPolyRingNTTGeneral<C, N, D>,
                len / D,
                cap / D,
            )
        })
    }
}

#[cfg(test)]
mod tests {
    use ark_ff::One;

    use crate::cyclotomic_ring::{
        flatten::Flatten,
        models::goldilocks::{Fq3, RqNTT},
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
}
