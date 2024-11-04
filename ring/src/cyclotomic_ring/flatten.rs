//!
//! This module introduces a function `flatten_to_coeffs` on vectors of cyclotomic ring elements
//! to cheaply cast them into vectors of corresponding base field coefficients.
//!

use std::ops::Deref;

use crate::PolyRing;

use super::{CyclotomicConfig, CyclotomicPolyRingNTTGeneral};

/// A trait to implement `flatten_to_coeffs` on the foreign type `Vec`.
pub trait Flatten<R: PolyRing>: Deref<Target = [R]> {
    fn flatten_to_coeffs(self) -> Vec<R::BaseRing>;
}

/// A trait to implement `promote_from_coeffs` on the foreign type `Vec`.
pub trait Promote<R: PolyRing>: Deref<Target = [R::BaseRing]> {
    fn promote_from_coeffs(self) -> Option<Vec<R>>;
}

impl<C: CyclotomicConfig<N>, const N: usize, const D: usize>
    Flatten<CyclotomicPolyRingNTTGeneral<C, N, D>> for Vec<CyclotomicPolyRingNTTGeneral<C, N, D>>
{
    fn flatten_to_coeffs(self) -> Vec<C::BaseCRTField> {
        let (ptr, len, cap) = self.into_raw_parts();

        unsafe { Vec::from_raw_parts(ptr as *mut C::BaseCRTField, len * D, cap * D) }
    }
}

impl<C: CyclotomicConfig<N>, const N: usize, const D: usize>
    Promote<CyclotomicPolyRingNTTGeneral<C, N, D>> for Vec<C::BaseCRTField>
{
    fn promote_from_coeffs(mut self) -> Option<Vec<CyclotomicPolyRingNTTGeneral<C, N, D>>> {
        if self.len() % D != 0 {
            return None;
        }

        if self.capacity() % D != 0 {
            self.shrink_to_fit();
        }

        let (ptr, len, cap) = self.into_raw_parts();

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
        flatten::*,
        models::goldilocks::{Fq3, RqNTT},
    };

    #[test]
    fn test_flatten_ntt() {
        let vec: Vec<RqNTT> = vec![RqNTT::one(), RqNTT::from(3u32), RqNTT::from(42u32)];
        let flattened = vec.flatten_to_coeffs();

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
        let promoted = vec.promote_from_coeffs().unwrap();

        assert_eq!(
            promoted,
            vec![RqNTT::one(), RqNTT::from(3u32), RqNTT::from(42u32)]
        )
    }
}
