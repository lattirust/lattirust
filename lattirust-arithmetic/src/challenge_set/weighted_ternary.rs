use std::ops::Neg;

use num_traits::{One, Zero};

use crate::ring::pow2_cyclotomic_poly_ring::Pow2CyclotomicPolyRing;
use crate::ring::pow2_cyclotomic_poly_ring_ntt::Pow2CyclotomicPolyRingNTT;
use crate::ring::{NttRing, Ring};
use crate::traits::FromRandomBytes;

/// Challenge set $\\{-1, 0, 1\\} \subset \mathbb{Z}_q$, where $Pr\[C = 0\] = \frac{1}{2}$ and $Pr\[C = 1\] = Pr\[C = -1\] = \frac{1}{4}$
pub struct WeightedTernaryChallengeSet<R> {
    _marker: std::marker::PhantomData<R>,
}

impl<T: Zero + One + Neg<Output = T>> FromRandomBytes<T> for WeightedTernaryChallengeSet<T> {
    fn has_no_bias() -> bool {
        true
    }

    fn needs_bytes() -> usize {
        1
    }

    fn try_from_random_bytes_inner(bytes: &[u8]) -> Option<T> {
        let val = bytes[0] & 0b11; // Technically a u4 now
        if val == 0 || val == 3 {
            Some(T::zero())
        } else if val == 1 {
            Some(T::one())
        } else if val == 2 {
            Some(-T::one())
        } else {
            unreachable!()
        }
    }
}

pub struct WeightedTernaryPolyChallengeSet<R> {
    _marker: std::marker::PhantomData<R>,
}

impl<BaseRing: Ring, const N: usize> FromRandomBytes<Pow2CyclotomicPolyRing<BaseRing, N>>
    for WeightedTernaryPolyChallengeSet<Pow2CyclotomicPolyRing<BaseRing, N>>
where
    WeightedTernaryChallengeSet<BaseRing>: FromRandomBytes<BaseRing>,
{
    fn has_no_bias() -> bool {
        true
    }

    fn needs_bytes() -> usize {
        N * WeightedTernaryChallengeSet::<BaseRing>::byte_size()
    }

    fn try_from_random_bytes_inner(bytes: &[u8]) -> Option<Pow2CyclotomicPolyRing<BaseRing, N>> {
        assert_eq!(bytes.len(), Self::byte_size());
        let b = WeightedTernaryChallengeSet::<BaseRing>::byte_size();
        Some(Pow2CyclotomicPolyRing::<BaseRing, N>::from_fn(|i| {
            WeightedTernaryChallengeSet::<BaseRing>::try_from_random_bytes(
                &bytes[i * b..(i + 1) * b],
            )
            .unwrap()
        }))
    }
}

impl<BaseRing: NttRing<N>, const N: usize> FromRandomBytes<Pow2CyclotomicPolyRingNTT<BaseRing, N>>
    for WeightedTernaryPolyChallengeSet<Pow2CyclotomicPolyRingNTT<BaseRing, N>>
where
    WeightedTernaryChallengeSet<BaseRing>: FromRandomBytes<BaseRing>,
{
    fn has_no_bias() -> bool {
        true
    }

    fn needs_bytes() -> usize {
        Pow2CyclotomicPolyRing::<BaseRing, N>::byte_size()
    }

    fn try_from_random_bytes_inner(bytes: &[u8]) -> Option<Pow2CyclotomicPolyRingNTT<BaseRing, N>> {
        Pow2CyclotomicPolyRing::<BaseRing, N>::try_from_random_bytes(bytes).map(|x| x.into())
    }
}
