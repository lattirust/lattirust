use num_traits::{One, Zero};

use crate::ring::pow2_cyclotomic_poly_ring::Pow2CyclotomicPolyRing;
use crate::ring::pow2_cyclotomic_poly_ring_ntt::Pow2CyclotomicPolyRingNTT;
use crate::ring::Zq;
use crate::traits::FromRandomBytes;

/// Challenge set $\\{-1, 0, 1\\} \subset \mathbb{Z}_q$, where $Pr\[C = 0\] = \frac{1}{2}$ and $Pr\[C = 1\] = Pr\[C = -1\] = \frac{1}{4}$
pub struct WeightedTernaryChallengeSet<R> {
    _marker: std::marker::PhantomData<R>,
}

impl<const Q: u64> FromRandomBytes<Zq<Q>> for WeightedTernaryChallengeSet<Zq<Q>> {
    fn has_no_bias() -> bool {
        true
    }
    
    fn needs_bytes() -> usize {
        1
    }

    fn try_from_random_bytes_inner(bytes: &[u8]) -> Option<Zq<Q>> {
        let val = bytes[0] & 0b11; // Technically a u4 now
        if val == 0 || val == 3 {
            Some(Zq::<Q>::zero())
        } else if val == 1 {
            Some(Zq::<Q>::one())
        } else if val == 2 {
            Some(-Zq::<Q>::one())
        } else {
            unreachable!()
        }
    }
}

impl<const Q: u64, const N: usize> FromRandomBytes<Pow2CyclotomicPolyRing<Zq<Q>, N>>
    for WeightedTernaryChallengeSet<Pow2CyclotomicPolyRing<Zq<Q>, N>>
{
    fn has_no_bias() -> bool {
        true
    }
    
    fn needs_bytes() -> usize {
        N * WeightedTernaryChallengeSet::<Zq<Q>>::byte_size()
    }

    fn try_from_random_bytes_inner(bytes: &[u8]) -> Option<Pow2CyclotomicPolyRing<Zq<Q>, N>> {
        assert_eq!(bytes.len(), Self::byte_size());
        let b = WeightedTernaryChallengeSet::<Zq<Q>>::byte_size();
        Some(Pow2CyclotomicPolyRing::<Zq<Q>, N>::from_fn(|i| {
            WeightedTernaryChallengeSet::<Zq<Q>>::try_from_random_bytes(&bytes[i * b..(i + 1) * b])
                .unwrap()
        }))
    }
}

impl<const Q: u64, const N: usize> FromRandomBytes<Pow2CyclotomicPolyRingNTT<Q, N>>
    for WeightedTernaryChallengeSet<Pow2CyclotomicPolyRingNTT<Q, N>>
{
    fn has_no_bias() -> bool {
        true
    }
    
    fn needs_bytes() -> usize {
        Pow2CyclotomicPolyRing::<Zq<Q>, N>::byte_size()
    }

    fn try_from_random_bytes_inner(bytes: &[u8]) -> Option<Pow2CyclotomicPolyRingNTT<Q, N>> {
        Pow2CyclotomicPolyRing::<Zq<Q>, N>::try_from_random_bytes(bytes).map(|x| x.into())
    }
}
