use num_traits::{One, Zero};

use crate::lattice_arithmetic::challenge_set::challenge::ChallengeSet;
use crate::lattice_arithmetic::poly_ring::PolyRing;
use crate::lattice_arithmetic::pow2_cyclotomic_poly_ring::Pow2CyclotomicPolyRing;
use crate::lattice_arithmetic::pow2_cyclotomic_poly_ring_ntt::Pow2CyclotomicPolyRingNTT;
use crate::lattice_arithmetic::ring::{Ring, Zq};
use crate::lattice_arithmetic::traits::FromRandomBytes;

pub struct WeightedTernaryChallengeSet<R: Ring> {
    _marker: std::marker::PhantomData<R>,
}

impl<const Q: u64> FromRandomBytes<Zq<Q>> for WeightedTernaryChallengeSet<Zq<Q>> {
    fn byte_size() -> usize {
        1
    }

    fn from_random_bytes(bytes: &[u8]) -> Option<Zq<Q>> {
        assert_eq!(bytes.len(), 1);
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

/// Challenge set {-1, 0, 1} for Zq, where Pr[0] = 1/2, Pr[1] = Pr[-1] = 1/4
impl<const Q: u64> ChallengeSet<Zq<Q>> for WeightedTernaryChallengeSet<Zq<Q>> {}

impl<R: PolyRing> WeightedTernaryChallengeSet<R> {
    pub type Ring = R;
    pub type BaseRing = R::BaseRing;
    pub type BaseChallengeSet = WeightedTernaryChallengeSet<R::BaseRing>;
}

impl<const Q: u64, const N: usize> FromRandomBytes<Pow2CyclotomicPolyRing<Zq<Q>, N>> for WeightedTernaryChallengeSet<Pow2CyclotomicPolyRing<Zq<Q>, N>> {
    fn byte_size() -> usize {
        N * WeightedTernaryChallengeSet::<Zq<Q>>::byte_size()
    }

    fn from_random_bytes(bytes: &[u8]) -> Option<Self::Ring> {
        assert_eq!(bytes.len(), Self::byte_size());
        let b = WeightedTernaryChallengeSet::<Zq<Q>>::byte_size();
        Some(
            Pow2CyclotomicPolyRing::<Zq<Q>, N>::from_fn(|i|
                WeightedTernaryChallengeSet::<Zq<Q>>::from_random_bytes(&bytes[i * b..(i + 1) * b]).unwrap()
            )
        )
    }
}

impl<const Q: u64, const N: usize> ChallengeSet<Pow2CyclotomicPolyRing<Zq<Q>, N>> for WeightedTernaryChallengeSet<Pow2CyclotomicPolyRing<Zq<Q>, N>> {}


impl<const Q: u64, const N: usize> WeightedTernaryChallengeSet<Pow2CyclotomicPolyRingNTT<Q, N>> {}

impl<const Q: u64, const N: usize> FromRandomBytes<Pow2CyclotomicPolyRingNTT<Q, N>> for WeightedTernaryChallengeSet<Pow2CyclotomicPolyRingNTT<Q, N>> {
    fn byte_size() -> usize { Pow2CyclotomicPolyRing::<Zq<Q>, N>::byte_size() }

    fn from_random_bytes(bytes: &[u8]) -> Option<Self::Ring> { Pow2CyclotomicPolyRing::<Zq<Q>, N>::from_random_bytes(bytes).map(|x| x.into()) }
}

impl<const Q: u64, const N: usize> ChallengeSet<Pow2CyclotomicPolyRingNTT<Q, N>> for WeightedTernaryChallengeSet<Pow2CyclotomicPolyRingNTT<Q, N>> {}
