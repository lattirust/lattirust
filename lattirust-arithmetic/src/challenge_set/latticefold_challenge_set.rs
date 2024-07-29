use ark_ff::Field;

use crate::ring::{ PolyRing, Pow2CyclotomicPolyRingNTT, Zq };

pub trait OverField: PolyRing {
    type F: Field;

    fn field_to_base_ring(f: &Self::F) -> Self::BaseRing;
}

impl<const Q: u64, const N: usize> OverField for Pow2CyclotomicPolyRingNTT<Q, N> {
    type F = Zq<Q>;

    fn field_to_base_ring(f: &Self::F) -> Self::BaseRing {
        *f
    }
}
pub trait LatticefoldChallengeSet<R: OverField + PolyRing> {
    fn small_challenge_coefficient_from_random_bytes(i: usize, bs: &[u8]) -> R::BaseRing;
    fn small_challenge_from_random_bytes(d: usize, bs: &[u8]) -> R {
        <R as From<Vec<R::BaseRing>>>::from(
            (0..d)
                .map(|i| Self::small_challenge_coefficient_from_random_bytes(i, &bs[i..]))
                .collect()
        )
    }

    #[inline]
    fn small_challenge(bs: &[u8]) -> R {
        Self::small_challenge_from_random_bytes(1, bs)
    }

    #[inline]
    fn big_challenge_from_field(f: &R::F) -> R::BaseRing {
        R::field_to_base_ring(f)
    }
}

pub struct BinarySmallSet<const Q: u64, const N: usize>;

impl<const Q: u64, const N: usize> LatticefoldChallengeSet<Pow2CyclotomicPolyRingNTT<Q, N>>
for BinarySmallSet<Q, N> {
    fn small_challenge_coefficient_from_random_bytes(_i: usize, bs: &[u8]) -> Zq<Q> {
        if bs[0] == 0 { <Zq<Q> as Field>::ZERO } else { <Zq<Q> as Field>::ONE }
    }

    fn small_challenge(bs: &[u8]) -> Pow2CyclotomicPolyRingNTT<Q, N> {
        Self::small_challenge_from_random_bytes(N, bs)
    }
}
