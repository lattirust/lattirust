use ark_ff::Field;

use crate::ring::{NttRing, PolyRing, Pow2CyclotomicPolyRingNTT};

pub trait OverField<F: Field>: PolyRing {
    fn field_to_base_ring(f: &F) -> Self::BaseRing;
}

impl<Base: NttRing<N> + Field, const N: usize> OverField<Base> for Pow2CyclotomicPolyRingNTT<Base, N> {
    fn field_to_base_ring(f: &Base) -> Self::BaseRing {
        *f
    }
}

pub trait LatticefoldChallengeSet<F: Field, R: OverField<F> + PolyRing> {
    fn small_challenge_coefficient_from_random_bytes(i: usize, bs: &[u8]) -> R::BaseRing;
    fn small_challenge_from_random_bytes(d: usize, bs: &[u8]) -> R {
        <R as From<Vec<R::BaseRing>>>::from(
            (0..d)
                .map(|i| Self::small_challenge_coefficient_from_random_bytes(i, &bs[i..]))
                .collect(),
        )
    }

    #[inline]
    fn small_challenge(bs: &[u8]) -> R {
        Self::small_challenge_from_random_bytes(1, bs)
    }

    #[inline]
    fn big_challenge_from_field(f: &F) -> R::BaseRing {
        R::field_to_base_ring(f)
    }
}

pub struct BinarySmallSet<Base: Field, const N: usize> {
    _marker: std::marker::PhantomData<Base>,
}

impl<Base: NttRing<N> + Field, const N: usize> LatticefoldChallengeSet<Base, Pow2CyclotomicPolyRingNTT<Base, N>>
    for BinarySmallSet<Base, N>
{
    fn small_challenge_coefficient_from_random_bytes(_i: usize, bs: &[u8]) -> Base {
        if bs[0] == 0 {
            <Base as Field>::ZERO
        } else {
            <Base as Field>::ONE
        }
    }

    fn small_challenge(bs: &[u8]) -> Pow2CyclotomicPolyRingNTT<Base, N> {
        Self::small_challenge_from_random_bytes(N, bs)
    }
}
