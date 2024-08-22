use crate::challenge_set::ternary::TernaryChallengeSet;
use crate::ring::{ConvertibleRing, SignedRepresentative};
use crate::traits::FromRandomBytes;

pub struct BinaryChallengeSet<R> {
    _marker: std::marker::PhantomData<R>,
}

impl<F: ConvertibleRing> FromRandomBytes<F> for BinaryChallengeSet<F> {
    fn byte_size() -> usize {
        1
    }

    fn try_from_random_bytes(bytes: &[u8]) -> Option<F> {
        assert_eq!(bytes.len(), Self::byte_size());
        let v = (bytes.last()? & 1) as i128;
        let s = SignedRepresentative(v);
        Some(Into::<F>::into(s))
    }
}