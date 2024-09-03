use crate::ring::Ring;
use crate::traits::FromRandomBytes;

pub struct BinaryChallengeSet<R> {
    _marker: std::marker::PhantomData<R>,
}

impl<F: Ring> FromRandomBytes<F> for BinaryChallengeSet<F> {
    fn has_no_bias() -> bool {
        true
    }

    fn needs_bytes() -> usize {
        1
    }

    fn try_from_random_bytes_inner(bytes: &[u8]) -> Option<F> {
        let v = (bytes.last()? & 1) != 0;
        Some(F::from(v))
    }
}
