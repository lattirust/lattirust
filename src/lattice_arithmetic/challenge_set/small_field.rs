use ark_ff::Field;
use crate::lattice_arithmetic::traits::FromRandomBytes;

pub struct SmallFieldChallengeSet<F: Field> {
    _marker: std::marker::PhantomData<F>,
}

impl<F: Field> FromRandomBytes<F> for SmallFieldChallengeSet<F> {
    fn byte_size() -> usize {
        todo!()
    }

    fn try_from_random_bytes(bytes: &[u8]) -> Option<F> {
        todo!()
    }
}