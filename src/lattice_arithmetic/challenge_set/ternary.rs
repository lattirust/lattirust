use crate::lattice_arithmetic::poly_ring::{ConvertibleField, SignedRepresentative};
use crate::lattice_arithmetic::traits::FromRandomBytes;

pub struct TernaryChallengeSet<R> {
    _marker: std::marker::PhantomData<R>,
}

const SEC_PARAM: usize = 128;

impl<F: ConvertibleField> FromRandomBytes<F> for TernaryChallengeSet<F> {
    fn byte_size() -> usize {
        SEC_PARAM.next_power_of_two().ilog2() as usize + 2
    }

    fn try_from_random_bytes(bytes: &[u8]) -> Option<F> {
        assert_eq!(bytes.len(), Self::byte_size());
        let v = (bytes.last()? % 3) as i128;
        let s = SignedRepresentative(v - 1);
        Some(Into::<F>::into(s))
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Trit {
    MinusOne,
    Zero,
    One,
}

impl FromRandomBytes<Trit> for TernaryChallengeSet<Trit> {
    fn byte_size() -> usize {
        SEC_PARAM.next_power_of_two().ilog2() as usize + 2
    }

    fn try_from_random_bytes(bytes: &[u8]) -> Option<Trit> {
        assert_eq!(bytes.len(), Self::byte_size());
        let v = (bytes.last()? % 3) as i8;
        let res = match v - 1 {
            -1 => Trit::MinusOne,
            0 => Trit::Zero,
            1 => Trit::One,
            _ => unreachable!(),
        };
        Some(res)
    }
}