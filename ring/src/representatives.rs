use ark_std::io::{Read, Write};
use ark_std::ops::Mul;

use ark_ff::{One, Zero};
use ark_serialize::Validate::Yes;
use ark_serialize::{
    CanonicalDeserialize, CanonicalSerialize, Compress, SerializationError, Valid, Validate,
};
use derive_more::{Add, AddAssign, Display, Mul, MulAssign, Product, Sub, SubAssign, Sum};
use serde::{Deserialize, Serialize};

// Work-around to allow us implementing From traits
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[repr(transparent)]
pub struct UnsignedRepresentative(pub u128);

#[derive(
    Clone,
    Copy,
    Debug,
    Display,
    PartialEq,
    Eq,
    Serialize,
    Deserialize,
    Add,
    AddAssign,
    Sub,
    SubAssign,
    Mul,
    MulAssign,
    Sum,
    Product,
)]
#[mul(forward)]
#[mul_assign(forward)]
#[repr(transparent)]
pub struct SignedRepresentative(pub i128);

impl Zero for SignedRepresentative {
    fn zero() -> Self {
        SignedRepresentative(0)
    }

    fn is_zero(&self) -> bool {
        self.0 == 0
    }
}

impl One for SignedRepresentative {
    fn one() -> Self {
        SignedRepresentative(1)
    }
}

impl<'a> Mul<&'a SignedRepresentative> for &'a SignedRepresentative {
    type Output = SignedRepresentative;

    fn mul(self, rhs: &'a SignedRepresentative) -> Self::Output {
        SignedRepresentative(self.0 * rhs.0)
    }
}

impl CanonicalSerialize for SignedRepresentative {
    fn serialize_with_mode<W: Write>(
        &self,
        mut writer: W,
        _compress: Compress,
    ) -> Result<(), SerializationError> {
        writer
            .write(&self.0.to_be_bytes())
            .map(|_| ())
            .map_err(SerializationError::IoError)
    }

    fn serialized_size(&self, _compress: Compress) -> usize {
        16
    }
}

impl Valid for SignedRepresentative {
    fn check(&self) -> Result<(), SerializationError> {
        Ok(())
    }
}

impl CanonicalDeserialize for SignedRepresentative {
    fn deserialize_with_mode<R: Read>(
        mut reader: R,
        _compress: Compress,
        validate: Validate,
    ) -> Result<Self, SerializationError> {
        let mut bytes = [0u8; 16];
        reader
            .read_exact(&mut bytes)
            .map_err(SerializationError::IoError)?;
        let value = i128::from_be_bytes(bytes);
        if validate == Yes {
            SignedRepresentative(value).check()?;
        }
        Ok(SignedRepresentative(value))
    }
}

#[cfg(test)]
mod tests {
    use crate::zn::z_q::Zq;
    use ark_ff::PrimeField;

    use super::*;

    const Q: u128 = 2u128.pow(61) - 1;
    const Q_HALF: i128 = (Q as i128 - 1) / 2;
    const TEST_STEP_SIZE: usize = (Q / 10) as usize;

    type F = Zq<{ Q as u64 }>;

    #[test]
    fn test_unsigned_representative() {
        for i in (0..Q).step_by(TEST_STEP_SIZE) {
            let f1 = F::from(i);
            let v1 = UnsignedRepresentative::from(f1);
            assert_eq!(i, v1.0);
        }
    }

    #[test]
    fn test_signed_representative() {
        assert_eq!(Q_HALF, F::MODULUS_MINUS_ONE_DIV_TWO.0[0] as i128);
        for i in (-Q_HALF..=Q_HALF).step_by(TEST_STEP_SIZE) {
            let v1 = SignedRepresentative(i);
            let f2 = F::from(v1);
            let v2 = SignedRepresentative::from(f2);
            assert_eq!(v1.0, v2.0);
        }
    }
}
