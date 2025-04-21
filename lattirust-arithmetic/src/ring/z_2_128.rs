use std::io::{Read, Write};
use std::iter::{Product, Sum};
use std::num::Wrapping;
use std::ops::{Add, AddAssign, Mul, MulAssign, Sub, SubAssign};

use ark_serialize::{
    CanonicalDeserialize, CanonicalSerialize, Compress, SerializationError, Valid, Validate,
};
use ark_std::rand::Rng;
use ark_std::UniformRand;
use derive_more::{
    Add, AddAssign, Display, From, Into, Mul, MulAssign, Neg, Product, Sub, SubAssign, Sum,
};
use i256::i256;
use num_bigint::BigUint;
use num_traits::{One, ToPrimitive, Zero};
use zeroize::Zeroize;

use crate::ring::representatives::WithSignedRepresentative;
use crate::ring::Ring;
use crate::traits::{FromRandomBytes, Modulus};

#[derive(
    Clone,
    Copy,
    Debug,
    Display,
    From,
    Into,
    PartialEq,
    Eq,
    Hash,
    Default,
    PartialOrd,
    Ord,
    Neg,
    Add,
    AddAssign,
    Sub,
    SubAssign,
    Mul,
    MulAssign,
    Sum,
    Product,
    Zeroize,
)]
#[mul(forward)]
#[mul_assign(forward)]
#[repr(transparent)]
pub struct Z2_128(Wrapping<i128>);

impl Z2_128 {
    // TODO: the following should work, but doesn't
    // const MODULUS_i256: i256 = <i256 as ConstOne>::ONE << 128;
}

impl Zero for Z2_128 {
    fn zero() -> Self {
        Self(Wrapping(0))
    }

    fn is_zero(&self) -> bool {
        self.0 == Wrapping(0)
    }
}
impl One for Z2_128 {
    fn one() -> Self {
        Self(Wrapping(1))
    }
}

/// Map `[0, MODULUS_HALF) <-> [0, MODULUS_HALF)` and `[-MODULUS_HALF, 0) <-> [MODULUS_HALF, MODULUS)`
macro_rules! from_primitive_type {
    ($($t:ty),*) => {
        $(
            impl From<$t> for Z2_128 {
                fn from(x: $t) -> Self {
                    let unsigned = x as u128;
                    let signed: i128 = if unsigned < (1u128<<127) {
                        unsigned as i128
                    } else {
                        let unsigned_256 = i256::from_u128(unsigned);
                        let modulus = i256::from_u8(1u8) << 128;
                        let res: i256 = unsigned_256.sub(modulus);
                        res.as_i128()
                    };
                    Self(Wrapping(signed as i128))
                }
            }
        )*
    }
}

/// Map `[0, MODULUS_HALF) <-> [0, MODULUS_HALF)` and `[-MODULUS_HALF, 0) <-> [MODULUS_HALF, MODULUS)`
macro_rules! into_primitive_type {
    ($($t:ty),*) => {
        $(
            impl From<Z2_128> for $t {
                fn from(x: Z2_128) -> Self {
                    let signed = x.0 .0;
                    let unsigned: u128 = if signed >= 0 {
                        signed as u128
                    } else {
                        let signed_256 = i256::from_i128(signed);
                        let modulus = i256::from_u8(1u8) << 128;
                        let res: i256 = signed_256.add(modulus);
                        res.as_u128()
                    };
                    unsigned as $t
                }
            }
        )*
    };
}

from_primitive_type!(bool, u8, u16, u32, u64, u128);
into_primitive_type!(u8, u16, u32, u64, u128);

impl Modulus for Z2_128 {
    fn modulus() -> BigUint {
        BigUint::from(2u8).pow(128)
    }
}

impl CanonicalSerialize for Z2_128 {
    #[inline]
    fn serialize_with_mode<W: Write>(
        &self,
        mut writer: W,
        _compress: Compress,
    ) -> Result<(), SerializationError> {
        Ok(writer.write_all(&self.0 .0.to_le_bytes())?)
    }

    #[inline]
    fn serialized_size(&self, _compress: Compress) -> usize {
        core::mem::size_of::<i128>()
    }
}

impl Valid for Z2_128 {
    #[inline]
    fn check(&self) -> Result<(), SerializationError> {
        Ok(())
    }

    #[inline]
    fn batch_check<'a>(_batch: impl Iterator<Item = &'a Self>) -> Result<(), SerializationError>
    where
        Self: 'a,
    {
        Ok(())
    }
}

impl CanonicalDeserialize for Z2_128 {
    fn deserialize_with_mode<R: Read>(
        mut reader: R,
        _compress: Compress,
        _validate: Validate,
    ) -> Result<Self, SerializationError> {
        let mut bytes = [0u8; core::mem::size_of::<i128>()];
        reader.read_exact(&mut bytes)?;
        Ok(Self(Wrapping(<i128>::from_le_bytes(bytes))))
    }
}

impl FromRandomBytes<Self> for Z2_128 {
    fn has_no_bias() -> bool {
        true
    }
    fn needs_bytes() -> usize {
        8
    }

    fn try_from_random_bytes_inner(bytes: &[u8]) -> Option<Self> {
        Some(Self(Wrapping(i128::from_be_bytes(bytes.try_into().ok()?))))
    }
}

macro_rules! impl_binop_ref {
    ($op: ident, $OpTrait: ident) => {
        impl<'a> $OpTrait<&'a Self> for Z2_128 {
            type Output = Self;

            fn $op(self, rhs: &'a Self) -> Self::Output {
                Self(self.0.$op(&rhs.0))
            }
        }

        impl<'a> $OpTrait<&'a mut Self> for Z2_128 {
            type Output = Self;

            fn $op(self, rhs: &'a mut Self) -> Self::Output {
                Self(self.0.$op(&rhs.0))
            }
        }
    };
}

impl_binop_ref!(add, Add);
impl_binop_ref!(sub, Sub);
impl_binop_ref!(mul, Mul);

macro_rules! impl_binop_assign_ref {
    ($op: ident, $OpTrait: ident) => {
        impl<'a> $OpTrait<&'a Self> for Z2_128 {
            fn $op(&mut self, rhs: &'a Self) {
                self.0.$op(&rhs.0)
            }
        }

        impl<'a> $OpTrait<&'a mut Self> for Z2_128 {
            fn $op(&mut self, rhs: &'a mut Self) {
                self.0.$op(&rhs.0)
            }
        }
    };
}

impl_binop_assign_ref!(add_assign, AddAssign);
impl_binop_assign_ref!(sub_assign, SubAssign);
impl_binop_assign_ref!(mul_assign, MulAssign);

impl<'a> Sum<&'a Self> for Z2_128 {
    fn sum<I: Iterator<Item = &'a Self>>(iter: I) -> Self {
        iter.fold(Self::zero(), |acc, x| acc + x)
    }
}

impl<'a> Product<&'a Self> for Z2_128 {
    fn product<I: Iterator<Item = &'a Self>>(iter: I) -> Self {
        iter.fold(Self::one(), |acc, x| acc * x)
    }
}

impl UniformRand for Z2_128 {
    fn rand<R: Rng + ?Sized>(rng: &mut R) -> Self {
        let mut bytes = [0u8; core::mem::size_of::<i128>()];
        rng.fill_bytes(&mut bytes);
        Self(Wrapping(i128::from_le_bytes(bytes)))
    }
}

impl Ring for Z2_128 {
    const ZERO: Self = Self(Wrapping(0));
    const ONE: Self = Self(Wrapping(1));

    fn inverse(&self) -> Option<Self> {
        if self.0 .0 % 2 == 0 {
            None
        } else {
            let self_unsigned = BigUint::from(Into::<u128>::into(*self));
            let inv_unsigned_bigint = self_unsigned.modinv(&Self::modulus()).unwrap();
            let inv_unsigned = inv_unsigned_bigint.to_u128().unwrap();
            Some(Self::from(inv_unsigned))
        }
    }
}

impl From<i128> for Z2_128 {
    fn from(value: i128) -> Self {
        Self(Wrapping(value))
    }
}

impl From<Z2_128> for i128 {
    fn from(value: Z2_128) -> Self {
        value.0 .0
    }
}

impl WithSignedRepresentative for Z2_128 {
    type SignedRepresentative = i128;

    fn as_signed_representative(&self) -> Self::SignedRepresentative {
        self.0 .0
    }
}

#[cfg(test)]
mod test {
    use crate::*;

    use super::*;

    test_ring!(Z2_128, 100);
}
