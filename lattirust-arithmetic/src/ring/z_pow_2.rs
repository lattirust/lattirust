use std::io::{Read, Write};
use std::iter::{Product, Sum};
use std::num::Wrapping;
use std::ops::{Add, AddAssign, Mul, MulAssign, Sub, SubAssign};

use ark_serialize::{
    CanonicalDeserialize, CanonicalSerialize, Compress, SerializationError, Valid, Validate,
};
use ark_std::rand::{Rng, RngCore};
use ark_std::UniformRand;
use derive_more::{
    Add, AddAssign, Display, From, Into, Mul, MulAssign, Neg, Product, Sub, SubAssign, Sum,
};
use num_traits::{One, Zero};
use zeroize::Zeroize;

use crate::ring::{ConvertibleRing, Ring, SignedRepresentative, UnsignedRepresentative};
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
pub struct Z2_64(Wrapping<i64>);

impl Z2_64 {
    const MODULUS: u128 = 1 << 64;
    const MODULUS_HALF: u128 = Self::MODULUS / 2;
}

impl Zero for Z2_64 {
    fn zero() -> Self {
        Self(Wrapping(0))
    }

    fn is_zero(&self) -> bool {
        self.0 == Wrapping(0)
    }
}
impl One for Z2_64 {
    fn one() -> Self {
        Self(Wrapping(1))
    }
}

macro_rules! from_primitive_type {
    ($($t:ty),*) => {
        $(
            impl From<$t> for Z2_64 {
                fn from(x: $t) -> Self {
                    Self(Wrapping(x as i64))
                }
            }
        )*
    };
}
from_primitive_type!(u8, u16, u32, u64, u128, bool);

impl Modulus for Z2_64 {
    fn modulus() -> u128 {
        1u128 << 64
    }
}

impl CanonicalSerialize for Z2_64 {
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
        core::mem::size_of::<i64>()
    }
}

impl Valid for Z2_64 {
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

impl CanonicalDeserialize for Z2_64 {
    fn deserialize_with_mode<R: Read>(
        mut reader: R,
        compress: Compress,
        validate: Validate,
    ) -> Result<Self, SerializationError> {
        let mut bytes = [0u8; core::mem::size_of::<i64>()];
        reader.read_exact(&mut bytes)?;
        Ok(Self(Wrapping(<i64>::from_le_bytes(bytes))))
    }
}

impl FromRandomBytes<Self> for Z2_64 {
    fn byte_size() -> usize {
        8
    }

    fn try_from_random_bytes(bytes: &[u8]) -> Option<Self> {
        Some(Self(Wrapping(i64::from_be_bytes(bytes.try_into().ok()?))))
    }
}

macro_rules! impl_binop_ref {
    ($op: ident, $OpTrait: ident) => {
        impl<'a> $OpTrait<&'a Self> for Z2_64 {
            type Output = Self;

            fn $op(self, rhs: &'a Self) -> Self::Output {
                Self(self.0.$op(&rhs.0))
            }
        }

        impl<'a> $OpTrait<&'a mut Self> for Z2_64 {
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
        impl<'a> $OpTrait<&'a Self> for Z2_64 {
            fn $op(&mut self, rhs: &'a Self) {
                self.0.$op(&rhs.0)
            }
        }

        impl<'a> $OpTrait<&'a mut Self> for Z2_64 {
            fn $op(&mut self, rhs: &'a mut Self) {
                self.0.$op(&rhs.0)
            }
        }
    };
}

impl_binop_assign_ref!(add_assign, AddAssign);
impl_binop_assign_ref!(sub_assign, SubAssign);
impl_binop_assign_ref!(mul_assign, MulAssign);

impl<'a> Sum<&'a Self> for Z2_64 {
    fn sum<I: Iterator<Item = &'a Self>>(iter: I) -> Self {
        iter.fold(Self::zero(), |acc, x| acc + x)
    }
}

impl<'a> Product<&'a Self> for Z2_64 {
    fn product<I: Iterator<Item = &'a Self>>(iter: I) -> Self {
        iter.fold(Self::one(), |acc, x| acc * x)
    }
}

impl UniformRand for Z2_64 {
    fn rand<R: Rng + ?Sized>(rng: &mut R) -> Self {
        let mut bytes = [0u8; core::mem::size_of::<i64>()];
        rng.fill_bytes(&mut bytes);
        Self(Wrapping(i64::from_le_bytes(bytes)))
    }
}

impl Ring for Z2_64 {
    const ZERO: Self = Self(Wrapping(0));
    const ONE: Self = Self(Wrapping(1));
}

/// Map `[0, MODULUS_HALF] -> [0, MODULUS_HALF]` and `(-MODULUS_HALF, 0) -> (MODULUS_HALF, MODULUS)`
impl From<SignedRepresentative> for Z2_64 {
    fn from(value: SignedRepresentative) -> Self {
        Self(Wrapping(value.0 as i64))
    }
}

/// Map `[0, MODULUS_HALF] -> [0, MODULUS_HALF]` and `(MODULUS_HALF, MODULUS) -> (-MODULUS_HALF, 0)`
impl Into<SignedRepresentative> for Z2_64 {
    fn into(self) -> SignedRepresentative {
        SignedRepresentative(self.0 .0 as i128)
    }
}

impl Into<UnsignedRepresentative> for Z2_64 {
    fn into(self) -> UnsignedRepresentative {
        unimplemented!()
    }
}

impl ConvertibleRing for Z2_64 {}
