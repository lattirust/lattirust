use std::io::{Read, Write};
use std::iter::{Product, Sum};
use std::num::Wrapping;
use std::ops::{Add, AddAssign, Mul, MulAssign, Sub, SubAssign};

use ark_serialize::{
    CanonicalDeserialize, CanonicalSerialize, Compress, SerializationError, Valid, Validate,
};
use ark_std::rand::Rng;
use ark_std::UniformRand;
use derive_more::{ Display, From, Into, Neg, Product, Sum,
};
use i256::i256;
use num_bigint::{BigInt, BigUint};
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
    Hash,
    Default,
    PartialOrd,
    Ord,
    Neg,
    Sum,
    Product,
    Zeroize,
)]
#[repr(transparent)]
pub struct Z2_k<const K: u32>(Wrapping<i128>);

impl<const K: u32> Z2_k<K> {
    //ensures values don't exceed the modulus
    const _ASSERT_K_VALID: () = {
        assert!(K > 0 && K <= 128, "K must be between 1 and 128");
    };
    const _VALIDATE: () = Self::_ASSERT_K_VALID;

    const MODULUS_HALF: u128 = 1u128 << (K - 1);

    // Reduce the value to fit within the range of the ring
    fn reduce(value: i128) -> i128 {
        if K >= 128 {
            value
        } else {
            let mask = (1u128 << K) - 1;
            let unsigned = value as u128 & mask;
            let modulus_half = Self::MODULUS_HALF;

            if unsigned < modulus_half {
                unsigned as i128
            } else {
                (unsigned as i128) - (1i128 << K)
            }
        }
    }
}

impl<const K: u32> Zero for Z2_k<K> {
    fn zero() -> Self {
        Self(Wrapping(0))
    }

    fn is_zero(&self) -> bool {
        self.0 == Wrapping(0)
    }
}

impl<const K: u32> One for Z2_k<K> {
    fn one() -> Self {
        Self(Wrapping(1))
    }
}

impl<const K: u32> Mul for Z2_k<K> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        let result = self.0 .0.wrapping_mul(rhs.0 .0);
        Self(Wrapping(Self::reduce(result)))
    }
}

impl<const K: u32> MulAssign for Z2_k<K> {
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}

impl<const K: u32> Add for Z2_k<K> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        let result = self.0.0.wrapping_add(rhs.0.0);
        Self(Wrapping(Self::reduce(result)))
    }
}

impl<const K: u32> Sub for Z2_k<K> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        let result = self.0.0.wrapping_sub(rhs.0.0);
        Self(Wrapping(Self::reduce(result)))
    }
}

impl<const K: u32> AddAssign for Z2_k<K> {
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}

impl<const K: u32> SubAssign for Z2_k<K> {
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs;
    }
}

impl<const K: u32> PartialEq for Z2_k<K> {
    fn eq(&self, other: &Self) -> bool {
        Self::reduce(self.0.0) == Self::reduce(other.0.0)
    }
}

impl<const K: u32> Eq for Z2_k<K> {}

/// Map `[0, MODULUS_HALF) <-> [0, MODULUS_HALF)` and `[-MODULUS_HALF, 0) <-> [MODULUS_HALF, MODULUS)`
macro_rules! from_primitive_type {
    ($($t:ty),*) => {
        $(
            impl<const K: u32> From<$t> for Z2_k<K> {
                fn from(x: $t) -> Self {
                    let unsigned = x as u128;
                    let modulus_half = Self::MODULUS_HALF;
                    let signed: i128 = if unsigned < modulus_half {
                        unsigned as i128
                    } else {
                        let unsigned_256 = i256::from_u128(unsigned);
                        let modulus = i256::from_u8(1u8) << K;
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
            impl<const K: u32> From<Z2_k<K>> for $t {
                fn from(x: Z2_k<K>) -> Self {
                    let signed = x.0 .0;
                    let unsigned: u128 = if signed >= 0 {
                        signed as u128
                    } else {
                        let signed_256 = i256::from_i128(signed);
                        let modulus = i256::from_u8(1u8) << K;
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

impl<const K: u32> Modulus for Z2_k<K> {
    fn modulus() -> BigUint {
        BigUint::from(2u8).pow(K)
    }
}

impl<const K: u32> CanonicalSerialize for Z2_k<K> {
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

impl<const K: u32> Valid for Z2_k<K> {
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

impl<const K: u32> CanonicalDeserialize for Z2_k<K> {
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

impl<const K: u32> FromRandomBytes<Self> for Z2_k<K> {
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
        impl<'a, const K: u32> $OpTrait<&'a Self> for Z2_k<K> {
            type Output = Self;

            fn $op(self, rhs: &'a Self) -> Self::Output {
                Self(self.0.$op(&rhs.0))
            }
        }

        impl<'a, const K: u32> $OpTrait<&'a mut Self> for Z2_k<K> {
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
        impl<'a, const K: u32> $OpTrait<&'a Self> for Z2_k<K> {
            fn $op(&mut self, rhs: &'a Self) {
                self.0.$op(&rhs.0)
            }
        }

        impl<'a, const K: u32> $OpTrait<&'a mut Self> for Z2_k<K> {
            fn $op(&mut self, rhs: &'a mut Self) {
                self.0.$op(&rhs.0)
            }
        }
    };
}

impl_binop_assign_ref!(add_assign, AddAssign);
impl_binop_assign_ref!(sub_assign, SubAssign);
impl_binop_assign_ref!(mul_assign, MulAssign);

impl<'a, const K: u32> Sum<&'a Self> for Z2_k<K> {
    fn sum<I: Iterator<Item = &'a Self>>(iter: I) -> Self {
        iter.fold(Self::zero(), |acc, x| acc + x)
    }
}

impl<'a, const K: u32> Product<&'a Self> for Z2_k<K> {
    fn product<I: Iterator<Item = &'a Self>>(iter: I) -> Self {
        iter.fold(Self::one(), |acc, x| acc * x)
    }
}

impl<const K: u32> UniformRand for Z2_k<K> {
    fn rand<R: Rng + ?Sized>(rng: &mut R) -> Self {
        let bytes_needed = ((K + 7) / 8) as usize;
        let bytes_needed = bytes_needed.min(16); // Cap at i128 size

        let mut bytes = [0u8; 16];
        rng.fill_bytes(&mut bytes[..bytes_needed]);
        let raw_value = i128::from_le_bytes(bytes);

        Self(Wrapping(Self::reduce(raw_value)))
    }
}

impl<const K: u32> Ring for Z2_k<K> {
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

impl<const K: u32> From<i128> for Z2_k<K> {
    fn from(value: i128) -> Self {
        Self(Wrapping(Self::reduce(value)))
    }
}

impl<const K: u32> From<Z2_k<K>> for i128 {
    fn from(value: Z2_k<K>) -> Self {
        value.0 .0
    }
}

impl<const K: u32> WithSignedRepresentative for Z2_k<K> {
    type SignedRepresentative = i128;

    fn as_signed_representative(&self) -> Self::SignedRepresentative {
        self.0 .0
    }

    fn signed_representative_to_bigint(repr: &Self::SignedRepresentative) -> BigInt {
        BigInt::from(*repr)
    }
    
    fn signed_representative_from_bigint(value: BigInt) -> Option<Self::SignedRepresentative> {
        value.to_i128()
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::*;

    // Test the generic version with different K values
    //test_ring!(Z2_k<8>, 100);
    //test_ring!(Z2_k<16>, 100);
    //test_ring!(Z2_k<32>, 100);
    test_ring!(Z2_k<37>, 100);
    //test_ring!(Z2_k<50>, 100);
    //test_ring!(Z2_k<64>, 100);
    //test_ring!(Z2_k<128>, 100);

    #[test]
    fn test_different_moduli() {
        // Test Z2_k<8> (modulus = 256)
        let a: Z2_k<8> = Z2_k::from(200u8);
        let b: Z2_k<8> = Z2_k::from(100u8);
        let c = a + b;

        // Test Z2_k<4> (modulus = 16)
        let x: Z2_k<4> = Z2_k::from(10u8);
        let y: Z2_k<4> = Z2_k::from(8u8);
        let z = x + y;

        // Test Z2_k<37> (modulus = 2^37)
        let a_37: Z2_k<37> = Z2_k::from(1000000000u32);
        let b_37: Z2_k<37> = Z2_k::from(2000000000u32);
        let c_37 = a_37 + b_37;

        println!("Z2_8: {} + {} = {}", a.0 .0, b.0 .0, c.0 .0);
        println!("Z2_4: {} + {} = {}", x.0 .0, y.0 .0, z.0 .0);
        println!("Z2_37: {} + {} = {}", a_37.0 .0, b_37.0 .0, c_37.0 .0);
    }
}
