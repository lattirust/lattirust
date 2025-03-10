use std::cmp::PartialEq;
use std::io::{Read, Write};
use std::iter;
use std::iter::{Product, Sum};
use std::ops::{
    Add, AddAssign, BitAnd, BitAndAssign, BitXor, BitXorAssign, Div, DivAssign, Mul, MulAssign,
    Neg, Sub, SubAssign,
};
use std::str::FromStr;

use ark_ff::{
    AdditiveGroup, BigInt, FftField, Field, LegendreSymbol, PrimeField, SqrtPrecomputation,
};
use ark_serialize::{
    CanonicalDeserialize, CanonicalDeserializeWithFlags, CanonicalSerialize,
    CanonicalSerializeWithFlags, Compress, EmptyFlags, Flags, SerializationError, Valid, Validate,
};
use ark_std::rand::Rng;
use ark_std::UniformRand;
use derive_more::{Display, From, Into};
use num_bigint::BigUint;
use num_traits::{One, Zero};
use zeroize::Zeroize;

use crate::nimue::serialization::{FromBytes, ToBytes};
use crate::ring::Ring;
use crate::traits::{FromRandomBytes, Modulus, WithL2Norm, WithLinfNorm};

#[derive(
    Clone, Copy, Debug, Display, From, Into, PartialEq, Eq, Hash, Default, PartialOrd, Ord, Zeroize,
)]
#[repr(transparent)]
pub struct Z2(pub(crate) bool);

impl Z2 {
    const ZERO: Self = Self(false);
    const ONE: Self = Self(true);
}

impl Zero for Z2 {
    fn zero() -> Self {
        Self(false)
    }

    fn is_zero(&self) -> bool {
        self.0 == false
    }
}

impl One for Z2 {
    fn one() -> Self {
        Self(true)
    }
}

impl Neg for Z2 {
    type Output = Self;

    fn neg(self) -> Self::Output {
        self
    }
}

macro_rules! impl_binop_ref {
    ($bitop: ident, $op: ident, $OpTrait: ident) => {
        impl $OpTrait<Self> for Z2 {
            type Output = Self;

            fn $op(self, rhs: Self) -> Self::Output {
                Self(self.0.$bitop(rhs.0))
            }
        }

        impl<'a> $OpTrait<&'a Self> for Z2 {
            type Output = Self;

            fn $op(self, rhs: &'a Self) -> Self::Output {
                Self(self.0.$bitop(&rhs.0))
            }
        }

        impl<'a> $OpTrait<&'a mut Self> for Z2 {
            type Output = Self;

            fn $op(self, rhs: &'a mut Self) -> Self::Output {
                Self(self.0.$bitop(&rhs.0))
            }
        }
    };
}

impl_binop_ref!(bitxor, add, Add);
impl_binop_ref!(bitxor, sub, Sub);
impl_binop_ref!(bitand, mul, Mul);

macro_rules! impl_binop_assign_ref {
    ($bitop: ident, $op: ident, $OpTrait: ident) => {
        impl $OpTrait<Self> for Z2 {
            fn $op(&mut self, rhs: Self) {
                self.0.$bitop(rhs.0)
            }
        }

        impl<'a> $OpTrait<&'a Self> for Z2 {
            fn $op(&mut self, rhs: &'a Self) {
                self.0.$bitop(&rhs.0)
            }
        }

        impl<'a> $OpTrait<&'a mut Self> for Z2 {
            fn $op(&mut self, rhs: &'a mut Self) {
                self.0.$bitop(&rhs.0)
            }
        }
    };
}

impl_binop_assign_ref!(bitxor_assign, add_assign, AddAssign);
impl_binop_assign_ref!(bitxor_assign, sub_assign, SubAssign);
impl_binop_assign_ref!(bitand_assign, mul_assign, MulAssign);

impl Modulus for Z2 {
    fn modulus() -> BigUint {
        BigUint::from(2u64)
    }
}

macro_rules! from_primitive_type {
    ($($t:ty),*) => {
        $(
        // This is not a good implementation of From, since we're panicking on invalid values, but we need to implement From for u8, ..., u64 to implement Field :(
        impl From<$t> for Z2 {
            fn from(x: $t) -> Self {
                if x == <$t>::zero() {
                    <Self as Ring>::ZERO
                } else if x == <$t>::one() {
                    <Self as Ring>::ONE
                } else {
                    panic!("Invalid value for Z2: {}", x)
                }
            }
        }
        )*
    };
}
from_primitive_type!(
    u8,
    u16,
    u32,
    u64,
    u128,
    i8,
    i16,
    i32,
    i64,
    i128,
    BigUint,
    BigInt::<1>
);

impl CanonicalSerialize for Z2 {
    #[inline]
    fn serialize_with_mode<W: Write>(
        &self,
        mut writer: W,
        _compress: Compress,
    ) -> Result<(), SerializationError> {
        Ok(writer.write_all(&self.0.to_bytes().unwrap())?)
    }

    #[inline]
    fn serialized_size(&self, _compress: Compress) -> usize {
        core::mem::size_of::<bool>()
    }
}

impl CanonicalSerializeWithFlags for Z2 {
    fn serialize_with_flags<W: Write, F: Flags>(
        &self,
        writer: W,
        flags: F,
    ) -> Result<(), SerializationError> {
        if flags.u8_bitmask() != EmptyFlags.u8_bitmask() {
            return Err(SerializationError::UnexpectedFlags);
        }
        self.serialize_with_mode(writer, Compress::No)
    }

    fn serialized_size_with_flags<F: Flags>(&self) -> usize {
        self.serialized_size(Compress::No)
    }
}

impl Valid for Z2 {
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

impl CanonicalDeserializeWithFlags for Z2 {
    fn deserialize_with_flags<R: Read, F: Flags>(
        _reader: R,
    ) -> Result<(Self, F), SerializationError> {
        unimplemented!();
    }
}

impl CanonicalDeserialize for Z2 {
    fn deserialize_with_mode<R: Read>(
        mut reader: R,
        _compress: Compress,
        _validate: Validate,
    ) -> Result<Self, SerializationError> {
        let mut bytes = [0u8; core::mem::size_of::<bool>()];
        reader.read_exact(&mut bytes)?;
        Ok(Self(<bool>::from_bytes(&bytes)?))
    }
}

impl FromRandomBytes<Self> for Z2 {
    fn has_no_bias() -> bool {
        true
    }

    fn needs_bytes() -> usize {
        1
    }

    fn try_from_random_bytes_inner(bytes: &[u8]) -> Option<Self> {
        Some(Self((bytes.last().unwrap() & 1) != 0))
    }
}

impl Sum<Self> for Z2 {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::zero(), |acc, x| acc + x)
    }
}

impl<'a> Sum<&'a Self> for Z2 {
    fn sum<I: Iterator<Item = &'a Self>>(iter: I) -> Self {
        iter.fold(Self::zero(), |acc, x| acc + x)
    }
}

impl Product<Self> for Z2 {
    fn product<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::one(), |acc, x| acc * x)
    }
}

impl<'a> Product<&'a Self> for Z2 {
    fn product<I: Iterator<Item = &'a Self>>(iter: I) -> Self {
        iter.fold(Self::one(), |acc, x| acc * x)
    }
}

impl UniformRand for Z2 {
    fn rand<R: Rng + ?Sized>(rng: &mut R) -> Self {
        Self(bool::rand(rng))
    }
}

impl WithLinfNorm for Z2 {
    fn linf_norm(&self) -> BigUint {
        if self.0 {
            BigUint::one()
        } else {
            BigUint::zero()
        }
    }
}

impl WithL2Norm for Z2 {
    fn l2_norm_squared(&self) -> BigUint {
        if self.0 {
            BigUint::one()
        } else {
            BigUint::zero()
        }
    }
}

// impl Ring for Z2 {
//     const ZERO: Self = Self(false);
//     const ONE: Self = Self(true);
//
//     fn inverse(&self) -> Option<Self> {
//         if self.is_zero() {
//             None
//         } else {
//             Some(self.clone())
//         }
//     }
// }

impl<'a> DivAssign<&'a Self> for Z2 {
    fn div_assign(&mut self, other: &Self) {
        self.mul_assign(Field::inverse(other).unwrap());
    }
}

impl<'a> DivAssign<&'a mut Self> for Z2 {
    fn div_assign(&mut self, other: &'a mut Self) {
        self.div_assign(&*other)
    }
}

impl DivAssign<Self> for Z2 {
    fn div_assign(&mut self, other: Self) {
        self.div_assign(&other)
    }
}

impl<'a> Div<&'a Self> for Z2 {
    type Output = Self;

    fn div(mut self, other: &Self) -> Self {
        self.mul_assign(Field::inverse(other).unwrap());
        self
    }
}

impl<'a> Div<&'a mut Self> for Z2 {
    type Output = Self;

    fn div(mut self, other: &mut Self) -> Self {
        self.mul_assign(Field::inverse(other).unwrap());
        self
    }
}

impl Div<Self> for Z2 {
    type Output = Self;

    fn div(mut self, other: Self) -> Self {
        self.mul_assign(Field::inverse(&other).unwrap());
        self
    }
}

impl FromStr for Z2 {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "0" => Ok(Self::ZERO),
            "1" => Ok(Self::ONE),
            _ => Err(format!("Invalid value for Z2: {}", s)),
        }
    }
}

impl From<Z2> for BigUint {
    fn from(x: Z2) -> BigUint {
        if x.0 {
            BigUint::one()
        } else {
            BigUint::zero()
        }
    }
}

impl From<Z2> for BigInt<1> {
    fn from(x: Z2) -> BigInt<1> {
        BigInt::<1>::from(x.0 as u64)
    }
}

impl AdditiveGroup for Z2 {
    type Scalar = Self;
    const ZERO: Self = Self(false);
}

impl Field for Z2 {
    type BasePrimeField = Self;
    const SQRT_PRECOMP: Option<SqrtPrecomputation<Self>> = None;
    const ONE: Self = Self(true);

    fn extension_degree() -> u64 {
        1
    }

    fn to_base_prime_field_elements(&self) -> impl Iterator<Item = Self::BasePrimeField> {
        iter::once(*self)
    }

    fn from_base_prime_field_elems(
        elems: impl IntoIterator<Item = Self::BasePrimeField>,
    ) -> Option<Self> {
        let mut elems = elems.into_iter();
        let elem = elems.next()?;
        if elems.next().is_some() {
            return None;
        }
        Some(elem)
    }

    fn from_base_prime_field(elem: Self::BasePrimeField) -> Self {
        elem
    }

    fn from_random_bytes_with_flags<F: Flags>(bytes: &[u8]) -> Option<(Self, F)> {
        unimplemented!("{:?}", bytes);
    }

    fn legendre(&self) -> LegendreSymbol {
        unimplemented!("Legendre symbol is not defined for Z2")
    }

    fn square(&self) -> Self {
        *self
    }

    fn square_in_place(&mut self) -> &mut Self {
        self
    }

    fn inverse(&self) -> Option<Self> {
        if self.is_zero() {
            None
        } else {
            Some(*self)
        }
    }

    fn inverse_in_place(&mut self) -> Option<&mut Self> {
        match Field::inverse(self) {
            Some(inverse) => {
                *self = inverse;
                Some(self)
            }
            None => None,
        }
    }

    fn frobenius_map_in_place(&mut self, _power: usize) {}

    fn mul_by_base_prime_field(&self, elem: &Self::BasePrimeField) -> Self {
        *self * elem
    }
}

impl PrimeField for Z2 {
    type BigInt = BigInt<1>;
    const MODULUS: Self::BigInt = BigInt::<1>([1]);
    const MODULUS_MINUS_ONE_DIV_TWO: Self::BigInt = BigInt::zero();
    const MODULUS_BIT_SIZE: u32 = 1;
    const TRACE: Self::BigInt = BigInt::one();
    const TRACE_MINUS_ONE_DIV_TWO: Self::BigInt = BigInt::zero();

    fn from_bigint(repr: Self::BigInt) -> Option<Self> {
        match repr.0[0] {
            0 => Some(Self::ZERO),
            1 => Some(Self::ZERO),
            _ => None,
        }
    }

    fn into_bigint(self) -> Self::BigInt {
        self.into()
    }
}

impl FftField for Z2 {
    const GENERATOR: Self = Self::ONE;
    const TWO_ADICITY: u32 = 0;
    const TWO_ADIC_ROOT_OF_UNITY: Self = Self::ONE;
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::*;

    test_ring!(Z2, 100);
}
