use std::fmt::{Debug, Display, Formatter};
use std::hash::Hash;
use std::io::{Read, Write};
use std::iter::Product;
use std::ops::{Mul, MulAssign, Neg};

use ark_serialize::{
    CanonicalDeserialize, CanonicalSerialize, Compress, SerializationError, Valid, Validate,
};
use ark_std::rand::Rng;
use ark_std::UniformRand;
use derive_more::{Add, AddAssign, From, Into, Sub, SubAssign, Sum};
use num_bigint::BigUint;
use num_traits::{One, Zero};

use crate::linear_algebra::SVector;
use crate::ring::{PolyRing, Ring};
use crate::traits::{
    FromRandomBytes, Modulus, WithConjugationAutomorphism, WithL2Norm, WithLinfNorm,
};

#[derive(
    Clone, Copy, Debug, Eq, PartialEq, Hash, Add, AddAssign, Sum, Sub, SubAssign, From, Into,
)]
pub struct Pow2CyclotomicPolyRing<BaseRing: Ring, const N: usize>(SVector<BaseRing, N>);

impl<BaseRing: Ring, const N: usize> Pow2CyclotomicPolyRing<BaseRing, N> {
    #[allow(dead_code)]
    pub(crate) type Inner = SVector<BaseRing, N>;
    pub fn from_fn<F>(f: F) -> Self
    where
        F: FnMut(usize) -> BaseRing,
    {
        let coeffs = core::array::from_fn(f);
        Self(Self::Inner::const_from_array(coeffs))
    }
    const fn const_from_element(elem: BaseRing) -> Self {
        let mut coeffs = [BaseRing::ZERO; N];
        coeffs[0] = elem;
        Self(Self::Inner::const_from_array(coeffs))
    }
}

impl<BaseRing: Ring, const N: usize> From<[BaseRing; N]> for Pow2CyclotomicPolyRing<BaseRing, N> {
    fn from(value: [BaseRing; N]) -> Self {
        Self(Self::Inner::const_from_array(value))
    }
}

impl<BaseRing: Ring, const N: usize> Modulus for Pow2CyclotomicPolyRing<BaseRing, N> {
    fn modulus() -> BigUint {
        BaseRing::modulus()
    }
}

macro_rules! impl_from_primitive_type {
    ($primitive_type: ty) => {
        impl<BaseRing: Ring, const N: usize> From<$primitive_type>
            for Pow2CyclotomicPolyRing<BaseRing, N>
        where
            BaseRing: From<$primitive_type>,
        {
            fn from(value: $primitive_type) -> Self {
                Self::from_scalar(BaseRing::from(value))
            }
        }
    };
}

macro_rules! impl_try_from_primitive_type {
    ($primitive_type: ty) => {
        impl<BaseRing: Ring, const N: usize> TryFrom<$primitive_type>
            for Pow2CyclotomicPolyRing<BaseRing, N>
        {
            type Error = <BaseRing as TryFrom<$primitive_type>>::Error;

            fn try_from(value: $primitive_type) -> Result<Self, Self::Error> {
                Ok(Self::from_scalar(BaseRing::try_from(value)?))
            }
        }
    };
}

impl_from_primitive_type!(BaseRing);
impl_from_primitive_type!(bool);
impl_try_from_primitive_type!(u8);
impl_try_from_primitive_type!(u16);
impl_try_from_primitive_type!(u32);
impl_try_from_primitive_type!(u64);
impl_try_from_primitive_type!(u128);

impl<'a, BaseRing: Ring, const N: usize> Sum<&'a Self> for Pow2CyclotomicPolyRing<BaseRing, N> {
    fn sum<I: Iterator<Item = &'a Self>>(iter: I) -> Self {
        iter.fold(Self::zero(), |acc, x| acc + x)
    }
}

impl<BaseRing: Ring, const N: usize> CanonicalSerialize for Pow2CyclotomicPolyRing<BaseRing, N> {
    fn serialize_with_mode<W: Write>(
        &self,
        writer: W,
        compress: Compress,
    ) -> Result<(), SerializationError> {
        self.0.serialize_with_mode(writer, compress)
    }

    fn serialized_size(&self, compress: Compress) -> usize {
        self.0.serialized_size(compress)
    }
}

impl<BaseRing: Ring, const N: usize> Valid for Pow2CyclotomicPolyRing<BaseRing, N> {
    fn check(&self) -> Result<(), SerializationError> {
        self.0.check()
    }
}

impl<BaseRing: Ring, const N: usize> CanonicalDeserialize for Pow2CyclotomicPolyRing<BaseRing, N> {
    fn deserialize_with_mode<R: Read>(
        reader: R,
        compress: Compress,
        validate: Validate,
    ) -> Result<Self, SerializationError> {
        Self::Inner::deserialize_with_mode(reader, compress, validate).map(Self)
    }
}

impl<BaseRing: Ring, const N: usize> Ring for Pow2CyclotomicPolyRing<BaseRing, N> {
    const ZERO: Self = Self::const_from_element(BaseRing::ZERO);
    const ONE: Self = Self::const_from_element(BaseRing::ONE);

    fn inverse(&self) -> Option<Self> {
        todo!()
    }
}

impl<BaseRing: Ring, const N: usize> FromRandomBytes<Self> for Pow2CyclotomicPolyRing<BaseRing, N> {
    fn needs_bytes() -> usize {
        N * BaseRing::byte_size()
    }

    fn try_from_random_bytes_inner(bytes: &[u8]) -> Option<Self> {
        let coeffs = core::array::from_fn(|i| {
            BaseRing::try_from_random_bytes(
                &bytes[i * BaseRing::byte_size()..(i + 1) * BaseRing::byte_size()],
            )
            .unwrap()
        });
        Some(Self::from(coeffs))
    }
}

// impl<BaseRing: ConvertibleField, const N: usize> Serialize for Pow2CyclotomicPolyRing<BaseRing, N> {
//     fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
//     where
//         S: Serializer,
//     {
//         // ark_se(&self.coeffs(), serializer)
//         // serialize(&self.0, serializer)
//         todo!()
//     }
// }

// impl<BaseRing: ConvertibleField, const N: usize> ToBytes for Pow2CyclotomicPolyRing<BaseRing, N>
// where
//     SVector<BaseRing, N>: ToBytes,
// {
//     type ToBytesError = <SVector<BaseRing, N> as ToBytes>::ToBytesError;
//
//     fn to_bytes(&self) -> Result<Vec<u8>, Self::ToBytesError> {
//         self.0.to_bytes()
//     }
// }

// impl<'de, BaseRing: ConvertibleField, const N: usize> Deserialize<'de>
//     for Pow2CyclotomicPolyRing<BaseRing, N>
// {
//     fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
//     where
//         D: Deserializer<'de>,
//     {
//         // ark_de(deserializer).map(|v: Vec<BaseRing>| Self::from(v))
//         // let coeffs: SVector<BaseRing, N> = Deserialize::deserialize(deserializer)?;
//         // Ok(Self(coeffs))
//         todo!()
//     }
// }

// impl<BaseRing: ConvertibleField, const N: usize> FromBytes for Pow2CyclotomicPolyRing<BaseRing, N>
// where
//     SVector<BaseRing, N>: FromBytes,
// {
//     type FromBytesError = <SVector<BaseRing, N> as FromBytes>::FromBytesError;
//
//     fn from_bytes(bytes: &[u8]) -> Result<Self, Self::FromBytesError> {
//         Self::Inner::from_bytes(bytes).map(Self)
//     }
// }

impl<BaseRing: Ring, const N: usize> Default for Pow2CyclotomicPolyRing<BaseRing, N> {
    #[inline(always)]
    fn default() -> Self {
        Self::zero()
    }
}

impl<BaseRing: Ring, const N: usize> Display for Pow2CyclotomicPolyRing<BaseRing, N> {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        std::fmt::Display::fmt(&self.0, f)
    }
}

impl<BaseRing: Ring, const N: usize> Zero for Pow2CyclotomicPolyRing<BaseRing, N> {
    #[inline(always)]
    fn zero() -> Self {
        Self::ZERO
    }

    #[inline(always)]
    fn is_zero(&self) -> bool {
        self.eq(&Self::ZERO)
    }
}

impl<BaseRing: Ring, const N: usize> One for Pow2CyclotomicPolyRing<BaseRing, N> {
    #[inline(always)]
    fn one() -> Self {
        Self::ONE
    }
}

impl<BaseRing: Ring, const N: usize> Mul<Self> for Pow2CyclotomicPolyRing<BaseRing, N> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        let mut out = vec![BaseRing::zero(); N];
        for i in 0..N {
            for j in 0..N - i {
                out[i + j] += self.0[i] * rhs.0[j];
            }
            for j in N - i..N {
                out[i + j - N] -= self.0[i] * rhs.0[j];
            }
        }
        Self::from(out)
    }
}

impl<BaseRing: Ring, const N: usize> Neg for Pow2CyclotomicPolyRing<BaseRing, N> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        self.0.neg().into()
    }
}

impl<BaseRing: Ring, const N: usize> UniformRand for Pow2CyclotomicPolyRing<BaseRing, N> {
    fn rand<R: Rng + ?Sized>(rng: &mut R) -> Self {
        Self::from_fn(|_| BaseRing::rand(rng))
    }
}

impl<BaseRing: Ring, const N: usize> MulAssign<Self> for Pow2CyclotomicPolyRing<BaseRing, N> {
    fn mul_assign(&mut self, rhs: Self) {
        let out = self.mul(rhs);
        self.0 = out.0;
    }
}

impl<'a, BaseRing: Ring, const N: usize> Add<&'a Self> for Pow2CyclotomicPolyRing<BaseRing, N> {
    type Output = Self;

    fn add(self, rhs: &'a Self) -> Self::Output {
        self.0.add(rhs.0).into()
    }
}

impl<'a, BaseRing: Ring, const N: usize> Sub<&'a Self> for Pow2CyclotomicPolyRing<BaseRing, N> {
    type Output = Self;

    fn sub(self, rhs: &'a Self) -> Self::Output {
        self.0.sub(rhs.0).into()
    }
}

impl<'a, BaseRing: Ring, const N: usize> Mul<&'a Self> for Pow2CyclotomicPolyRing<BaseRing, N> {
    type Output = Self;

    fn mul(self, rhs: &'a Self) -> Self::Output {
        self.mul(*rhs)
    }
}

impl<'a, BaseRing: Ring, const N: usize> AddAssign<&'a Self>
    for Pow2CyclotomicPolyRing<BaseRing, N>
{
    fn add_assign(&mut self, rhs: &'a Self) {
        self.0.add_assign(rhs.0)
    }
}

impl<'a, BaseRing: Ring, const N: usize> SubAssign<&'a Self>
    for Pow2CyclotomicPolyRing<BaseRing, N>
{
    fn sub_assign(&mut self, rhs: &'a Self) {
        self.0.sub_assign(rhs.0)
    }
}

impl<'a, BaseRing: Ring, const N: usize> MulAssign<&'a Self>
    for Pow2CyclotomicPolyRing<BaseRing, N>
{
    fn mul_assign(&mut self, rhs: &'a Self) {
        let out = self.mul(rhs);
        self.0 = out.0;
    }
}

impl<'a, BaseRing: Ring, const N: usize> Add<&'a mut Self> for Pow2CyclotomicPolyRing<BaseRing, N> {
    type Output = Self;

    fn add(self, rhs: &'a mut Self) -> Self::Output {
        self.0.add(rhs.0).into()
    }
}

impl<'a, BaseRing: Ring, const N: usize> Sub<&'a mut Self> for Pow2CyclotomicPolyRing<BaseRing, N> {
    type Output = Self;

    fn sub(self, rhs: &'a mut Self) -> Self::Output {
        self.0.sub(rhs.0).into()
    }
}

impl<'a, BaseRing: Ring, const N: usize> Mul<&'a mut Self> for Pow2CyclotomicPolyRing<BaseRing, N> {
    type Output = Self;

    fn mul(self, rhs: &'a mut Self) -> Self::Output {
        self.mul(*rhs)
    }
}

impl<'a, BaseRing: Ring, const N: usize> AddAssign<&'a mut Self>
    for Pow2CyclotomicPolyRing<BaseRing, N>
{
    fn add_assign(&mut self, rhs: &'a mut Self) {
        self.0.add_assign(rhs.0)
    }
}

impl<'a, BaseRing: Ring, const N: usize> SubAssign<&'a mut Self>
    for Pow2CyclotomicPolyRing<BaseRing, N>
{
    fn sub_assign(&mut self, rhs: &'a mut Self) {
        self.0.sub_assign(rhs.0)
    }
}

impl<'a, BaseRing: Ring, const N: usize> MulAssign<&'a mut Self>
    for Pow2CyclotomicPolyRing<BaseRing, N>
{
    fn mul_assign(&mut self, rhs: &'a mut Self) {
        let out = self.mul(rhs);
        self.0 = out.0;
    }
}

impl<BaseRing: Ring, const N: usize> Product<Self> for Pow2CyclotomicPolyRing<BaseRing, N> {
    fn product<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::one(), |a, b| a * b)
    }
}

impl<'a, BaseRing: Ring, const N: usize> Product<&'a Self> for Pow2CyclotomicPolyRing<BaseRing, N> {
    fn product<I: Iterator<Item = &'a Self>>(iter: I) -> Self {
        iter.fold(Self::one(), |a, b| a * b)
    }
}

impl<BaseRing: Ring, const N: usize> Mul<BaseRing> for Pow2CyclotomicPolyRing<BaseRing, N> {
    type Output = Self;

    fn mul(self, rhs: BaseRing) -> Self::Output {
        self.mul(Self::from_scalar(rhs))
    }
}

impl<BaseRing: Ring, const N: usize> PolyRing for Pow2CyclotomicPolyRing<BaseRing, N> {
    type BaseRing = BaseRing;
    fn coeffs(&self) -> Vec<Self::BaseRing> {
        self.0.into_iter().copied().collect()
    }
    fn dimension() -> usize {
        N
    }

    fn from_scalar(v: Self::BaseRing) -> Self {
        Self::from_fn(|i| if i == 0 { v } else { BaseRing::zero() })
    }
}

impl<BaseRing: Ring, const N: usize> From<Vec<BaseRing>> for Pow2CyclotomicPolyRing<BaseRing, N> {
    fn from(value: Vec<BaseRing>) -> Self {
        Self(SVector::<BaseRing, N>::try_from(value).unwrap())
    }
}

impl<BaseRing: Ring, const N: usize> WithConjugationAutomorphism
    for Pow2CyclotomicPolyRing<BaseRing, N>
{
    fn sigma(&self) -> Self {
        let coeffs = self.0 .0.as_slice();
        let mut new_coeffs = Vec::<BaseRing>::with_capacity(N);
        new_coeffs.push(coeffs[0]);
        new_coeffs.extend(
            coeffs[1..]
                .iter()
                .rev()
                .map(|v_i| -*v_i)
                .collect::<Vec<BaseRing>>(),
        );
        Self::from(new_coeffs)
    }
}

impl<BaseRing: Ring, const N: usize> WithL2Norm for Pow2CyclotomicPolyRing<BaseRing, N> {
    fn l2_norm_squared(&self) -> BigUint {
        self.coeffs().l2_norm_squared()
    }
}

impl<BaseRing: Ring, const N: usize> WithLinfNorm for Pow2CyclotomicPolyRing<BaseRing, N> {
    fn linf_norm(&self) -> BigUint {
        self.coeffs().linf_norm()
    }
}
