use std::cmp::Ordering;
use std::fmt::{Debug, Display, Formatter};
use std::hash::Hash;
use std::io::{Read, Write};
use std::iter::Product;
use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

use ark_serialize::{CanonicalDeserialize, CanonicalDeserializeWithFlags, CanonicalSerialize, CanonicalSerializeWithFlags, Compress, Flags, SerializationError, Valid, Validate};
use ark_std::UniformRand;
use derive_more::{Add, AddAssign, From, Into, Sub, SubAssign, Sum};
use nalgebra::{ArrayStorage, SVector};
use num_traits::{One, Zero};
use rand::Rng;
use serde::{Deserialize, Deserializer, Serialize, Serializer};

use crate::lattice_arithmetic::poly_ring::PolyRing;
use crate::lattice_arithmetic::pow2_cyclotomic_poly_ring::Pow2CyclotomicPolyRing;
use crate::lattice_arithmetic::ring::Ring;
use crate::lattice_arithmetic::traits::{FromRandomBytes, IntegerDiv, Modulus, Normed, WithConjugationAutomorphism, WithLog2};

#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash, Add, AddAssign, Sum, Sub, SubAssign, From, Into)]
pub struct Pow2CyclotomicPolyRingNTT<BaseRing: Ring, const N: usize>(SVector<BaseRing, N>);

impl<BaseRing: Ring, const N: usize> Pow2CyclotomicPolyRingNTT<BaseRing, N> {
    pub fn from_fn<F>(mut f: F) -> Self
        where F: FnMut(usize) -> BaseRing {
        Self { 0: SVector::<BaseRing, N>::from_fn(|i, _| f(i)) }
    }
    pub fn from_value(v: BaseRing) -> Self {
        Self { 0: SVector::<BaseRing, N>::from_element(v) }
    }
}

impl<BaseRing: Ring, const N: usize> Valid for Pow2CyclotomicPolyRingNTT<BaseRing, N> {
    fn check(&self) -> Result<(), SerializationError> {
        todo!()
    }
}

const fn vec_from_element<BaseRing: Ring, const N: usize>(elem: BaseRing) -> SVector<BaseRing, N> {
    SVector::<BaseRing, N>::from_array_storage(ArrayStorage::<BaseRing, { N }, 1> { 0: [[elem; N]; 1] })
}

impl<BaseRing: Ring, const N: usize> Ring for Pow2CyclotomicPolyRingNTT<BaseRing, N> {
    const ZERO: Self = Self { 0: vec_from_element(BaseRing::ZERO) };
    const ONE: Self = Self { 0: vec_from_element(BaseRing::ONE) };
}

impl<BaseRing: Ring, const N: usize> FromRandomBytes<Self> for Pow2CyclotomicPolyRingNTT<BaseRing, N> {
    fn byte_size() -> usize {
        N * BaseRing::byte_size()
    }

    fn from_random_bytes(bytes: &[u8]) -> Option<Self> {
        assert_eq!(bytes.len(), Self::byte_size());
        let coeffs = SVector::<BaseRing, N>::from_fn(|i, _|
            BaseRing::from_random_bytes(&bytes[i * BaseRing::byte_size()..(i + 1) * BaseRing::byte_size()]).unwrap()
        );
        Some(Self::from(coeffs))
    }
}


impl<BaseRing: Ring, const N: usize> Serialize for Pow2CyclotomicPolyRingNTT<BaseRing, N> {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error> where S: Serializer {
        self.0.serialize(serializer)
    }
}

impl<'a, BaseRing: Ring, const N: usize> Deserialize<'a> for Pow2CyclotomicPolyRingNTT<BaseRing, N> {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error> where D: Deserializer<'a> {
        Ok(Self { 0: Deserialize::deserialize(deserializer)? })
    }
}


impl<BaseRing: Ring, const N: usize> CanonicalSerialize for Pow2CyclotomicPolyRingNTT<BaseRing, N> {
    fn serialize_with_mode<W: Write>(&self, writer: W, compress: Compress) -> Result<(), SerializationError> {
        todo!()
    }

    fn serialized_size(&self, compress: Compress) -> usize {
        self.0.map(|v_i| v_i.serialized_size(compress)).sum()
    }
}

impl<BaseRing: Ring, const N: usize> CanonicalDeserialize for Pow2CyclotomicPolyRingNTT<BaseRing, N> {
    fn deserialize_with_mode<R: Read>(reader: R, compress: Compress, validate: Validate) -> Result<Self, SerializationError> {
        todo!()
    }
}


impl<BaseRing: Ring, const N: usize> Default for Pow2CyclotomicPolyRingNTT<BaseRing, N> {
    fn default() -> Self {
        todo!()
    }
}

impl<BaseRing: Ring, const N: usize> Display for Pow2CyclotomicPolyRingNTT<BaseRing, N> {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        todo!()
    }
}


impl<BaseRing: Ring, const N: usize> Zero for Pow2CyclotomicPolyRingNTT<BaseRing, N> {
    fn zero() -> Self {
        Self::ZERO
    }

    fn is_zero(&self) -> bool {
        self.eq(&Self::ZERO)
    }
}

impl<BaseRing: Ring, const N: usize> One for Pow2CyclotomicPolyRingNTT<BaseRing, N> {
    fn one() -> Self {
        Self::ONE
    }
}

impl<BaseRing: Ring, const N: usize> Mul<Self, > for Pow2CyclotomicPolyRingNTT<BaseRing, N> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        Self(self.0.component_mul(&rhs.0))
    }
}

impl<BaseRing: Ring, const N: usize> Ord for Pow2CyclotomicPolyRingNTT<BaseRing, N> {
    fn cmp(&self, other: &Self) -> Ordering {
        unimplemented!()
    }
}

impl<BaseRing: Ring, const N: usize> PartialOrd<Self> for Pow2CyclotomicPolyRingNTT<BaseRing, N> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        self.0.partial_cmp(&other.0)
    }
}

impl<BaseRing: Ring, const N: usize> Neg<> for Pow2CyclotomicPolyRingNTT<BaseRing, N> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        todo!()
    }
}

impl<BaseRing: Ring, const N: usize> UniformRand for Pow2CyclotomicPolyRingNTT<BaseRing, N> {
    fn rand<R: Rng + ?Sized>(rng: &mut R) -> Self {
        Self::from_fn(|_| BaseRing::rand(rng))
    }
}

impl<BaseRing: Ring, const N: usize> CanonicalSerializeWithFlags for Pow2CyclotomicPolyRingNTT<BaseRing, N> {
    fn serialize_with_flags<W: Write, F: Flags>(&self, writer: W, flags: F) -> Result<(), SerializationError> {
        todo!()
    }

    fn serialized_size_with_flags<F: Flags>(&self) -> usize {
        todo!()
    }
}

impl<BaseRing: Ring, const N: usize> CanonicalDeserializeWithFlags for Pow2CyclotomicPolyRingNTT<BaseRing, N> {
    fn deserialize_with_flags<R: Read, F: Flags>(reader: R) -> Result<(Self, F), SerializationError> {
        todo!()
    }
}

// impl<BaseRing: Ring, const N: usize> Sub<Self, > for Pow2CyclotomicPolyRingNTT<BaseRing, N> {
//     type Output = Self;
//
//     fn sub(self, rhs: Self) -> Self::Output {
//         todo!()
//     }
// }

// impl<BaseRing: Ring, const N: usize> AddAssign<Self> for Pow2CyclotomicPolyRingNTT<BaseRing, N> {
//     fn add_assign(&mut self, rhs: Self) {
//         todo!()
//     }
// }

// impl<BaseRing: Ring, const N: usize> SubAssign<Self> for Pow2CyclotomicPolyRingNTT<BaseRing, N> {
//     fn sub_assign(&mut self, rhs: Self) {
//         todo!()
//     }
// }

impl<BaseRing: Ring, const N: usize> MulAssign<Self> for Pow2CyclotomicPolyRingNTT<BaseRing, N> {
    fn mul_assign(&mut self, rhs: Self) {
        self.0.component_mul_assign(&rhs.0)
    }
}

impl<'a, BaseRing: Ring, const N: usize> Add<&'a Self, > for Pow2CyclotomicPolyRingNTT<BaseRing, N> {
    type Output = Self;

    fn add(self, rhs: &'a Self) -> Self::Output { self.0.add(&rhs.0).into() }
}

impl<'a, BaseRing: Ring, const N: usize> Sub<&'a Self, > for Pow2CyclotomicPolyRingNTT<BaseRing, N> {
    type Output = Self;

    fn sub(self, rhs: &'a Self) -> Self::Output { self.0.sub(rhs.0).into() }
}

impl<'a, BaseRing: Ring, const N: usize> Mul<&'a Self, > for Pow2CyclotomicPolyRingNTT<BaseRing, N> {
    type Output = Self;

    fn mul(self, rhs: &'a Self) -> Self::Output {
        self.0.component_mul(&rhs.0).into()
    }
}

impl<'a, BaseRing: Ring, const N: usize> AddAssign<&'a Self> for Pow2CyclotomicPolyRingNTT<BaseRing, N> {
    fn add_assign(&mut self, rhs: &'a Self) {
        todo!()
    }
}

impl<'a, BaseRing: Ring, const N: usize> SubAssign<&'a Self> for Pow2CyclotomicPolyRingNTT<BaseRing, N> {
    fn sub_assign(&mut self, rhs: &'a Self) {
        todo!()
    }
}

impl<'a, BaseRing: Ring, const N: usize> MulAssign<&'a Self> for Pow2CyclotomicPolyRingNTT<BaseRing, N> {
    fn mul_assign(&mut self, rhs: &'a Self) {
        todo!()
    }
}

impl<'a, BaseRing: Ring, const N: usize> Add<&'a mut Self, > for Pow2CyclotomicPolyRingNTT<BaseRing, N> {
    type Output = Self;

    fn add(self, rhs: &'a mut Self) -> Self::Output {
        todo!()
    }
}

impl<'a, BaseRing: Ring, const N: usize> Sub<&'a mut Self, > for Pow2CyclotomicPolyRingNTT<BaseRing, N> {
    type Output = Self;

    fn sub(self, rhs: &'a mut Self) -> Self::Output {
        todo!()
    }
}

impl<'a, BaseRing: Ring, const N: usize> Mul<&'a mut Self, > for Pow2CyclotomicPolyRingNTT<BaseRing, N> {
    type Output = Self;

    fn mul(self, rhs: &'a mut Self) -> Self::Output {
        todo!()
    }
}

impl<'a, BaseRing: Ring, const N: usize> AddAssign<&'a mut Self> for Pow2CyclotomicPolyRingNTT<BaseRing, N> {
    fn add_assign(&mut self, rhs: &'a mut Self) {
        todo!()
    }
}

impl<'a, BaseRing: Ring, const N: usize> SubAssign<&'a mut Self> for Pow2CyclotomicPolyRingNTT<BaseRing, N> {
    fn sub_assign(&mut self, rhs: &'a mut Self) {
        todo!()
    }
}

impl<'a, BaseRing: Ring, const N: usize> MulAssign<&'a mut Self> for Pow2CyclotomicPolyRingNTT<BaseRing, N> {
    fn mul_assign(&mut self, rhs: &'a mut Self) {
        todo!()
    }
}

/*
impl<BaseRing: Ring, const N: usize> Sum<Self> for Pow2CyclotomicPolyRingNTT<BaseRing, N> {
    fn sum<I: Iterator<Item=Self>>(iter: I) -> Self {
        todo!()
    }
}

impl<'a, BaseRing: Ring, const N: usize> Sum<&'a Self> for Pow2CyclotomicPolyRingNTT<BaseRing, N> {
    fn sum<I: Iterator<Item=&'a Self>>(iter: I) -> Self {
        todo!()
    }
}

 */

impl<BaseRing: Ring, const N: usize> Product<Self> for Pow2CyclotomicPolyRingNTT<BaseRing, N> {
    fn product<I: Iterator<Item=Self>>(iter: I) -> Self {
        todo!()
    }
}

impl<'a, BaseRing: Ring, const N: usize> Product<&'a Self> for Pow2CyclotomicPolyRingNTT<BaseRing, N> {
    fn product<I: Iterator<Item=&'a Self>>(iter: I) -> Self {
        todo!()
    }
}

impl<BaseRing: Ring, const N: usize> From<u128> for Pow2CyclotomicPolyRingNTT<BaseRing, N> {
    fn from(value: u128) -> Self { Self::from_value(BaseRing::from(value)) }
}

impl<BaseRing: Ring, const N: usize> From<u64> for Pow2CyclotomicPolyRingNTT<BaseRing, N> {
    fn from(value: u64) -> Self { Self::from_value(BaseRing::from(value)) }
}

impl<BaseRing: Ring, const N: usize> From<u32> for Pow2CyclotomicPolyRingNTT<BaseRing, N> {
    fn from(value: u32) -> Self { Self::from_value(BaseRing::from(value)) }
}

impl<BaseRing: Ring, const N: usize> From<u16> for Pow2CyclotomicPolyRingNTT<BaseRing, N> {
    fn from(value: u16) -> Self { Self::from_value(BaseRing::from(value)) }
}

impl<BaseRing: Ring, const N: usize> From<u8> for Pow2CyclotomicPolyRingNTT<BaseRing, N> {
    fn from(value: u8) -> Self { Self::from_value(BaseRing::from(value)) }
}

impl<BaseRing: Ring, const N: usize> From<bool> for Pow2CyclotomicPolyRingNTT<BaseRing, N> {
    fn from(value: bool) -> Self { Self::from_value(BaseRing::from(value)) }
}

impl<BaseRing: Ring, const N: usize> Mul<BaseRing> for Pow2CyclotomicPolyRingNTT<BaseRing, N> {
    type Output = Self;

    fn mul(self, rhs: BaseRing) -> Self::Output {
        self.mul(Self::from_value(rhs))
    }
}

impl<BaseRing: Ring + IntegerDiv + WithLog2 + Modulus, const N: usize> PolyRing for Pow2CyclotomicPolyRingNTT<BaseRing, N>
    where i64: From<BaseRing> {
    type BaseRing = BaseRing;
    fn coeffs(&self) -> Vec<Self::BaseRing> {
        self.0.iter().map(|v_i| BaseRing::from(*v_i)).collect()
    }
}

impl<BaseRing: Ring, const N: usize> Normed<u64> for Pow2CyclotomicPolyRingNTT<BaseRing, N> where i64: From<BaseRing> {
    fn norm(&self) -> f64 {
        (self.norm_squared() as f64).sqrt()
    }

    fn norm_squared(&self) -> u64 {
        let vec_i128 = SVector::<i128, N>::from_fn(|i, _| i64::from(self.0[i]) as i128);
        vec_i128.dot(&vec_i128) as u64
    }
}

impl<BaseRing: Ring, const N: usize> From<Vec<BaseRing>> for Pow2CyclotomicPolyRingNTT<BaseRing, N> {
    fn from(value: Vec<BaseRing>) -> Self {
        Self { 0: SVector::<BaseRing, N>::from_vec(value) }
    }
}

impl<BaseRing: Ring, const N: usize> From<BaseRing> for Pow2CyclotomicPolyRingNTT<BaseRing, N> {
    fn from(value: BaseRing) -> Self { Self::from_value(value) }
}

impl<BaseRing: Ring, const N: usize> WithConjugationAutomorphism for Pow2CyclotomicPolyRingNTT<BaseRing, N> {
    fn sigma(&self) -> Self {
        todo!()
    }
}