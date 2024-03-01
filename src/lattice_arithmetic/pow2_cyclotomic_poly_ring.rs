use std::fmt::{Debug, Display, Formatter};
use std::hash::Hash;
use std::iter::Product;
use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};
use ark_ff::Field;

use ark_serialize::{SerializationError, Valid};
use ark_std::UniformRand;
use derive_more::{Add, AddAssign, From, Into, Sub, SubAssign, Sum};
use nalgebra::{ArrayStorage, SVector};
use num_traits::{One, Zero};
use rand::Rng;
use serde::{Deserialize, Deserializer, Serialize, Serializer};

use crate::lattice_arithmetic::poly_ring::{ConvertibleField, PolyRing};
use crate::lattice_arithmetic::ring::Ring;
use crate::lattice_arithmetic::traits::{FromRandomBytes, IntegerDiv, Modulus, WithL2Norm, WithConjugationAutomorphism, WithLinfNorm};

#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash, Add, AddAssign, Sum, Sub, SubAssign, From, Into)]
pub struct Pow2CyclotomicPolyRing<BaseRing: ConvertibleField, const N: usize>(SVector<BaseRing, N>);

impl<BaseRing: ConvertibleField, const N: usize> Pow2CyclotomicPolyRing<BaseRing, N> {
    pub fn from_fn<F>(mut f: F) -> Self
        where F: FnMut(usize) -> BaseRing {
        Self { 0: SVector::<BaseRing, N>::from_fn(|i, _| f(i)) }
    }
}

impl<BaseRing: ConvertibleField, const N: usize> Valid for Pow2CyclotomicPolyRing<BaseRing, N> {
    fn check(&self) -> Result<(), SerializationError> {
        todo!()
    }
}

const fn vec_from_element<BaseRing: ConvertibleField, const N: usize>(elem: BaseRing) -> SVector<BaseRing, N> {
    let mut coeffs = [BaseRing::ZERO; N];
    coeffs[0] = elem;
    SVector::<BaseRing, N>::from_array_storage(ArrayStorage::<BaseRing, { N }, 1> { 0: [coeffs; 1] })
}

impl<BaseRing: ConvertibleField, const N: usize> Modulus for Pow2CyclotomicPolyRing<BaseRing, N> {
    fn modulus() -> u64 {
       todo!()
    }
}

impl<BaseRing: ConvertibleField, const N: usize> Ring for Pow2CyclotomicPolyRing<BaseRing, N> {
    const ZERO: Self = Self { 0: vec_from_element(BaseRing::ZERO) };
    const ONE: Self = Self { 0: vec_from_element(BaseRing::ONE) };

    fn inverse(&self) -> Option<Self> { None } // TODO: should we move this into a separate trait?
}

impl<BaseRing: ConvertibleField, const N: usize> FromRandomBytes<Self> for Pow2CyclotomicPolyRing<BaseRing, N> {
    fn byte_size() -> usize {
        N * BaseRing::byte_size()
    }

    fn try_from_random_bytes(bytes: &[u8]) -> Option<Self> {
        assert_eq!(bytes.len(), Self::byte_size());
        let coeffs = SVector::<BaseRing, N>::from_fn(|i, _|
            BaseRing::try_from_random_bytes(&bytes[i * BaseRing::byte_size()..(i + 1) * BaseRing::byte_size()]).unwrap()
        );
        Some(Self::from(coeffs))
    }
}

impl<BaseRing: ConvertibleField, const N: usize> Serialize for Pow2CyclotomicPolyRing<BaseRing, N> {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error> where S: Serializer {
        // self.0.serialize(serializer)
        todo!()
    }
}

impl<'a, BaseRing: ConvertibleField, const N: usize> Deserialize<'a> for Pow2CyclotomicPolyRing<BaseRing, N> {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error> where D: Deserializer<'a> {
        // Ok(Self { 0: Deserialize::deserialize(deserializer)? })
        todo!()
    }
}

impl<BaseRing: ConvertibleField, const N: usize> Default for Pow2CyclotomicPolyRing<BaseRing, N> {
    #[inline(always)]
    fn default() -> Self {
        Self::zero()
    }
}

impl<BaseRing: ConvertibleField, const N: usize> Display for Pow2CyclotomicPolyRing<BaseRing, N> {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result { std::fmt::Display::fmt(&self.0, f) }
}

impl<BaseRing: ConvertibleField, const N: usize> Zero for Pow2CyclotomicPolyRing<BaseRing, N> {
    #[inline(always)]
    fn zero() -> Self {
        Self::ZERO
    }

    #[inline(always)]
    fn is_zero(&self) -> bool {
        self.eq(&Self::ZERO)
    }
}

impl<BaseRing: ConvertibleField, const N: usize> One for Pow2CyclotomicPolyRing<BaseRing, N> {
    #[inline(always)]
    fn one() -> Self {
        Self::ONE
    }
}

impl<BaseRing: ConvertibleField, const N: usize> Mul<Self, > for Pow2CyclotomicPolyRing<BaseRing, N> {
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

impl<BaseRing: ConvertibleField, const N: usize> Neg<> for Pow2CyclotomicPolyRing<BaseRing, N> {
    type Output = Self;

    fn neg(self) -> Self::Output { self.0.neg().into() }
}

impl<BaseRing: ConvertibleField, const N: usize> UniformRand for Pow2CyclotomicPolyRing<BaseRing, N> {
    fn rand<R: Rng + ?Sized>(rng: &mut R) -> Self {
        Self::from_fn(|_| BaseRing::rand(rng))
    }
}

impl<BaseRing: ConvertibleField, const N: usize> MulAssign<Self> for Pow2CyclotomicPolyRing<BaseRing, N> {
    fn mul_assign(&mut self, rhs: Self) {
        let out = self.mul(rhs);
        self.0 = out.0;
    }
}

impl<'a, BaseRing: ConvertibleField, const N: usize> Add<&'a Self, > for Pow2CyclotomicPolyRing<BaseRing, N> {
    type Output = Self;

    fn add(self, rhs: &'a Self) -> Self::Output { self.0.add(&rhs.0).into() }
}

impl<'a, BaseRing: ConvertibleField, const N: usize> Sub<&'a Self, > for Pow2CyclotomicPolyRing<BaseRing, N> {
    type Output = Self;

    fn sub(self, rhs: &'a Self) -> Self::Output { self.0.sub(rhs.0).into() }
}

impl<'a, BaseRing: ConvertibleField, const N: usize> Mul<&'a Self, > for Pow2CyclotomicPolyRing<BaseRing, N> {
    type Output = Self;

    fn mul(self, rhs: &'a Self) -> Self::Output {
        self.mul(*rhs)
    }
}

impl<'a, BaseRing: ConvertibleField, const N: usize> AddAssign<&'a Self> for Pow2CyclotomicPolyRing<BaseRing, N> {
    fn add_assign(&mut self, rhs: &'a Self) { self.0.add_assign(&rhs.0) }
}

impl<'a, BaseRing: ConvertibleField, const N: usize> SubAssign<&'a Self> for Pow2CyclotomicPolyRing<BaseRing, N> {
    fn sub_assign(&mut self, rhs: &'a Self) { self.0.sub_assign(&rhs.0) }
}

impl<'a, BaseRing: ConvertibleField, const N: usize> MulAssign<&'a Self> for Pow2CyclotomicPolyRing<BaseRing, N> {
    fn mul_assign(&mut self, rhs: &'a Self) {
        let out = self.mul(rhs);
        self.0 = out.0;
    }
}

impl<'a, BaseRing: ConvertibleField, const N: usize> Add<&'a mut Self, > for Pow2CyclotomicPolyRing<BaseRing, N> {
    type Output = Self;

    fn add(self, rhs: &'a mut Self) -> Self::Output { self.0.add(&rhs.0).into() }
}

impl<'a, BaseRing: ConvertibleField, const N: usize> Sub<&'a mut Self, > for Pow2CyclotomicPolyRing<BaseRing, N> {
    type Output = Self;

    fn sub(self, rhs: &'a mut Self) -> Self::Output { self.0.sub(&rhs.0).into() }
}

impl<'a, BaseRing: ConvertibleField, const N: usize> Mul<&'a mut Self, > for Pow2CyclotomicPolyRing<BaseRing, N> {
    type Output = Self;

    fn mul(self, rhs: &'a mut Self) -> Self::Output {
        self.mul(*rhs)
    }
}

impl<'a, BaseRing: ConvertibleField, const N: usize> AddAssign<&'a mut Self> for Pow2CyclotomicPolyRing<BaseRing, N> {
    fn add_assign(&mut self, rhs: &'a mut Self) {
        self.0.add_assign(&rhs.0).into()
    }
}

impl<'a, BaseRing: ConvertibleField, const N: usize> SubAssign<&'a mut Self> for Pow2CyclotomicPolyRing<BaseRing, N> {
    fn sub_assign(&mut self, rhs: &'a mut Self) {
        self.0.sub_assign(&rhs.0).into()
    }
}

impl<'a, BaseRing: ConvertibleField, const N: usize> MulAssign<&'a mut Self> for Pow2CyclotomicPolyRing<BaseRing, N> {
    fn mul_assign(&mut self, rhs: &'a mut Self) {
        let out = self.mul(rhs);
        self.0 = out.0;
    }
}

impl<BaseRing: ConvertibleField, const N: usize> Product<Self> for Pow2CyclotomicPolyRing<BaseRing, N> {
    fn product<I: Iterator<Item=Self>>(iter: I) -> Self { iter.fold(Self::one(), |a, b| a * b) }
}

impl<'a, BaseRing: ConvertibleField, const N: usize> Product<&'a Self> for Pow2CyclotomicPolyRing<BaseRing, N> {
    fn product<I: Iterator<Item=&'a Self>>(iter: I) -> Self {
        iter.fold(Self::one(), |a, b| a * b)
    }
}

impl<BaseRing: ConvertibleField, const N: usize> Mul<BaseRing> for Pow2CyclotomicPolyRing<BaseRing, N> {
    type Output = Self;

    fn mul(self, rhs: BaseRing) -> Self::Output {
        self.mul(Self::from_scalar(rhs))
    }
}

impl<BaseRing: ConvertibleField, const N: usize> PolyRing for Pow2CyclotomicPolyRing<BaseRing, N> {
    type BaseRing = BaseRing;
    fn coeffs(&self) -> Vec<Self::BaseRing> {
        self.0.iter().map(|v_i| BaseRing::from(*v_i)).collect()
    }
    fn dimension() -> usize { N }

    fn from_scalar(v: Self::BaseRing) -> Self {
        Self { 0: SVector::<BaseRing, N>::from_fn(|i, _| if i == 0 { v } else { BaseRing::zero() }) }
    }

    fn signed_repr(x: &Self::BaseRing) -> i128 {
        todo!()
    }

    fn unsigned_repr(x: &Self::BaseRing) -> u128 {
        todo!()
    }
}

impl<BaseRing: ConvertibleField, const N: usize> From<Vec<BaseRing>> for Pow2CyclotomicPolyRing<BaseRing, N> {
    fn from(value: Vec<BaseRing>) -> Self {
        Self { 0: SVector::<BaseRing, N>::from_vec(value) }
    }
}

impl<BaseRing: ConvertibleField, const N: usize> From<BaseRing> for Pow2CyclotomicPolyRing<BaseRing, N> {
    fn from(value: BaseRing) -> Self {Self::from_scalar(value)}
}

impl<BaseRing: ConvertibleField, const N: usize> From<u128> for Pow2CyclotomicPolyRing<BaseRing, N> {
    fn from(value: u128) -> Self {Self::from_scalar(value.into())}
}

impl<BaseRing: ConvertibleField, const N: usize> WithConjugationAutomorphism for Pow2CyclotomicPolyRing<BaseRing, N> {
    fn sigma(&self) -> Self {
        let coeffs = self.0.as_slice();
        let mut new_coeffs = Vec::<BaseRing>::with_capacity(N);
        new_coeffs.push(coeffs[0]);
        new_coeffs.extend(coeffs[1..].iter().rev().map(|v_i| -*v_i).collect::<Vec<BaseRing>>());
        Self::from(new_coeffs)
    }
}

impl<BaseRing: ConvertibleField, const N: usize>  WithL2Norm  for Pow2CyclotomicPolyRing<BaseRing, N> {
    fn l2_norm_squared(&self) -> u64 {
        self.coeffs().l2_norm_squared()
    }
}

impl<BaseRing: ConvertibleField, const N: usize>  WithLinfNorm  for Pow2CyclotomicPolyRing<BaseRing, N> {
    fn linf_norm(&self) -> u128 {
        self.coeffs().linf_norm()
    }
}