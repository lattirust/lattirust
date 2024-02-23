use std::fmt::{Debug, Display, Formatter};
use std::hash::Hash;
use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

use ark_serialize::{SerializationError, Valid};
use ark_std::UniformRand;
use derive_more::{Add, AddAssign, From, Into, Sub, SubAssign, Sum};
use nalgebra::{ArrayStorage, SVector};
use num_traits::{One, Zero};
use rand::Rng;
use serde::{Deserialize, Deserializer, Serialize, Serializer};

use crate::lattice_arithmetic::ntt::NTT;
use crate::lattice_arithmetic::poly_ring::PolyRing;
use crate::lattice_arithmetic::pow2_cyclotomic_poly_ring::Pow2CyclotomicPolyRing;
use crate::lattice_arithmetic::ring::{Fq, Ring, Zq};
use crate::lattice_arithmetic::traits::{FromRandomBytes, Modulus, Normed, WithConjugationAutomorphism};

#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash, Add, AddAssign, Sum, Sub, SubAssign, From, Into)]
pub struct Pow2CyclotomicPolyRingNTT<const Q: u64, const N: usize>(SVector<Zq<Q>, N>);

impl<const Q: u64, const N: usize> Pow2CyclotomicPolyRingNTT<Q, N> {
    pub fn from_slice(coeffs: &[Zq<Q>]) -> Self {
        assert_eq!(coeffs.len(), N);
        Self { 0: SVector::<Zq<Q>, N>::from_row_slice(coeffs) }
    }

    pub fn from_fn<F>(mut f: F) -> Self
        where F: FnMut(usize) -> Zq<Q> {
        let mut coeffs = (0..N).map(|i| f(i)).collect::<Vec<Zq<Q>>>();
        Self::ntt(coeffs.as_mut_slice());
        Self::from_slice(&coeffs)
    }

    pub fn from_value(v: Zq<Q>) -> Self {
        Self { 0: SVector::<Zq<Q>, N>::from_fn(|_i, _| v) }
    }
}

impl<const Q: u64, const N: usize> const NTT<Q, N> for Pow2CyclotomicPolyRingNTT<Q, N> {
    fn ntt_coeffs(&self) -> Vec<Fq<Q>> {
        self.coeffs().into_iter().map(|x| Into::into(x)).collect()
    }
}

impl<const Q: u64, const N: usize> Modulus for Pow2CyclotomicPolyRingNTT<Q, N> {
    fn modulus() -> u64 {
        Zq::<Q>::modulus()
    }
}

impl<const Q: u64, const N: usize> Valid for Pow2CyclotomicPolyRingNTT<Q, N> {
    fn check(&self) -> Result<(), SerializationError> {
        todo!()
    }
}

const fn vec_from_element<const Q: u64, const N: usize>(elem: Zq<Q>) -> SVector<Zq<Q>, N> {
    let coeffs = [elem; N];
    SVector::<Zq<Q>, N>::from_array_storage(ArrayStorage::<Zq<Q>, { N }, 1> { 0: [coeffs; 1] })
}

impl<const Q: u64, const N: usize> Ring for Pow2CyclotomicPolyRingNTT<Q, N> {
    const ZERO: Self = Self { 0: vec_from_element(Zq::<Q>::ZERO) };
    const ONE: Self = Self { 0: vec_from_element(Zq::<Q>::ONE) };

    fn inverse(&self) -> Option<Self> { None } // TODO: should we move this into a separate trait?
}

impl<const Q: u64, const N: usize> FromRandomBytes<Self> for Pow2CyclotomicPolyRingNTT<Q, N> {
    fn byte_size() -> usize {
        N * Zq::<Q>::byte_size()
    }

    fn from_random_bytes(bytes: &[u8]) -> Option<Self> {
        assert_eq!(bytes.len(), Self::byte_size());
        let coeffs = SVector::<Zq::<Q>, N>::from_fn(|i, _|
            Zq::<Q>::from_random_bytes(&bytes[i * Zq::<Q>::byte_size()..(i + 1) * Zq::<Q>::byte_size()]).unwrap()
        );
        Some(Self::from(coeffs))
    }
}

impl<const Q: u64, const N: usize> Serialize for Pow2CyclotomicPolyRingNTT<Q, N> {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error> where S: Serializer {
        self.0.serialize(serializer)
    }
}

impl<'a, const Q: u64, const N: usize> Deserialize<'a> for Pow2CyclotomicPolyRingNTT<Q, N> {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error> where D: Deserializer<'a> {
        Ok(Self { 0: Deserialize::deserialize(deserializer)? })
    }
}

impl<const Q: u64, const N: usize> Default for Pow2CyclotomicPolyRingNTT<Q, N> {
    #[inline(always)]
    fn default() -> Self {
        Self::zero()
    }
}

impl<const Q: u64, const N: usize> Display for Pow2CyclotomicPolyRingNTT<Q, N> {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result { std::fmt::Display::fmt(&self.0, f) }
}

impl<const Q: u64, const N: usize> Zero for Pow2CyclotomicPolyRingNTT<Q, N> {
    #[inline(always)]
    fn zero() -> Self {
        Self::ZERO
    }

    #[inline(always)]
    fn is_zero(&self) -> bool {
        self.eq(&Self::ZERO)
    }
}

impl<const Q: u64, const N: usize> One for Pow2CyclotomicPolyRingNTT<Q, N> {
    #[inline(always)]
    fn one() -> Self {
        Self::ONE
    }
}

impl<const Q: u64, const N: usize> Mul<Self, > for Pow2CyclotomicPolyRingNTT<Q, N> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        self.0.component_mul(&rhs.0).into()
    }
}

impl<const Q: u64, const N: usize> Neg<> for Pow2CyclotomicPolyRingNTT<Q, N> {
    type Output = Self;

    fn neg(self) -> Self::Output { self.0.neg().into() }
}

impl<const Q: u64, const N: usize> UniformRand for Pow2CyclotomicPolyRingNTT<Q, N> {
    fn rand<R: Rng + ?Sized>(rng: &mut R) -> Self {
        Self::from_fn(|_| Zq::<Q>::rand(rng))
    }
}

impl<const Q: u64, const N: usize> MulAssign<Self> for Pow2CyclotomicPolyRingNTT<Q, N> {
    fn mul_assign(&mut self, rhs: Self) {
        self.0.component_mul_assign(&rhs.0)
    }
}

impl<'a, const Q: u64, const N: usize> Add<&'a Self, > for Pow2CyclotomicPolyRingNTT<Q, N> {
    type Output = Self;

    fn add(self, rhs: &'a Self) -> Self::Output { self.0.add(&rhs.0).into() }
}

impl<'a, const Q: u64, const N: usize> Sub<&'a Self, > for Pow2CyclotomicPolyRingNTT<Q, N> {
    type Output = Self;

    fn sub(self, rhs: &'a Self) -> Self::Output { self.0.sub(rhs.0).into() }
}

impl<'a, const Q: u64, const N: usize> Mul<&'a Self, > for Pow2CyclotomicPolyRingNTT<Q, N> {
    type Output = Self;

    fn mul(self, rhs: &'a Self) -> Self::Output {
        self.0.component_mul(&rhs.0).into()
    }
}

impl<'a, const Q: u64, const N: usize> AddAssign<&'a Self> for Pow2CyclotomicPolyRingNTT<Q, N> {
    fn add_assign(&mut self, rhs: &'a Self) { self.0.add_assign(&rhs.0) }
}

impl<'a, const Q: u64, const N: usize> SubAssign<&'a Self> for Pow2CyclotomicPolyRingNTT<Q, N> {
    fn sub_assign(&mut self, rhs: &'a Self) { self.0.sub_assign(&rhs.0) }
}

impl<'a, const Q: u64, const N: usize> MulAssign<&'a Self> for Pow2CyclotomicPolyRingNTT<Q, N> {
    fn mul_assign(&mut self, rhs: &'a Self) { self.0.component_mul_assign(&rhs.0) }
}

impl<'a, const Q: u64, const N: usize> Add<&'a mut Self, > for Pow2CyclotomicPolyRingNTT<Q, N> {
    type Output = Self;

    fn add(self, rhs: &'a mut Self) -> Self::Output { self.0.add(&rhs.0).into() }
}

impl<'a, const Q: u64, const N: usize> Sub<&'a mut Self, > for Pow2CyclotomicPolyRingNTT<Q, N> {
    type Output = Self;

    fn sub(self, rhs: &'a mut Self) -> Self::Output { self.0.sub(&rhs.0).into() }
}

impl<'a, const Q: u64, const N: usize> Mul<&'a mut Self, > for Pow2CyclotomicPolyRingNTT<Q, N> {
    type Output = Self;

    fn mul(self, rhs: &'a mut Self) -> Self::Output { self.0.component_mul(&rhs.0).into() }
}

impl<'a, const Q: u64, const N: usize> AddAssign<&'a mut Self> for Pow2CyclotomicPolyRingNTT<Q, N> {
    fn add_assign(&mut self, rhs: &'a mut Self) {
        self.0.add_assign(&rhs.0).into()
    }
}

impl<'a, const Q: u64, const N: usize> SubAssign<&'a mut Self> for Pow2CyclotomicPolyRingNTT<Q, N> {
    fn sub_assign(&mut self, rhs: &'a mut Self) {
        self.0.sub_assign(&rhs.0).into()
    }
}

impl<'a, const Q: u64, const N: usize> MulAssign<&'a mut Self> for Pow2CyclotomicPolyRingNTT<Q, N> {
    fn mul_assign(&mut self, rhs: &'a mut Self) { self.0.component_mul_assign(&rhs.0) }
}

impl<const Q: u64, const N: usize> From<u128> for Pow2CyclotomicPolyRingNTT<Q, N> {
    fn from(value: u128) -> Self { Self::from_value(Zq::<Q>::from(value)) }
}

impl<const Q: u64, const N: usize> From<u64> for Pow2CyclotomicPolyRingNTT<Q, N> {
    fn from(value: u64) -> Self { Self::from_value(Zq::<Q>::from(value)) }
}

impl<const Q: u64, const N: usize> From<u32> for Pow2CyclotomicPolyRingNTT<Q, N> {
    fn from(value: u32) -> Self { Self::from_value(Zq::<Q>::from(value)) }
}

impl<const Q: u64, const N: usize> From<u16> for Pow2CyclotomicPolyRingNTT<Q, N> {
    fn from(value: u16) -> Self { Self::from_value(Zq::<Q>::from(value)) }
}

impl<const Q: u64, const N: usize> From<u8> for Pow2CyclotomicPolyRingNTT<Q, N> {
    fn from(value: u8) -> Self { Self::from_value(Zq::<Q>::from(value)) }
}

impl<const Q: u64, const N: usize> From<bool> for Pow2CyclotomicPolyRingNTT<Q, N> {
    fn from(value: bool) -> Self { Self::from_value(Zq::<Q>::from(value)) }
}

impl<const Q: u64, const N: usize> From<Pow2CyclotomicPolyRing<Zq<Q>, N>> for Pow2CyclotomicPolyRingNTT<Q, N> {
    fn from(value: Pow2CyclotomicPolyRing<Zq<Q>, N>) -> Self {
        let mut coeffs = value.coeffs();
        Self::ntt(coeffs.as_mut_slice());
        Self::from_slice(coeffs.as_slice())
    }
}

impl<const Q: u64, const N: usize> Mul<Zq<Q>> for Pow2CyclotomicPolyRingNTT<Q, N> {
    type Output = Self;

    fn mul(self, rhs: Zq::<Q>) -> Self::Output {
        self.mul(Self::from_value(rhs))
    }
}

impl<const Q: u64, const N: usize> PolyRing for Pow2CyclotomicPolyRingNTT<Q, N>
    where u64: From<Zq<Q>>, i64: From<Zq<Q>> {
    type BaseRing = Zq<Q>;
    fn coeffs(&self) -> Vec<Zq::<Q>> {
        self.0.iter().map(|v_i| Zq::<Q>::from(*v_i)).collect()
    }
    fn dimension() -> usize { N }
}

impl<const Q: u64, const N: usize> Normed<u64> for Pow2CyclotomicPolyRingNTT<Q, N> where u64: From<Zq<Q>>, i64: From<Zq<Q>> {
    fn norm(&self) -> f64 {
        (self.norm_squared() as f64).sqrt()
    }

    fn norm_squared(&self) -> u64 {
        let vec_i128 = SVector::<i128, N>::from_fn(|i, _| i64::from(self.0[i]) as i128);
        vec_i128.dot(&vec_i128) as u64
    }
}

impl<const Q: u64, const N: usize> From<Vec<Zq<Q>>> for Pow2CyclotomicPolyRingNTT<Q, N> {
    fn from(value: Vec<Zq<Q>>) -> Self {
        Self { 0: SVector::<Zq::<Q>, N>::from_vec(value) }
    }
}

impl<const Q: u64, const N: usize> From<Zq<Q>> for Pow2CyclotomicPolyRingNTT<Q, N> {
    fn from(value: Zq::<Q>) -> Self { Self::from_value(value) }
}

impl<const Q: u64, const N: usize> WithConjugationAutomorphism for Pow2CyclotomicPolyRingNTT<Q, N> {
    fn sigma(&self) -> Self {
        let coeffs = self.0.as_slice();
        let mut new_coeffs = Vec::<Zq<Q>>::with_capacity(N);
        new_coeffs.push(coeffs[0]);
        new_coeffs.extend(coeffs[1..].iter().rev().map(|v_i| -*v_i).collect::<Vec<Zq<Q>>>());
        Self::from(new_coeffs)
    }
}