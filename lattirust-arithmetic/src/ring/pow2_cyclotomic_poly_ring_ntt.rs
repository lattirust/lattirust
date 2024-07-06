use std::fmt::{Debug, Display, Formatter};
use std::hash::Hash;
use std::io::{Read, Write};
use std::iter::{Product, Sum};
use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

use ark_serialize::{
    CanonicalDeserialize, CanonicalSerialize, Compress, SerializationError, Valid, Validate,
};
use ark_std::rand::Rng;
use ark_std::UniformRand;
use derive_more::{Add, AddAssign, From, Into, Sub, SubAssign, Sum};
use num_bigint::BigUint;
use num_traits::{One, Zero};

use crate::linear_algebra::{Matrix, SVector, Vector};
use crate::ntt::NTT;
use crate::ring::pow2_cyclotomic_poly_ring::Pow2CyclotomicPolyRing;
use crate::ring::PolyRing;
use crate::ring::{Ring, Zq};
use crate::traits::{
    FromRandomBytes, Modulus, WithConjugationAutomorphism, WithL2Norm, WithLinfNorm,
};

use super::poly_ring::WithRot;

#[derive(
    Clone, Copy, Debug, Eq, PartialEq, Hash, Add, AddAssign, Sum, Sub, SubAssign, From, Into,
)]
pub struct Pow2CyclotomicPolyRingNTT<const Q: u64, const N: usize>(SVector<Zq<Q>, N>);

impl<const Q: u64, const N: usize> Pow2CyclotomicPolyRingNTT<Q, N> {
    #[allow(dead_code)]
    pub(crate) type Inner = SVector<Zq<Q>, N>;

    /// Constructs a polynomial from an array of coefficients in NTT form.
    /// This function is private since we can't enforce coeffs being in NTT form if called from outside.
    fn from_array(coeffs_ntt: [Zq<Q>; N]) -> Self {
        Self(Self::Inner::const_from_array(coeffs_ntt))
    }

    /// Constructs a polynomial from a function specifying coefficients in non-NTT form.
    pub fn from_fn<F>(f: F) -> Self
    where
        F: FnMut(usize) -> Zq<Q>,
    {
        let mut coeffs = core::array::from_fn(f);
        Self::ntt(&mut coeffs);
        Self::from_array(coeffs)
    }
}

impl<const Q: u64, const N: usize> const NTT<Q, N> for Pow2CyclotomicPolyRingNTT<Q, N> {
    fn ntt_coeffs(&self) -> Vec<Zq<Q>> {
        self.0.iter().cloned().collect()
    }
}

impl<const Q: u64, const N: usize> Modulus for Pow2CyclotomicPolyRingNTT<Q, N> {
    fn modulus() -> BigUint {
        Zq::<Q>::modulus()
    }
}

const fn vec_from_element<const Q: u64, const N: usize>(elem: Zq<Q>) -> SVector<Zq<Q>, N> {
    SVector::<Zq<Q>, N>::const_from_array([elem; N])
}

impl<const Q: u64, const N: usize> CanonicalSerialize for Pow2CyclotomicPolyRingNTT<Q, N> {
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

impl<const Q: u64, const N: usize> Valid for Pow2CyclotomicPolyRingNTT<Q, N> {
    fn check(&self) -> Result<(), SerializationError> {
        self.0.check()
    }
}

impl<const Q: u64, const N: usize> CanonicalDeserialize for Pow2CyclotomicPolyRingNTT<Q, N> {
    fn deserialize_with_mode<R: Read>(
        reader: R,
        compress: Compress,
        validate: Validate,
    ) -> Result<Self, SerializationError> {
        Self::Inner::deserialize_with_mode(reader, compress, validate).map(Self)
    }
}

impl<const Q: u64, const N: usize> Ring for Pow2CyclotomicPolyRingNTT<Q, N> {
    const ZERO: Self = Self {
        0: vec_from_element(<Zq<Q> as Ring>::ZERO),
    };
    const ONE: Self = Self {
        0: vec_from_element(<Zq<Q> as Ring>::ONE),
    };
}

impl<const Q: u64, const N: usize> FromRandomBytes<Self> for Pow2CyclotomicPolyRingNTT<Q, N> {
    fn byte_size() -> usize {
        N * Zq::<Q>::byte_size()
    }

    fn try_from_random_bytes(bytes: &[u8]) -> Option<Self> {
        assert_eq!(bytes.len(), Self::byte_size());
        let coeffs = core::array::from_fn(|i| {
            Zq::<Q>::try_from_random_bytes(
                &bytes[i * Zq::<Q>::byte_size()..(i + 1) * Zq::<Q>::byte_size()],
            )
            .unwrap()
        });
        Some(Self::from_array(coeffs))
    }
}

// impl<const Q: u64, const N: usize> Serialize for Pow2CyclotomicPolyRingNTT<Q, N> {
//     fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
//     where
//         S: Serializer,
//     {
//         ark_se(&self.coeffs(), serializer)
//     }
// }
//
// impl<'a, const Q: u64, const N: usize> Deserialize<'a> for Pow2CyclotomicPolyRingNTT<Q, N> {
//     fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
//     where
//         D: Deserializer<'a>,
//     {
//         ark_de(deserializer).map(|v: Vec<Fq<Q>>| Self::from(v))
//     }
// }

impl<const Q: u64, const N: usize> Default for Pow2CyclotomicPolyRingNTT<Q, N> {
    #[inline(always)]
    fn default() -> Self {
        Self::zero()
    }
}

impl<const Q: u64, const N: usize> Display for Pow2CyclotomicPolyRingNTT<Q, N> {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        std::fmt::Display::fmt(&self.0, f)
    }
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

impl<const Q: u64, const N: usize> Mul<Self> for Pow2CyclotomicPolyRingNTT<Q, N> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        self.0.component_mul(&rhs.0).into()
    }
}

impl<const Q: u64, const N: usize> Neg for Pow2CyclotomicPolyRingNTT<Q, N> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        self.0.neg().into()
    }
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

impl<'a, const Q: u64, const N: usize> Add<&'a Self> for Pow2CyclotomicPolyRingNTT<Q, N> {
    type Output = Self;

    fn add(self, rhs: &'a Self) -> Self::Output {
        self.0.add(rhs.0).into()
    }
}

impl<'a, const Q: u64, const N: usize> Sub<&'a Self> for Pow2CyclotomicPolyRingNTT<Q, N> {
    type Output = Self;

    fn sub(self, rhs: &'a Self) -> Self::Output {
        self.0.sub(rhs.0).into()
    }
}

impl<'a, const Q: u64, const N: usize> Mul<&'a Self> for Pow2CyclotomicPolyRingNTT<Q, N> {
    type Output = Self;

    fn mul(self, rhs: &'a Self) -> Self::Output {
        self.0.component_mul(&rhs.0).into()
    }
}

impl<'a, const Q: u64, const N: usize> AddAssign<&'a Self> for Pow2CyclotomicPolyRingNTT<Q, N> {
    fn add_assign(&mut self, rhs: &'a Self) {
        self.0.add_assign(rhs.0)
    }
}

impl<'a, const Q: u64, const N: usize> SubAssign<&'a Self> for Pow2CyclotomicPolyRingNTT<Q, N> {
    fn sub_assign(&mut self, rhs: &'a Self) {
        self.0.sub_assign(rhs.0)
    }
}

impl<'a, const Q: u64, const N: usize> MulAssign<&'a Self> for Pow2CyclotomicPolyRingNTT<Q, N> {
    fn mul_assign(&mut self, rhs: &'a Self) {
        self.0.component_mul_assign(&rhs.0)
    }
}

impl<'a, const Q: u64, const N: usize> Add<&'a mut Self> for Pow2CyclotomicPolyRingNTT<Q, N> {
    type Output = Self;

    fn add(self, rhs: &'a mut Self) -> Self::Output {
        self.0.add(rhs.0).into()
    }
}

impl<'a, const Q: u64, const N: usize> Sub<&'a mut Self> for Pow2CyclotomicPolyRingNTT<Q, N> {
    type Output = Self;

    fn sub(self, rhs: &'a mut Self) -> Self::Output {
        self.0.sub(rhs.0).into()
    }
}

impl<'a, const Q: u64, const N: usize> Mul<&'a mut Self> for Pow2CyclotomicPolyRingNTT<Q, N> {
    type Output = Self;

    fn mul(self, rhs: &'a mut Self) -> Self::Output {
        self.0.component_mul(&rhs.0).into()
    }
}

impl<'a, const Q: u64, const N: usize> AddAssign<&'a mut Self> for Pow2CyclotomicPolyRingNTT<Q, N> {
    fn add_assign(&mut self, rhs: &'a mut Self) {
        self.0.add_assign(rhs.0).into()
    }
}

impl<'a, const Q: u64, const N: usize> SubAssign<&'a mut Self> for Pow2CyclotomicPolyRingNTT<Q, N> {
    fn sub_assign(&mut self, rhs: &'a mut Self) {
        self.0.sub_assign(rhs.0).into()
    }
}

impl<'a, const Q: u64, const N: usize> MulAssign<&'a mut Self> for Pow2CyclotomicPolyRingNTT<Q, N> {
    fn mul_assign(&mut self, rhs: &'a mut Self) {
        self.0.component_mul_assign(&rhs.0)
    }
}

impl<const Q: u64, const N: usize> From<Pow2CyclotomicPolyRing<Zq<Q>, N>>
    for Pow2CyclotomicPolyRingNTT<Q, N>
{
    fn from(value: Pow2CyclotomicPolyRing<Zq<Q>, N>) -> Self {
        let mut coeffs: [Zq<Q>; N] = value.coeffs().try_into().unwrap();
        Self::ntt(&mut coeffs);
        Self::from_array(coeffs)
    }
}

impl<const Q: u64, const N: usize> Into<Pow2CyclotomicPolyRing<Zq<Q>, N>>
    for Pow2CyclotomicPolyRingNTT<Q, N>
{
    fn into(self) -> Pow2CyclotomicPolyRing<Zq<Q>, N> {
        Pow2CyclotomicPolyRing::<Zq<Q>, N>::from(self.coeffs())
    }
}

impl<const Q: u64, const N: usize> Mul<Zq<Q>> for Pow2CyclotomicPolyRingNTT<Q, N> {
    type Output = Self;

    fn mul(self, rhs: Zq<Q>) -> Self::Output {
        self.mul(Self::from_scalar(rhs))
    }
}

macro_rules! impl_from_primitive_type {
    ($primitive_type: ty) => {
        impl<const Q: u64, const N: usize> From<$primitive_type>
            for Pow2CyclotomicPolyRingNTT<Q, N>
        {
            fn from(value: $primitive_type) -> Self {
                Self::from_scalar(Zq::<Q>::from(value))
            }
        }
    };
}

impl_from_primitive_type!(u128);
impl_from_primitive_type!(u64);
impl_from_primitive_type!(u32);
impl_from_primitive_type!(u16);
impl_from_primitive_type!(u8);
impl_from_primitive_type!(bool);

impl<'a, const Q: u64, const N: usize> Sum<&'a Self> for Pow2CyclotomicPolyRingNTT<Q, N> {
    fn sum<I: Iterator<Item = &'a Self>>(iter: I) -> Self {
        iter.fold(Self::zero(), |acc, x| acc + x)
    }
}

impl<const Q: u64, const N: usize> Product<Self> for Pow2CyclotomicPolyRingNTT<Q, N> {
    fn product<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::one(), |acc, x| acc * x)
    }
}

impl<'a, const Q: u64, const N: usize> Product<&'a Self> for Pow2CyclotomicPolyRingNTT<Q, N> {
    fn product<I: Iterator<Item = &'a Self>>(iter: I) -> Self {
        iter.fold(Self::one(), |acc, x| acc * x)
    }
}

impl<const Q: u64, const N: usize> PolyRing for Pow2CyclotomicPolyRingNTT<Q, N> {
    type BaseRing = Zq<Q>;
    fn coeffs(&self) -> Vec<Zq<Q>> {
        let mut coeffs = self.ntt_coeffs().try_into().unwrap();
        Self::intt(&mut coeffs);
        coeffs.to_vec()
    }
    fn dimension() -> usize {
        N
    }

    fn from_scalar(v: Self::BaseRing) -> Self {
        // NTT([v, 0, ..., 0]) = ([v, ..., v])
        Self::from_array([v; N])
    }
}

impl<const Q: u64, const N: usize> From<Vec<Zq<Q>>> for Pow2CyclotomicPolyRingNTT<Q, N> {
    fn from(value: Vec<Zq<Q>>) -> Self {
        let mut array = TryInto::<[Zq<Q>; N]>::try_into(value).unwrap();
        Self::ntt(&mut array);
        Self::from_array(array)
    }
}

impl<const Q: u64, const N: usize> From<Zq<Q>> for Pow2CyclotomicPolyRingNTT<Q, N> {
    fn from(value: Zq<Q>) -> Self {
        Self::from_scalar(value)
    }
}

impl<const Q: u64, const N: usize> WithConjugationAutomorphism for Pow2CyclotomicPolyRingNTT<Q, N> {
    fn sigma(&self) -> Self {
        // TODO: can we implement the automorphism directly in NTT form?
        Into::<Pow2CyclotomicPolyRing<Zq<Q>, N>>::into(*self)
            .sigma()
            .into()
    }
}

impl<const Q: u64, const N: usize> WithL2Norm for Pow2CyclotomicPolyRingNTT<Q, N> {
    fn l2_norm_squared(&self) -> u128 {
        self.coeffs().l2_norm_squared()
    }
}

impl<const Q: u64, const N: usize> WithLinfNorm for Pow2CyclotomicPolyRingNTT<Q, N> {
    fn linf_norm(&self) -> u128 {
        self.coeffs().linf_norm()
    }
}

impl<const Q: u64, const N: usize> WithRot for Pow2CyclotomicPolyRingNTT<Q, N> {
    fn rot(&self) -> Matrix<Self::BaseRing> {
        let degree = Self::dimension();
        let coeffs = self.coeffs();
        let mut columns = Vec::with_capacity(degree);

        for i in 0..degree {
            let vec_xi_a = if i == 0 {
                Vector::from_vec(coeffs.clone())
            } else {
                Vector::from_vec(Self::multiply_by_xi(&coeffs, i))
            };
            columns.push(vec_xi_a);
        }

        Matrix::from_columns(columns.as_slice())
    }

    fn multiply_by_xi(bs: &Vec<Self::BaseRing>, i: usize) -> Vec<Self::BaseRing> {
        let len = bs.len();
        let mut result = vec![Self::BaseRing::ZERO; len];
        for (j, &coeff) in bs.iter().enumerate() {
            if j + i < len {
                result[(j + i) % len] += coeff;
            } else {
                result[(j + i) % len] -= coeff;
            }
        }
        result
    }
}
