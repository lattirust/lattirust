use std::fmt::Debug;
use std::hash::Hash;
use std::io::{Read, Write};
use std::iter::Product;
use std::ops::{Mul, MulAssign, Neg};

use ark_serialize::{
    CanonicalDeserialize, CanonicalSerialize, Compress, SerializationError, Valid, Validate,
};
use ark_std::rand::Rng;
use ark_std::UniformRand;
use derive_more::{Add, AddAssign, Display, From, Into, Sub, SubAssign, Sum};
use num_bigint::BigUint;
use num_traits::{One, Zero};

use crate::linear_algebra::SVector;
use crate::ring::ntt::NttRing;
use crate::ring::pow2_cyclotomic_poly_ring::Pow2CyclotomicPolyRing;
use crate::ring::PolyRing;
use crate::ring::Ring;
use crate::traits::{
    FromRandomBytes, Modulus, WithConjugationAutomorphism, WithL2Norm, WithLinfNorm,
};

#[derive(
    Clone,
    Copy,
    Debug,
    Eq,
    PartialEq,
    Hash,
    Add,
    AddAssign,
    Sum,
    Sub,
    SubAssign,
    From,
    Into,
    Display,
)]
#[display("NTT({})", self.0)]
pub struct Pow2CyclotomicPolyRingNTT<BaseRing: NttRing<N>, const N: usize>(SVector<BaseRing, N>);

impl<BaseRing: NttRing<N>, const N: usize> Pow2CyclotomicPolyRingNTT<BaseRing, N> {
    #[allow(dead_code)]
    pub(crate) type Inner = SVector<BaseRing, N>;

    pub type BaseRing = BaseRing;

    pub fn from_coefficient_array(mut coeffs: [BaseRing; N]) -> Self {
        Self::ntt(&mut coeffs);
        Self(Self::Inner::const_from_array(coeffs))
    }

    pub fn from_coefficients(coeffs: &[BaseRing]) -> Self {
        let mut evals_array: [BaseRing; N] = coeffs.try_into().expect(
            format!(
                "Invalid number of coefficients: expected {}, got {}",
                N,
                coeffs.len()
            )
            .as_str(),
        );
        Self::ntt(&mut evals_array);
        Self(Self::Inner::const_from_array(evals_array))
    }

    /// Constructs a polynomial from an array of coefficients in NTT form.
    /// This function is private since we can't enforce coeffs being in NTT form if called from outside.
    fn from_array(coeffs_ntt: [BaseRing; N]) -> Self {
        Self(Self::Inner::const_from_array(coeffs_ntt))
    }

    /// Constructs a polynomial from a function specifying coefficients in non-NTT form.
    pub fn from_fn<F>(f: F) -> Self
    where
        F: FnMut(usize) -> BaseRing,
    {
        let mut coeffs = core::array::from_fn(f);
        Self::ntt(&mut coeffs);
        Self::from_array(coeffs)
    }

    fn ntt(coeffs: &mut [BaseRing; N]) {
        BaseRing::ntt(coeffs);
    }

    fn intt(evals: &mut [BaseRing; N]) {
        BaseRing::intt(evals);
    }
    fn ntt_coeffs(&self) -> Vec<BaseRing> {
        self.0.iter().cloned().collect()
    }
}

impl<BaseRing: NttRing<N>, const N: usize> Modulus for Pow2CyclotomicPolyRingNTT<BaseRing, N> {
    fn modulus() -> BigUint {
        BaseRing::modulus()
    }
}

const fn vec_from_element<BaseRing: NttRing<N>, const N: usize>(
    elem: BaseRing,
) -> SVector<BaseRing, N> {
    SVector::<BaseRing, N>::const_from_array([elem; N])
}

impl<BaseRing: NttRing<N>, const N: usize> CanonicalSerialize
    for Pow2CyclotomicPolyRingNTT<BaseRing, N>
{
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

impl<BaseRing: NttRing<N>, const N: usize> Valid for Pow2CyclotomicPolyRingNTT<BaseRing, N> {
    fn check(&self) -> Result<(), SerializationError> {
        self.0.check()
    }
}

impl<BaseRing: NttRing<N>, const N: usize> CanonicalDeserialize
    for Pow2CyclotomicPolyRingNTT<BaseRing, N>
{
    fn deserialize_with_mode<R: Read>(
        reader: R,
        compress: Compress,
        validate: Validate,
    ) -> Result<Self, SerializationError> {
        Self::Inner::deserialize_with_mode(reader, compress, validate).map(Self)
    }
}

impl<BaseRing: NttRing<N>, const N: usize> Ring for Pow2CyclotomicPolyRingNTT<BaseRing, N> {
    const ZERO: Self = Self(vec_from_element(<BaseRing as Ring>::ZERO));
    const ONE: Self = Self(vec_from_element(<BaseRing as Ring>::ONE));

    fn inverse(&self) -> Option<Self> {
        todo!()
    }
}

impl<BaseRing: NttRing<N>, const N: usize> FromRandomBytes<Self>
    for Pow2CyclotomicPolyRingNTT<BaseRing, N>
{
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
        Some(Self::from_array(coeffs))
    }
}

impl<BaseRing: NttRing<N>, const N: usize> Default for Pow2CyclotomicPolyRingNTT<BaseRing, N> {
    #[inline(always)]
    fn default() -> Self {
        Self::zero()
    }
}

impl<BaseRing: NttRing<N>, const N: usize> Zero for Pow2CyclotomicPolyRingNTT<BaseRing, N> {
    #[inline(always)]
    fn zero() -> Self {
        Self::ZERO
    }

    #[inline(always)]
    fn is_zero(&self) -> bool {
        self.eq(&Self::ZERO)
    }
}

impl<BaseRing: NttRing<N>, const N: usize> One for Pow2CyclotomicPolyRingNTT<BaseRing, N> {
    #[inline(always)]
    fn one() -> Self {
        Self::ONE
    }
}

impl<BaseRing: NttRing<N>, const N: usize> Mul<Self> for Pow2CyclotomicPolyRingNTT<BaseRing, N> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        self.0.component_mul(&rhs.0).into()
    }
}

impl<BaseRing: NttRing<N>, const N: usize> Neg for Pow2CyclotomicPolyRingNTT<BaseRing, N> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        self.0.neg().into()
    }
}

impl<BaseRing: NttRing<N>, const N: usize> UniformRand for Pow2CyclotomicPolyRingNTT<BaseRing, N> {
    fn rand<R: Rng + ?Sized>(rng: &mut R) -> Self {
        Self::from_fn(|_| BaseRing::rand(rng))
    }
}

impl<BaseRing: NttRing<N>, const N: usize> MulAssign<Self>
    for Pow2CyclotomicPolyRingNTT<BaseRing, N>
{
    fn mul_assign(&mut self, rhs: Self) {
        self.0.component_mul_assign(&rhs.0)
    }
}

impl<'a, BaseRing: NttRing<N>, const N: usize> Add<&'a Self>
    for Pow2CyclotomicPolyRingNTT<BaseRing, N>
{
    type Output = Self;

    fn add(self, rhs: &'a Self) -> Self::Output {
        self.0.add(rhs.0).into()
    }
}

impl<'a, BaseRing: NttRing<N>, const N: usize> Sub<&'a Self>
    for Pow2CyclotomicPolyRingNTT<BaseRing, N>
{
    type Output = Self;

    fn sub(self, rhs: &'a Self) -> Self::Output {
        self.0.sub(rhs.0).into()
    }
}

impl<'a, BaseRing: NttRing<N>, const N: usize> Mul<&'a Self>
    for Pow2CyclotomicPolyRingNTT<BaseRing, N>
{
    type Output = Self;

    fn mul(self, rhs: &'a Self) -> Self::Output {
        self.0.component_mul(&rhs.0).into()
    }
}

impl<'a, BaseRing: NttRing<N>, const N: usize> AddAssign<&'a Self>
    for Pow2CyclotomicPolyRingNTT<BaseRing, N>
{
    fn add_assign(&mut self, rhs: &'a Self) {
        self.0.add_assign(rhs.0)
    }
}

impl<'a, BaseRing: NttRing<N>, const N: usize> SubAssign<&'a Self>
    for Pow2CyclotomicPolyRingNTT<BaseRing, N>
{
    fn sub_assign(&mut self, rhs: &'a Self) {
        self.0.sub_assign(rhs.0)
    }
}

impl<'a, BaseRing: NttRing<N>, const N: usize> MulAssign<&'a Self>
    for Pow2CyclotomicPolyRingNTT<BaseRing, N>
{
    fn mul_assign(&mut self, rhs: &'a Self) {
        self.0.component_mul_assign(&rhs.0)
    }
}

impl<'a, BaseRing: NttRing<N>, const N: usize> Add<&'a mut Self>
    for Pow2CyclotomicPolyRingNTT<BaseRing, N>
{
    type Output = Self;

    fn add(self, rhs: &'a mut Self) -> Self::Output {
        self.0.add(rhs.0).into()
    }
}

impl<'a, BaseRing: NttRing<N>, const N: usize> Sub<&'a mut Self>
    for Pow2CyclotomicPolyRingNTT<BaseRing, N>
{
    type Output = Self;

    fn sub(self, rhs: &'a mut Self) -> Self::Output {
        self.0.sub(rhs.0).into()
    }
}

impl<'a, BaseRing: NttRing<N>, const N: usize> Mul<&'a mut Self>
    for Pow2CyclotomicPolyRingNTT<BaseRing, N>
{
    type Output = Self;

    fn mul(self, rhs: &'a mut Self) -> Self::Output {
        self.0.component_mul(&rhs.0).into()
    }
}

impl<'a, BaseRing: NttRing<N>, const N: usize> AddAssign<&'a mut Self>
    for Pow2CyclotomicPolyRingNTT<BaseRing, N>
{
    fn add_assign(&mut self, rhs: &'a mut Self) {
        self.0.add_assign(rhs.0)
    }
}

impl<'a, BaseRing: NttRing<N>, const N: usize> SubAssign<&'a mut Self>
    for Pow2CyclotomicPolyRingNTT<BaseRing, N>
{
    fn sub_assign(&mut self, rhs: &'a mut Self) {
        self.0.sub_assign(rhs.0)
    }
}

impl<'a, BaseRing: NttRing<N>, const N: usize> MulAssign<&'a mut Self>
    for Pow2CyclotomicPolyRingNTT<BaseRing, N>
{
    fn mul_assign(&mut self, rhs: &'a mut Self) {
        self.0.component_mul_assign(&rhs.0)
    }
}

impl<BaseRing: NttRing<N>, const N: usize> From<Pow2CyclotomicPolyRing<BaseRing, N>>
    for Pow2CyclotomicPolyRingNTT<BaseRing, N>
{
    fn from(value: Pow2CyclotomicPolyRing<BaseRing, N>) -> Self {
        let mut coeffs: [BaseRing; N] = value.coeffs().try_into().unwrap();
        Self::ntt(&mut coeffs);
        Self::from_array(coeffs)
    }
}

impl<BaseRing: NttRing<N>, const N: usize> From<Pow2CyclotomicPolyRingNTT<BaseRing, N>>
    for Pow2CyclotomicPolyRing<BaseRing, N>
{
    fn from(val: Pow2CyclotomicPolyRingNTT<BaseRing, N>) -> Self {
        Pow2CyclotomicPolyRing::<BaseRing, N>::from(val.coeffs())
    }
}

impl<BaseRing: NttRing<N>, const N: usize> Mul<BaseRing>
    for Pow2CyclotomicPolyRingNTT<BaseRing, N>
{
    type Output = Self;

    fn mul(self, rhs: BaseRing) -> Self::Output {
        self.mul(Self::from_scalar(rhs))
    }
}

macro_rules! impl_try_from_primitive_type {
    ($primitive_type: ty) => {
        impl<BaseRing: NttRing<N>, const N: usize> TryFrom<$primitive_type>
            for Pow2CyclotomicPolyRingNTT<BaseRing, N>
        {
            type Error = <BaseRing as TryFrom<$primitive_type>>::Error;

            fn try_from(value: $primitive_type) -> Result<Self, Self::Error> {
                Ok(Self::from_scalar(BaseRing::try_from(value)?))
            }
        }
    };
}

macro_rules! impl_from_primitive_type {
    ($primitive_type: ty) => {
        impl<BaseRing: NttRing<N>, const N: usize> From<$primitive_type>
            for Pow2CyclotomicPolyRingNTT<BaseRing, N>
        {
            fn from(value: $primitive_type) -> Self {
                Self::from_scalar(BaseRing::from(value))
            }
        }
    };
}

impl_from_primitive_type!(bool);
impl_try_from_primitive_type!(u8);
impl_try_from_primitive_type!(u16);
impl_try_from_primitive_type!(u32);
impl_try_from_primitive_type!(u64);
impl_try_from_primitive_type!(u128);

impl<'a, BaseRing: NttRing<N>, const N: usize> Sum<&'a Self>
    for Pow2CyclotomicPolyRingNTT<BaseRing, N>
{
    fn sum<I: Iterator<Item = &'a Self>>(iter: I) -> Self {
        iter.fold(Self::zero(), |acc, x| acc + x)
    }
}

impl<BaseRing: NttRing<N>, const N: usize> Product<Self>
    for Pow2CyclotomicPolyRingNTT<BaseRing, N>
{
    fn product<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::one(), |acc, x| acc * x)
    }
}

impl<'a, BaseRing: NttRing<N>, const N: usize> Product<&'a Self>
    for Pow2CyclotomicPolyRingNTT<BaseRing, N>
{
    fn product<I: Iterator<Item = &'a Self>>(iter: I) -> Self {
        iter.fold(Self::one(), |acc, x| acc * x)
    }
}

impl<BaseRing: NttRing<N>, const N: usize> PolyRing for Pow2CyclotomicPolyRingNTT<BaseRing, N> {
    type BaseRing = BaseRing;
    fn coeffs(&self) -> Vec<Self::BaseRing> {
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

impl<BaseRing: NttRing<N>, const N: usize> From<Vec<BaseRing>>
    for Pow2CyclotomicPolyRingNTT<BaseRing, N>
{
    fn from(value: Vec<BaseRing>) -> Self {
        let n = value.len();
        let mut array = TryInto::<[BaseRing; N]>::try_into(value).expect(
            format!(
                "Invalid vector length {} for polynomial ring dimension N={}",
                n, N
            )
            .as_str(),
        );
        Self::ntt(&mut array);
        Self::from_array(array)
    }
}

impl<BaseRing: NttRing<N>, const N: usize> From<BaseRing>
    for Pow2CyclotomicPolyRingNTT<BaseRing, N>
{
    fn from(value: BaseRing) -> Self {
        Self::from_scalar(value)
    }
}

impl<BaseRing: NttRing<N>, const N: usize> WithConjugationAutomorphism
    for Pow2CyclotomicPolyRingNTT<BaseRing, N>
{
    fn sigma(&self) -> Self {
        // TODO: can we implement the automorphism directly in NTT form?
        Into::<Pow2CyclotomicPolyRing<BaseRing, N>>::into(*self)
            .sigma()
            .into()
    }
}

impl<BaseRing: NttRing<N>, const N: usize> WithL2Norm for Pow2CyclotomicPolyRingNTT<BaseRing, N>
where
    Vec<BaseRing>: WithL2Norm,
{
    fn l2_norm_squared(&self) -> BigUint {
        self.coeffs().l2_norm_squared()
    }
}

impl<BaseRing: NttRing<N>, const N: usize> WithLinfNorm for Pow2CyclotomicPolyRingNTT<BaseRing, N>
where
    Vec<BaseRing>: WithLinfNorm,
{
    fn linf_norm(&self) -> BigUint {
        self.coeffs().linf_norm()
    }
}
