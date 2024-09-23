use ark_ff::{Field, Fp, FpConfig};
use ark_std::fmt::{Debug, Display, Formatter};
use ark_std::hash::Hash;
use ark_std::io::{Read, Write};
use ark_std::iter::{Product, Sum};
use ark_std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

use ark_serialize::{
    CanonicalDeserialize, CanonicalSerialize, Compress, SerializationError, Valid, Validate,
};
use ark_std::rand::Rng;
use ark_std::UniformRand;
use derive_more::{From, Into};
use num_traits::{One, Zero};

use crate::balanced_decomposition::{decompose_balanced_polyring, Decompose};
use crate::pow2_cyclotomic_poly_ring_ntt::Fp64Pow2;
use crate::ring_config::{Pow2Rp64Config, RpConfig};
use crate::traits::{FromRandomBytes, WithL2Norm, WithLinfNorm};
use crate::Ring;
use crate::{CyclotomicPolyRingNTTGeneral, PolyRing};
use lattirust_linear_algebra::SVector;

use super::poly_ring::WithRot;

#[derive(From, Into)]
pub struct CyclotomicPolyRingGeneral<C: RpConfig<N>, const N: usize, const PHI_D: usize>(
    SVector<Fp<C::FpConfig, N>, PHI_D>,
);

pub type Pow2CyclotomicPolyRing<const Q: u64, const PHI_D: usize> =
    CyclotomicPolyRingGeneral<Pow2Rp64Config<Q, PHI_D>, 1, PHI_D>;

impl<C: RpConfig<N> + 'static, const N: usize, const PHI_D: usize>
    CyclotomicPolyRingGeneral<C, N, PHI_D>
{
    pub fn from_coeffs_vec(coeffs: Vec<Fp<C::FpConfig, N>>) -> Self {
        let mut coeffs = coeffs;
        C::reduce_in_place(&mut coeffs);
        Self::from_array(coeffs.try_into().map_err(|v: Vec<_>| v).unwrap())
    }
    pub fn from_fn<F>(f: F) -> Self
    where
        F: FnMut(usize) -> Fp<C::FpConfig, N>,
    {
        let coeffs = core::array::from_fn::<_, PHI_D, _>(f);
        Self::from_array(coeffs)
    }
    pub fn from_array(coeffs: [Fp<C::FpConfig, N>; PHI_D]) -> Self {
        Self(SVector::const_from_array(coeffs))
    }

    fn poly_mul(&self, rhs: &Self) -> Self {
        let lhs_coeffs = self.coeffs();
        let rhs_coeffs = rhs.coeffs();
        let mut coeffs = vec![<Fp::<C::FpConfig, N> as Field>::ZERO; 2 * PHI_D - 1];
        for i in 0..PHI_D {
            for j in 0..PHI_D {
                coeffs[i + j] += lhs_coeffs[i] * rhs_coeffs[j];
            }
        }
        C::reduce_in_place(&mut coeffs);
        Self::from_coeffs_vec(coeffs)
    }
    fn poly_mul_in_place(&mut self, rhs: &Self) {
        // we need a reduce function for SVector to properly do a multiplication in place
        let res = *self * rhs;
        self.0 = res.0;
    }
    // TODO: test mul in place
}
impl<C: RpConfig<N>, const N: usize, const PHI_D: usize> PartialEq
    for CyclotomicPolyRingGeneral<C, N, PHI_D>
{
    fn eq(&self, other: &Self) -> bool {
        self.0 == other.0
    }
}

impl<C: RpConfig<N>, const N: usize, const PHI_D: usize> Eq
    for CyclotomicPolyRingGeneral<C, N, PHI_D>
{
}

impl<C: RpConfig<N>, const N: usize, const PHI_D: usize> Clone
    for CyclotomicPolyRingGeneral<C, N, PHI_D>
{
    fn clone(&self) -> Self {
        Self(self.0.clone())
    }
}

impl<C: RpConfig<N>, const N: usize, const PHI_D: usize> Copy
    for CyclotomicPolyRingGeneral<C, N, PHI_D>
{
}

impl<C: RpConfig<N>, const N: usize, const PHI_D: usize> Debug
    for CyclotomicPolyRingGeneral<C, N, PHI_D>
{
    fn fmt(&self, f: &mut Formatter<'_>) -> ark_std::fmt::Result {
        Debug::fmt(&self.0, f)
    }
}

impl<C: RpConfig<N>, const N: usize, const PHI_D: usize> Display
    for CyclotomicPolyRingGeneral<C, N, PHI_D>
{
    fn fmt(&self, f: &mut Formatter<'_>) -> ark_std::fmt::Result {
        write!(f, "CyclotomicPolyRingGeneral(")?;
        let mut iter = self.0.iter();
        if let Some(first) = iter.next() {
            write!(f, "{}", first)?;
            for field_element in iter {
                write!(f, ", {}", field_element)?;
            }
        }
        write!(f, ")")
    }
}

impl<C: RpConfig<N>, const N: usize, const PHI_D: usize> Hash
    for CyclotomicPolyRingGeneral<C, N, PHI_D>
{
    fn hash<H: ark_std::hash::Hasher>(&self, state: &mut H) {
        self.0.hash(state);
    }
}

const fn vec_from_element<FP: FpConfig<N>, const N: usize, const PHI_D: usize>(
    elem: Fp<FP, N>,
) -> SVector<Fp<FP, N>, PHI_D> {
    SVector::const_from_array([elem; PHI_D])
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
//
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
//
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
//
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

impl<C: RpConfig<N>, const N: usize, const PHI_D: usize> CanonicalSerialize
    for CyclotomicPolyRingGeneral<C, N, PHI_D>
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

impl<C: RpConfig<N>, const N: usize, const PHI_D: usize> Valid
    for CyclotomicPolyRingGeneral<C, N, PHI_D>
{
    fn check(&self) -> Result<(), SerializationError> {
        self.0.check()
    }
}

impl<C: RpConfig<N>, const N: usize, const PHI_D: usize> CanonicalDeserialize
    for CyclotomicPolyRingGeneral<C, N, PHI_D>
{
    fn deserialize_with_mode<R: Read>(
        reader: R,
        compress: Compress,
        validate: Validate,
    ) -> Result<Self, SerializationError> {
        SVector::<Fp<C::FpConfig, N>, PHI_D>::deserialize_with_mode(reader, compress, validate)
            .map(Self)
    }
}

impl<C: RpConfig<N> + 'static, const N: usize, const PHI_D: usize> Ring
    for CyclotomicPolyRingGeneral<C, N, PHI_D>
{
    const ZERO: Self = Self(vec_from_element(<Fp<C::FpConfig, N> as Field>::ZERO));
    const ONE: Self = Self(vec_from_element(<Fp<C::FpConfig, N> as Field>::ONE));
}

impl<C: RpConfig<N> + 'static, const N: usize, const PHI_D: usize> FromRandomBytes<Self>
    for CyclotomicPolyRingGeneral<C, N, PHI_D>
{
    fn byte_size() -> usize {
        PHI_D * Fp::<C::FpConfig, N>::byte_size()
    }

    fn try_from_random_bytes(bytes: &[u8]) -> Option<Self> {
        assert_eq!(bytes.len(), Self::byte_size());
        let coeffs = core::array::from_fn(|i| {
            Fp::<C::FpConfig, N>::try_from_random_bytes(
                &bytes[i * Fp::<C::FpConfig, N>::byte_size()
                    ..(i + 1) * Fp::<C::FpConfig, N>::byte_size()],
            )
            .unwrap()
        });
        Some(Self::from_array(coeffs))
    }
}

impl<C: RpConfig<N> + 'static, const N: usize, const PHI_D: usize> Default
    for CyclotomicPolyRingGeneral<C, N, PHI_D>
{
    #[inline(always)]
    fn default() -> Self {
        Self::zero()
    }
}

impl<C: RpConfig<N> + 'static, const N: usize, const PHI_D: usize> Zero
    for CyclotomicPolyRingGeneral<C, N, PHI_D>
{
    #[inline(always)]
    fn zero() -> Self {
        Self::ZERO
    }

    #[inline(always)]
    fn is_zero(&self) -> bool {
        self.eq(&Self::ZERO)
    }
}

impl<C: RpConfig<N> + 'static, const N: usize, const PHI_D: usize> One
    for CyclotomicPolyRingGeneral<C, N, PHI_D>
{
    #[inline(always)]
    fn one() -> Self {
        Self::ONE
    }
}

impl<C: RpConfig<N> + 'static, const N: usize, const PHI_D: usize> Mul<Self>
    for CyclotomicPolyRingGeneral<C, N, PHI_D>
{
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        self.poly_mul(&rhs)
    }
}

impl<C: RpConfig<N>, const N: usize, const PHI_D: usize> Neg
    for CyclotomicPolyRingGeneral<C, N, PHI_D>
{
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self(self.0.neg())
    }
}

impl<C: RpConfig<N> + 'static, const N: usize, const PHI_D: usize> UniformRand
    for CyclotomicPolyRingGeneral<C, N, PHI_D>
{
    fn rand<R: Rng + ?Sized>(rng: &mut R) -> Self {
        Self::from_fn(|_| Fp::<C::FpConfig, N>::rand(rng))
    }
}

impl<C: RpConfig<N> + 'static, const N: usize, const PHI_D: usize> MulAssign<Self>
    for CyclotomicPolyRingGeneral<C, N, PHI_D>
{
    fn mul_assign(&mut self, rhs: Self) {
        self.poly_mul_in_place(&rhs);
    }
}

impl<'a, C: RpConfig<N>, const N: usize, const PHI_D: usize> Add<&'a Self>
    for CyclotomicPolyRingGeneral<C, N, PHI_D>
{
    type Output = Self;

    fn add(self, rhs: &'a Self) -> Self::Output {
        Self(self.0 + rhs.0)
    }
}

impl<'a, C: RpConfig<N>, const N: usize, const PHI_D: usize> Sub<&'a Self>
    for CyclotomicPolyRingGeneral<C, N, PHI_D>
{
    type Output = Self;

    fn sub(self, rhs: &'a Self) -> Self::Output {
        Self(self.0 - rhs.0)
    }
}

impl<'a, C: RpConfig<N>, const N: usize, const PHI_D: usize> Add<&'a mut Self>
    for CyclotomicPolyRingGeneral<C, N, PHI_D>
{
    type Output = Self;

    fn add(self, rhs: &'a mut Self) -> Self::Output {
        Self(self.0 + rhs.0)
    }
}

impl<'a, C: RpConfig<N>, const N: usize, const PHI_D: usize> Sub<&'a mut Self>
    for CyclotomicPolyRingGeneral<C, N, PHI_D>
{
    type Output = Self;

    fn sub(self, rhs: &'a mut Self) -> Self::Output {
        Self(self.0 - rhs.0)
    }
}

impl<'a, C: RpConfig<N> + 'static, const N: usize, const PHI_D: usize> Mul<&'a mut Self>
    for CyclotomicPolyRingGeneral<C, N, PHI_D>
{
    type Output = Self;

    fn mul(self, rhs: &'a mut Self) -> Self::Output {
        self.poly_mul(rhs)
    }
}

impl<'a, C: RpConfig<N>, const N: usize, const PHI_D: usize> AddAssign<&'a mut Self>
    for CyclotomicPolyRingGeneral<C, N, PHI_D>
{
    fn add_assign(&mut self, rhs: &'a mut Self) {
        self.0 += rhs.0;
    }
}

impl<'a, C: RpConfig<N>, const N: usize, const PHI_D: usize> SubAssign<&'a mut Self>
    for CyclotomicPolyRingGeneral<C, N, PHI_D>
{
    fn sub_assign(&mut self, rhs: &'a mut Self) {
        self.0 -= rhs.0;
    }
}

impl<'a, C: RpConfig<N> + 'static, const N: usize, const PHI_D: usize> MulAssign<&'a mut Self>
    for CyclotomicPolyRingGeneral<C, N, PHI_D>
{
    fn mul_assign(&mut self, rhs: &'a mut Self) {
        self.poly_mul_in_place(&rhs);
    }
}

impl<C: RpConfig<N> + 'static, const N: usize, const PHI_D: usize>
    From<CyclotomicPolyRingNTTGeneral<C, N, PHI_D>> for CyclotomicPolyRingGeneral<C, N, PHI_D>
{
    fn from(value: CyclotomicPolyRingNTTGeneral<C, N, PHI_D>) -> Self {
        let coeffs: Vec<Fp<C::FpConfig, N>> = value.coeffs();
        Self(
            coeffs
                .try_into()
                .map_err(|v: Vec<Fp<C::FpConfig, N>>| v)
                .unwrap(),
        )
    }
}

impl<C: RpConfig<N> + 'static, const N: usize, const PHI_D: usize> Mul<Fp<C::FpConfig, N>>
    for CyclotomicPolyRingGeneral<C, N, PHI_D>
{
    type Output = Self;

    fn mul(self, rhs: Fp<C::FpConfig, N>) -> Self::Output {
        self.poly_mul(&Self::from_scalar(rhs))
    }
}

macro_rules! impl_from_primitive_type {
    ($primitive_type: ty) => {
        impl<C: RpConfig<N> + 'static, const N: usize, const PHI_D: usize> From<$primitive_type>
            for CyclotomicPolyRingGeneral<C, N, PHI_D>
        {
            fn from(value: $primitive_type) -> Self {
                Self::from_scalar(Fp::<C::FpConfig, N>::from(value))
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

impl<'a, C: RpConfig<N> + 'static, const N: usize, const PHI_D: usize> Mul<&'a Self>
    for CyclotomicPolyRingGeneral<C, N, PHI_D>
{
    type Output = Self;

    fn mul(self, rhs: &'a Self) -> Self::Output {
        self.poly_mul(rhs)
    }
}

impl<'a, C: RpConfig<N>, const N: usize, const PHI_D: usize> AddAssign<&'a Self>
    for CyclotomicPolyRingGeneral<C, N, PHI_D>
{
    fn add_assign(&mut self, rhs: &'a Self) {
        self.0 += rhs.0;
    }
}

impl<'a, C: RpConfig<N>, const N: usize, const PHI_D: usize> SubAssign<&'a Self>
    for CyclotomicPolyRingGeneral<C, N, PHI_D>
{
    fn sub_assign(&mut self, rhs: &'a Self) {
        self.0 -= rhs.0;
    }
}

impl<'a, C: RpConfig<N> + 'static, const N: usize, const PHI_D: usize> MulAssign<&'a Self>
    for CyclotomicPolyRingGeneral<C, N, PHI_D>
{
    fn mul_assign(&mut self, rhs: &'a Self) {
        self.poly_mul_in_place(rhs)
    }
}

impl<'a, C: RpConfig<N>, const N: usize, const PHI_D: usize> Add<Self>
    for CyclotomicPolyRingGeneral<C, N, PHI_D>
{
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self(self.0 + rhs.0)
    }
}

impl<'a, C: RpConfig<N>, const N: usize, const PHI_D: usize> Sub<Self>
    for CyclotomicPolyRingGeneral<C, N, PHI_D>
{
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self(self.0 - rhs.0)
    }
}

impl<'a, C: RpConfig<N>, const N: usize, const PHI_D: usize> AddAssign<Self>
    for CyclotomicPolyRingGeneral<C, N, PHI_D>
{
    fn add_assign(&mut self, rhs: Self) {
        self.0 += rhs.0
    }
}

impl<'a, C: RpConfig<N>, const N: usize, const PHI_D: usize> SubAssign<Self>
    for CyclotomicPolyRingGeneral<C, N, PHI_D>
{
    fn sub_assign(&mut self, rhs: Self) {
        self.0 -= rhs.0
    }
}

impl<'a, C: RpConfig<N> + 'static, const N: usize, const PHI_D: usize> Sum<Self>
    for CyclotomicPolyRingGeneral<C, N, PHI_D>
{
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::zero(), |acc, x| acc + x)
    }
}

impl<'a, C: RpConfig<N> + 'static, const N: usize, const PHI_D: usize> Sum<&'a Self>
    for CyclotomicPolyRingGeneral<C, N, PHI_D>
{
    fn sum<I: Iterator<Item = &'a Self>>(iter: I) -> Self {
        iter.fold(Self::zero(), |acc, x| acc + x)
    }
}

impl<C: RpConfig<N> + 'static, const N: usize, const PHI_D: usize> Product<Self>
    for CyclotomicPolyRingGeneral<C, N, PHI_D>
{
    fn product<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::one(), |acc, x| acc * x)
    }
}

impl<'a, C: RpConfig<N> + 'static, const N: usize, const PHI_D: usize> Product<&'a Self>
    for CyclotomicPolyRingGeneral<C, N, PHI_D>
{
    fn product<I: Iterator<Item = &'a Self>>(iter: I) -> Self {
        iter.fold(Self::one(), |acc, x| acc * x)
    }
}

impl<C: RpConfig<N> + 'static, const N: usize, const PHI_D: usize> PolyRing
    for CyclotomicPolyRingGeneral<C, N, PHI_D>
{
    type BaseRing = Fp<C::FpConfig, N>;
    fn coeffs(&self) -> Vec<Fp<C::FpConfig, N>> {
        self.0.as_slice().to_vec()
    }
    fn dimension() -> usize {
        PHI_D
    }

    fn from_scalar(v: Self::BaseRing) -> Self {
        // NTT([v, 0, ..., 0]) = ([v, ..., v])
        Self::from_array([v; PHI_D])
    }
}

impl<C: RpConfig<N>, const N: usize, const PHI_D: usize> From<Vec<Fp<C::FpConfig, N>>>
    for CyclotomicPolyRingGeneral<C, N, PHI_D>
{
    fn from(mut value: Vec<Fp<C::FpConfig, N>>) -> Self {
        value.resize_with(PHI_D, Fp::zero);
        C::crt_in_place(&mut value);
        Self(
            value
                .try_into()
                .map_err(|v: Vec<Fp<C::FpConfig, N>>| v)
                .unwrap(),
        )
    }
}

impl<C: RpConfig<N> + 'static, const N: usize, const PHI_D: usize> From<Fp<C::FpConfig, N>>
    for CyclotomicPolyRingGeneral<C, N, PHI_D>
{
    fn from(value: Fp<C::FpConfig, N>) -> Self {
        Self::from_scalar(value)
    }
}

// impl<BaseRing: Ring, const N: usize> WithConjugationAutomorphism
//     for Pow2CyclotomicPolyRing<BaseRing, N>
// {
//     fn sigma(&self) -> Self {
//         let coeffs = self.0 .0.as_slice();
//         let mut new_coeffs = Vec::<BaseRing>::with_capacity(N);
//         new_coeffs.push(coeffs[0]);
//         new_coeffs.extend(
//             coeffs[1..]
//                 .iter()
//                 .rev()
//                 .map(|v_i| -*v_i)
//                 .collect::<Vec<BaseRing>>(),
//         );
//         Self::from(new_coeffs)
//     }
// }

impl<const Q: u64, const PHI_D: usize> WithRot for Pow2CyclotomicPolyRing<Q, PHI_D> {
    fn multiply_by_xi(&self, i: usize) -> Vec<Self::BaseRing> {
        let bs = self.0;
        let len = bs.ncols();
        assert_eq!(len, PHI_D);
        let mut result = vec![<Fp64Pow2<Q, PHI_D> as Field>::ZERO; len];
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

impl<const Q: u64, const PHI_D: usize> Decompose for Pow2CyclotomicPolyRing<Q, PHI_D> {
    fn decompose(&self, b: u128, padding_size: Option<usize>) -> Vec<Self> {
        decompose_balanced_polyring(self, b, padding_size)
    }
}
impl<const Q: u64, const PHI_D: usize> WithL2Norm for Pow2CyclotomicPolyRing<Q, PHI_D> {
    fn l2_norm_squared(&self) -> u128 {
        self.coeffs().l2_norm_squared()
    }
}

impl<const Q: u64, const PHI_D: usize> WithLinfNorm for Pow2CyclotomicPolyRing<Q, PHI_D> {
    fn linf_norm(&self) -> u128 {
        self.coeffs().linf_norm()
    }
}
