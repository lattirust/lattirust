use ark_ff::{Field, Fp, FpConfig};
use ark_serialize::{
    CanonicalDeserialize, CanonicalSerialize, Compress, SerializationError, Valid, Validate,
};
use ark_std::{
    array,
    convert::TryInto,
    fmt::{Debug, Display, Formatter},
    hash::Hash,
    io::{Read, Write},
    iter::{Iterator, Product, Sum},
    ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign},
    rand::Rng,
    One, UniformRand, Zero,
};
use derive_more::{From, Into};

use super::{ring_config::CyclotomicConfig, CyclotomicPolyRingNTTGeneral};
use crate::{
    balanced_decomposition::{decompose_balanced_polyring, Decompose},
    traits::FromRandomBytes,
    PolyRing, Ring,
};

/// A cyclotomic ring Fp[X]/(Phi_m(X)) in the coefficient form.
/// * `C` is the configuration of the cyclotomic ring.
/// * `N` is the byte size of the underlying prime field.
/// * `D` is the degree of the cyclotomic polynomial.
#[derive(From, Into)]
pub struct CyclotomicPolyRingGeneral<C: CyclotomicConfig<N>, const N: usize, const D: usize>(
    pub(crate) [Fp<C::BaseFieldConfig, N>; D],
);

impl<C: CyclotomicConfig<N>, const N: usize, const D: usize> CyclotomicPolyRingGeneral<C, N, D> {
    fn from_coeffs_vec(mut coeffs: Vec<Fp<C::BaseFieldConfig, N>>) -> Self {
        C::reduce_in_place(&mut coeffs);

        Self(coeffs.try_into().expect("Wrong length"))
    }

    fn from_fn<F>(f: F) -> Self
    where
        F: FnMut(usize) -> Fp<C::BaseFieldConfig, N>,
    {
        let coeffs = core::array::from_fn::<_, D, _>(f);
        Self::from_array(coeffs)
    }

    fn from_array(coeffs: [Fp<C::BaseFieldConfig, N>; D]) -> Self {
        Self(coeffs)
    }

    fn poly_mul(&self, rhs: &Self) -> Self {
        let lhs_coeffs = self.coeffs();
        let rhs_coeffs = rhs.coeffs();
        let mut coeffs = vec![<Fp::<C::BaseFieldConfig, N> as Field>::ZERO; 2 * D - 1];
        for i in 0..D {
            for j in 0..D {
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
impl<C: CyclotomicConfig<N>, const N: usize, const D: usize> PartialEq
    for CyclotomicPolyRingGeneral<C, N, D>
{
    fn eq(&self, other: &Self) -> bool {
        self.0 == other.0
    }
}

impl<C: CyclotomicConfig<N>, const N: usize, const D: usize> Eq
    for CyclotomicPolyRingGeneral<C, N, D>
{
}

impl<C: CyclotomicConfig<N>, const N: usize, const D: usize> Clone
    for CyclotomicPolyRingGeneral<C, N, D>
{
    fn clone(&self) -> Self {
        *self
    }
}

impl<C: CyclotomicConfig<N>, const N: usize, const D: usize> Copy
    for CyclotomicPolyRingGeneral<C, N, D>
{
}

impl<C: CyclotomicConfig<N>, const N: usize, const D: usize> Debug
    for CyclotomicPolyRingGeneral<C, N, D>
{
    fn fmt(&self, f: &mut Formatter<'_>) -> ark_std::fmt::Result {
        Debug::fmt(&self.0, f)
    }
}

impl<C: CyclotomicConfig<N>, const N: usize, const D: usize> Display
    for CyclotomicPolyRingGeneral<C, N, D>
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

impl<C: CyclotomicConfig<N>, const N: usize, const D: usize> Hash
    for CyclotomicPolyRingGeneral<C, N, D>
{
    fn hash<H: ark_std::hash::Hasher>(&self, state: &mut H) {
        self.0.hash(state);
    }
}

const fn array_one<FP: FpConfig<N>, const N: usize, const D: usize>() -> [Fp<FP, N>; D] {
    let mut array = [FP::ZERO; D];
    array[0] = FP::ONE;
    array
}

impl<C: CyclotomicConfig<N>, const N: usize, const D: usize> CanonicalSerialize
    for CyclotomicPolyRingGeneral<C, N, D>
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

impl<C: CyclotomicConfig<N>, const N: usize, const D: usize> Valid
    for CyclotomicPolyRingGeneral<C, N, D>
{
    fn check(&self) -> Result<(), SerializationError> {
        self.0.check()
    }
}

impl<C: CyclotomicConfig<N>, const N: usize, const D: usize> CanonicalDeserialize
    for CyclotomicPolyRingGeneral<C, N, D>
{
    fn deserialize_with_mode<R: Read>(
        reader: R,
        compress: Compress,
        validate: Validate,
    ) -> Result<Self, SerializationError> {
        <[Fp<C::BaseFieldConfig, N>; D]>::deserialize_with_mode(reader, compress, validate)
            .map(Self)
    }
}

impl<C: CyclotomicConfig<N>, const N: usize, const D: usize> Ring
    for CyclotomicPolyRingGeneral<C, N, D>
{
    const ZERO: Self = Self([<Fp<C::BaseFieldConfig, N> as Field>::ZERO; D]);
    const ONE: Self = Self(array_one());
}

impl<C: CyclotomicConfig<N>, const N: usize, const D: usize> FromRandomBytes<Self>
    for CyclotomicPolyRingGeneral<C, N, D>
{
    fn byte_size() -> usize {
        D * Fp::<C::BaseFieldConfig, N>::byte_size()
    }

    fn try_from_random_bytes(bytes: &[u8]) -> Option<Self> {
        assert_eq!(bytes.len(), Self::byte_size());
        let coeffs = core::array::from_fn(|i| {
            Fp::<C::BaseFieldConfig, N>::try_from_random_bytes(
                &bytes[i * Fp::<C::BaseFieldConfig, N>::byte_size()
                    ..(i + 1) * Fp::<C::BaseFieldConfig, N>::byte_size()],
            )
            .unwrap()
        });
        Some(Self::from_array(coeffs))
    }
}

impl<C: CyclotomicConfig<N>, const N: usize, const D: usize> Default
    for CyclotomicPolyRingGeneral<C, N, D>
{
    #[inline(always)]
    fn default() -> Self {
        Self::zero()
    }
}

impl<C: CyclotomicConfig<N>, const N: usize, const D: usize> Zero
    for CyclotomicPolyRingGeneral<C, N, D>
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

impl<C: CyclotomicConfig<N>, const N: usize, const D: usize> One
    for CyclotomicPolyRingGeneral<C, N, D>
{
    #[inline(always)]
    fn one() -> Self {
        Self::ONE
    }
}

impl<C: CyclotomicConfig<N>, const N: usize, const D: usize> Mul<Self>
    for CyclotomicPolyRingGeneral<C, N, D>
{
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        self.poly_mul(&rhs)
    }
}

impl<C: CyclotomicConfig<N>, const N: usize, const D: usize> Neg
    for CyclotomicPolyRingGeneral<C, N, D>
{
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self(self.0.map(|x| -x))
    }
}

impl<C: CyclotomicConfig<N>, const N: usize, const D: usize> UniformRand
    for CyclotomicPolyRingGeneral<C, N, D>
{
    fn rand<R: Rng + ?Sized>(rng: &mut R) -> Self {
        Self::from_fn(|_| Fp::<C::BaseFieldConfig, N>::rand(rng))
    }
}

impl<C: CyclotomicConfig<N>, const N: usize, const D: usize> MulAssign<Self>
    for CyclotomicPolyRingGeneral<C, N, D>
{
    fn mul_assign(&mut self, rhs: Self) {
        self.poly_mul_in_place(&rhs);
    }
}

impl<'a, C: CyclotomicConfig<N>, const N: usize, const D: usize> Add<&'a Self>
    for CyclotomicPolyRingGeneral<C, N, D>
{
    type Output = Self;

    fn add(self, rhs: &'a Self) -> Self::Output {
        Self(array::from_fn(|i| self.0[i] + rhs.0[i]))
    }
}

impl<'a, C: CyclotomicConfig<N>, const N: usize, const D: usize> Sub<&'a Self>
    for CyclotomicPolyRingGeneral<C, N, D>
{
    type Output = Self;

    fn sub(self, rhs: &'a Self) -> Self::Output {
        Self(array::from_fn(|i| self.0[i] - rhs.0[i]))
    }
}

impl<'a, C: CyclotomicConfig<N>, const N: usize, const D: usize> Add<&'a mut Self>
    for CyclotomicPolyRingGeneral<C, N, D>
{
    type Output = Self;

    fn add(self, rhs: &'a mut Self) -> Self::Output {
        for (a, b) in rhs.0.iter_mut().zip(self.0.iter()) {
            *a += *b;
        }
        *rhs
    }
}

impl<'a, C: CyclotomicConfig<N>, const N: usize, const D: usize> Sub<&'a mut Self>
    for CyclotomicPolyRingGeneral<C, N, D>
{
    type Output = Self;

    fn sub(self, rhs: &'a mut Self) -> Self::Output {
        for (a, b) in rhs.0.iter_mut().zip(self.0.iter()) {
            *a -= *b;
        }
        *rhs
    }
}

impl<'a, C: CyclotomicConfig<N>, const N: usize, const D: usize> Mul<&'a mut Self>
    for CyclotomicPolyRingGeneral<C, N, D>
{
    type Output = Self;

    fn mul(self, rhs: &'a mut Self) -> Self::Output {
        self.poly_mul(rhs)
    }
}

impl<'a, C: CyclotomicConfig<N>, const N: usize, const D: usize> AddAssign<&'a mut Self>
    for CyclotomicPolyRingGeneral<C, N, D>
{
    fn add_assign(&mut self, rhs: &'a mut Self) {
        for (a, b) in self.0.iter_mut().zip(rhs.0.iter()) {
            *a += *b;
        }
    }
}

impl<'a, C: CyclotomicConfig<N>, const N: usize, const D: usize> SubAssign<&'a mut Self>
    for CyclotomicPolyRingGeneral<C, N, D>
{
    fn sub_assign(&mut self, rhs: &'a mut Self) {
        for (a, b) in self.0.iter_mut().zip(rhs.0.iter()) {
            *a -= *b;
        }
    }
}

impl<'a, C: CyclotomicConfig<N>, const N: usize, const D: usize> MulAssign<&'a mut Self>
    for CyclotomicPolyRingGeneral<C, N, D>
{
    fn mul_assign(&mut self, rhs: &'a mut Self) {
        self.poly_mul_in_place(rhs);
    }
}

impl<C: CyclotomicConfig<N>, const N: usize, const D: usize>
    From<CyclotomicPolyRingNTTGeneral<C, N, D>>
    for CyclotomicPolyRingGeneral<C, N, { D * C::CRT_FIELD_EXTENSION_DEGREE }>
{
    fn from(value: CyclotomicPolyRingNTTGeneral<C, N, D>) -> Self {
        let coeffs: Vec<Fp<C::BaseFieldConfig, N>> = C::icrt(value.into_coeffs());

        Self(coeffs.try_into().unwrap())
    }
}

impl<C: CyclotomicConfig<N>, const N: usize, const D: usize> Mul<Fp<C::BaseFieldConfig, N>>
    for CyclotomicPolyRingGeneral<C, N, D>
{
    type Output = Self;

    fn mul(self, rhs: Fp<C::BaseFieldConfig, N>) -> Self::Output {
        self.poly_mul(&Self::from_scalar(rhs))
    }
}

macro_rules! impl_from_primitive_type {
    ($primitive_type: ty) => {
        impl<C: CyclotomicConfig<N>, const N: usize, const D: usize> From<$primitive_type>
            for CyclotomicPolyRingGeneral<C, N, D>
        {
            fn from(value: $primitive_type) -> Self {
                Self::from_scalar(Fp::<C::BaseFieldConfig, N>::from(value))
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

impl<'a, C: CyclotomicConfig<N>, const N: usize, const D: usize> Mul<&'a Self>
    for CyclotomicPolyRingGeneral<C, N, D>
{
    type Output = Self;

    fn mul(self, rhs: &'a Self) -> Self::Output {
        self.poly_mul(rhs)
    }
}

impl<'a, C: CyclotomicConfig<N>, const N: usize, const D: usize> AddAssign<&'a Self>
    for CyclotomicPolyRingGeneral<C, N, D>
{
    fn add_assign(&mut self, rhs: &'a Self) {
        for (a, b) in self.0.iter_mut().zip(rhs.0.iter()) {
            *a += *b;
        }
    }
}

impl<'a, C: CyclotomicConfig<N>, const N: usize, const D: usize> SubAssign<&'a Self>
    for CyclotomicPolyRingGeneral<C, N, D>
{
    fn sub_assign(&mut self, rhs: &'a Self) {
        for (a, b) in self.0.iter_mut().zip(rhs.0.iter()) {
            *a -= *b;
        }
    }
}

impl<'a, C: CyclotomicConfig<N>, const N: usize, const D: usize> MulAssign<&'a Self>
    for CyclotomicPolyRingGeneral<C, N, D>
{
    fn mul_assign(&mut self, rhs: &'a Self) {
        self.poly_mul_in_place(rhs)
    }
}

impl<C: CyclotomicConfig<N>, const N: usize, const D: usize> Add<Self>
    for CyclotomicPolyRingGeneral<C, N, D>
{
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self(array::from_fn(|i| self.0[i] + rhs.0[i]))
    }
}

impl<C: CyclotomicConfig<N>, const N: usize, const D: usize> Sub<Self>
    for CyclotomicPolyRingGeneral<C, N, D>
{
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self(array::from_fn(|i| self.0[i] - rhs.0[i]))
    }
}

impl<C: CyclotomicConfig<N>, const N: usize, const D: usize> AddAssign<Self>
    for CyclotomicPolyRingGeneral<C, N, D>
{
    fn add_assign(&mut self, rhs: Self) {
        for (a, b) in self.0.iter_mut().zip(rhs.0.iter()) {
            *a += *b;
        }
    }
}

impl<C: CyclotomicConfig<N>, const N: usize, const D: usize> SubAssign<Self>
    for CyclotomicPolyRingGeneral<C, N, D>
{
    fn sub_assign(&mut self, rhs: Self) {
        for (a, b) in self.0.iter_mut().zip(rhs.0.iter()) {
            *a -= *b;
        }
    }
}

impl<C: CyclotomicConfig<N>, const N: usize, const D: usize> Sum<Self>
    for CyclotomicPolyRingGeneral<C, N, D>
{
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::zero(), |acc, x| acc + x)
    }
}

impl<'a, C: CyclotomicConfig<N>, const N: usize, const D: usize> Sum<&'a Self>
    for CyclotomicPolyRingGeneral<C, N, D>
{
    fn sum<I: Iterator<Item = &'a Self>>(iter: I) -> Self {
        iter.fold(Self::zero(), |acc, x| acc + x)
    }
}

impl<C: CyclotomicConfig<N>, const N: usize, const D: usize> Product<Self>
    for CyclotomicPolyRingGeneral<C, N, D>
{
    fn product<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::one(), |acc, x| acc * x)
    }
}

impl<'a, C: CyclotomicConfig<N>, const N: usize, const D: usize> Product<&'a Self>
    for CyclotomicPolyRingGeneral<C, N, D>
{
    fn product<I: Iterator<Item = &'a Self>>(iter: I) -> Self {
        iter.fold(Self::one(), |acc, x| acc * x)
    }
}

impl<C: CyclotomicConfig<N>, const N: usize, const D: usize> PolyRing
    for CyclotomicPolyRingGeneral<C, N, D>
{
    type BaseRing = Fp<C::BaseFieldConfig, N>;

    fn coeffs(&self) -> &[Fp<C::BaseFieldConfig, N>] {
        &self.0
    }

    fn coeffs_mut(&mut self) -> &mut [Self::BaseRing] {
        &mut self.0
    }

    fn dimension() -> usize {
        D
    }

    fn from_scalar(v: Self::BaseRing) -> Self {
        // NTT([v, 0, ..., 0]) = ([v, ..., v])
        let mut coeffs = [Self::BaseRing::zero(); D];
        coeffs[0] = v;
        Self::from_array(coeffs)
    }

    fn into_coeffs(self) -> Vec<Self::BaseRing> {
        self.0.into()
    }
}

impl<C: CyclotomicConfig<N>, const N: usize, const D: usize> From<Vec<Fp<C::BaseFieldConfig, N>>>
    for CyclotomicPolyRingGeneral<C, N, D>
{
    fn from(mut value: Vec<Fp<C::BaseFieldConfig, N>>) -> Self {
        if value.len() < D {
            value.resize(D, Fp::zero())
        }

        Self::from_coeffs_vec(value)
    }
}

impl<C: CyclotomicConfig<N>, const N: usize, const D: usize> From<Fp<C::BaseFieldConfig, N>>
    for CyclotomicPolyRingGeneral<C, N, D>
{
    fn from(value: Fp<C::BaseFieldConfig, N>) -> Self {
        Self::from_scalar(value)
    }
}

impl<C: CyclotomicConfig<N>, const N: usize, const D: usize> Decompose
    for CyclotomicPolyRingGeneral<C, N, D>
where
    Fp<C::BaseFieldConfig, N>: Decompose,
{
    fn decompose(&self, b: u128, padding_size: Option<usize>) -> Vec<Self> {
        decompose_balanced_polyring(self, b, padding_size)
    }
}

macro_rules! impl_add_mul_primitive_type {
    ($primitive_type: ty) => {
        impl<C: CyclotomicConfig<N>, const N: usize, const D: usize> Mul<$primitive_type>
            for CyclotomicPolyRingGeneral<C, N, D>
        {
            type Output = Self;

            fn mul(mut self, rhs: $primitive_type) -> Self::Output {
                self.0.iter_mut().for_each(|lhs| *lhs *= Fp::from(rhs));

                self
            }
        }

        impl<C: CyclotomicConfig<N>, const N: usize, const D: usize> Add<$primitive_type>
            for CyclotomicPolyRingGeneral<C, N, D>
        {
            type Output = Self;

            fn add(mut self, rhs: $primitive_type) -> Self::Output {
                self.0[0] += Fp::from(rhs);

                self
            }
        }
    };
}

// only works for types that can be cast to Field
impl_add_mul_primitive_type!(u128);
impl_add_mul_primitive_type!(u64);
impl_add_mul_primitive_type!(u32);
impl_add_mul_primitive_type!(u16);
impl_add_mul_primitive_type!(u8);
impl_add_mul_primitive_type!(bool);

#[cfg(test)]
mod tests {
    use crate::{
        cyclotomic_ring::models::pow2_debug::{Pow2CyclotomicPolyRing, Pow2CyclotomicPolyRingNTT},
        Ring,
    };

    #[test]
    fn test_ntt_one() {
        let one = Pow2CyclotomicPolyRing::<131072, 1024>::ONE;

        assert_eq!(
            Pow2CyclotomicPolyRingNTT::from(one),
            Pow2CyclotomicPolyRingNTT::<131072, 1024>::ONE
        )
    }

    #[test]
    fn test_intt_one() {
        let one = Pow2CyclotomicPolyRingNTT::<131072, 1024>::ONE;

        assert_eq!(
            Pow2CyclotomicPolyRing::<131072, 1024>::from(one),
            Pow2CyclotomicPolyRing::<131072, 1024>::ONE
        )
    }
}
