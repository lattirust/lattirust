use std::hash::Hash;

use ark_std::UniformRand;
use std::fmt::Display;
use derive_more::{From, Into, Debug};
use crate::ring::Fq;
use crate::ring::z_q::{Zq1, Zq2, Zq3, Zq4, Zq5};

macro_rules! pow2_cyclotomic_poly_ring_ntt_crt {
    ($L:literal $(,$l:literal)+) => {
            paste::expr! {
#[derive(
    Clone,
    Copy,
    Debug,
    Eq,
    PartialEq,
    Hash,
    From,
    Into,
    Debug,
)]
pub struct [< Pow2CyclotomicPolyRingNTT $L >]<N, $([< Q $l >],)*>(($([Fq<[< Q $l >]>; N],)*));
    
// impl<const N: usize, Q1: u64, ..., QL: u64> Pow2CyclotomicPolyRingNTTL<const N: usize, Q1: u64, ..., QL: u64> {
impl<const N: usize, $(const [< Q $l >]: u64,)*> [< Pow2CyclotomicPolyRingNTT $L >]<N, $([< Q $l >],)*>
    where $(Fq::<[< Q $l >]>: Ntt,)*
{
    #[allow(dead_code)]
    pub(crate) type Inner = ($([Fq<[< Q $l >]>; N],)*); // ([Fq<Q1>; N], ..., [Fq<QL>; N])

    pub type BaseRing = [< Zq $L >]<$([< Q $l >],)*>;
    
    /// Convert from array-of-structs to struct-of-arrays
    fn aos_to_soa(aos: [Self::BaseRing; N]) -> Inner {
        ($(core::array::from_fn(|i| coeffs[i].$l),)*)
    }
    
    /// Convert from struct-of-arrays to array-of-structs
    fn soa_to_aos(soa: Inner) -> [Self::BaseRing; N] {
        core::array::from_fn(|i| ($(soa.$l[i],)*))
    }

    pub fn from_coefficient_array(coeffs: &[Self::BaseRing; N]) -> Self {
        let mut struct_of_arrays = Self::aos_to_soa(coeffs);
        Self::ntt_inplace(&mut struct_of_arrays);
        Self(struct_of_arrays)
    }

    /// Constructs a polynomial from the values of the polynomial in NTT form.
    pub fn from_ntt_array(coeffs_ntt: [Self::BaseRing; N]) -> Self {
        Self(Self::aos_to_soa(coeffs_ntt))
    }

    /// Constructs a polynomial from a function specifying coefficients in non-NTT form.
    pub fn from_fn<F>(f: F) -> Self
    where
        F: FnMut(usize) -> Self::BaseRing,
    {
        Self::from_coefficient_array(core::array::from_fn(f))
    }
    
    fn ntt_inplace(coeffs: &mut Self::Inner) {
        $(
            Fq::<[< Q $l >]>::ntt_inplace(&mut coeffs.$l);
        )*
    }
    
    fn intt_inplace(evals: &mut Self::Inner) {
        $(
            Fq::<[< Q $l >]>::intt_inplace(&mut evals.$l);
        )*
    }

    pub fn ntt_values(&self) -> [Self::BaseRing; N] {
        Self::soa_to_aos(self.0)
    }
}
    
impl<const N: usize, $(const [< Q $l >]: u64,)*> Display for [< Pow2CyclotomicPolyRingNTT $L >]<N, $([< Q $l >],)*>
    where Self::BaseRing: Display    
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "NTT({})", Self::soa_to_aos(self.0))
    }
}
    
impl<const N: usize, $(const [< Q $l >]: u64,)*> Modulus for [< Pow2CyclotomicPolyRingNTT $L >]<N, $([< Q $l >],)*>
{
    fn modulus() -> BigUint {
        BaseRing::modulus()
    }
}

impl<const N: usize, $(const [< Q $l >]: u64,)*> CanonicalSerialize for [< Pow2CyclotomicPolyRingNTT $L >]<N, $([< Q $l >],)*>
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

impl<const N: usize, $(const [< Q $l >]: u64,)*> Valid for [< Pow2CyclotomicPolyRingNTT $L >]<N, $([< Q $l >],)*>
{
    fn check(&self) -> Result<(), SerializationError> {
        self.0.check()
    }
}

impl<const N: usize, $(const [< Q $l >]: u64,)*> CanonicalDeserialize for [< Pow2CyclotomicPolyRingNTT $L >]<N, $([< Q $l >],)*>
{
    fn deserialize_with_mode<R: Read>(
        reader: R,
        compress: Compress,
        validate: Validate,
    ) -> Result<Self, SerializationError> {
        Self::Inner::deserialize_with_mode(reader, compress, validate).map(Self)
    }
}

impl<const N: usize, $(const [< Q $l >]: u64,)*> Ring for [< Pow2CyclotomicPolyRingNTT $L >]<N, $([< Q $l >],)*>
{
    const ZERO: Self = Self($([Fq::<[< Q $l >]>::ZERO; N],)*);
    const ONE: Self = Self($([Fq::<[< Q $l >]>::ONE; N],)*);

    fn inverse(&self) -> Option<Self> {
        let mut inv = Self::zero();
        $(
            match self.0.$l.try_map(|x| x.inverse()) {
                Some(v) => inv.0.$l = v,
                None => return None,
            }
        )*
    }
}

impl<const N: usize, $(const [< Q $l >]: u64,)*> FromRandomBytes for [< Pow2CyclotomicPolyRingNTT $L >]<N, $([< Q $l >],)*>
    where Self::Inner: FromRandomBytes,
{
    fn needs_bytes() -> usize {
        Self::Inner::needs_bytes()
    //N * Self::BaseRing::byte_size()
    }

    fn try_from_random_bytes_inner(bytes: &[u8]) -> Option<Self> {
        Self::Inner::try_from_random_bytes_inner(bytes)
    }
}

impl<const N: usize, $(const [< Q $l >]: u64,)*> Default for [< Pow2CyclotomicPolyRingNTT $L >]<N, $([< Q $l >],)*>
{
    #[inline(always)]
    fn default() -> Self {
        Self::zero()
    }
}

impl<const N: usize, $(const [< Q $l >]: u64,)*> Zero for [< Pow2CyclotomicPolyRingNTT $L >]<N, $([< Q $l >],)*>
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

impl<const N: usize, $(const [< Q $l >]: u64,)*> One for [< Pow2CyclotomicPolyRingNTT $L >]<N, $([< Q $l >],)*>
{
    #[inline(always)]
    fn one() -> Self {
        Self::ONE
    }
}

impl<const N: usize, $(const [< Q $l >]: u64,)*> Mul<Self> for [< Pow2CyclotomicPolyRingNTT $L >]<N, $([< Q $l >],)*>
{
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        Self($(
            core::array::from_fn(|i| self.0.$l[i].mul(rhs.0.$l[i]))
        ,)*)
    }
}

impl<const N: usize, $(const [< Q $l >]: u64,)*> Neg for [< Pow2CyclotomicPolyRingNTT $L >]<N, $([< Q $l >],)*>
{
    type Output = Self;

    fn neg(self) -> Self::Output {
        self.0.neg().into()
    }
}

impl<const N: usize, $(const [< Q $l >]: u64,)*> UniformRand for [< Pow2CyclotomicPolyRingNTT $L >]<N, $([< Q $l >],)*>
{
    fn rand<R: Rng + ?Sized>(rng: &mut R) -> Self {
        // Self::from_fn(|_| BaseRing::rand(rng))
        Pow2CyclotomicPolyRing::rand(rng).into()
    }
}

impl<const N: usize, $(const [< Q $l >]: u64,)*> MulAssign<Self> for [< Pow2CyclotomicPolyRingNTT $L >]<N, $([< Q $l >],)*>
{
    fn mul_assign(&mut self, rhs: Self) {
        $(
        self.0.$l.as_mut_slice().zip(rhs.0.$l.iter()).for_each(|(a, b)| { a.mul_assign(b); });
        )*
    }
}

impl<'a, const N: usize, $(const [< Q $l >]: u64,)*> Add<&'a Self> for [< Pow2CyclotomicPolyRingNTT $L >]<N, $([< Q $l >],)*>
{
    type Output = Self;

    fn add(self, rhs: &'a Self) -> Self::Output {
        Self($(
            core::array::from_fn(|i| self.0.$l[i].add(rhs.0.$l[i]))
        ,)*)
    }
}

impl<'a, const N: usize, $(const [< Q $l >]: u64,)*> Sub<&'a Self> for [< Pow2CyclotomicPolyRingNTT $L >]<N, $([< Q $l >],)*>
{
    type Output = Self;

    fn sub(self, rhs: &'a Self) -> Self::Output {
        Self($(
            core::array::from_fn(|i| self.0.$l[i].sub(rhs.0.$l[i]))
        ,)*)
    }
}

impl<'a, const N: usize, $(const [< Q $l >]: u64,)*> Mul<&'a Self> for [< Pow2CyclotomicPolyRingNTT $L >]<N, $([< Q $l >],)*>
{
    type Output = Self;

    fn mul(self, rhs: &'a Self) -> Self::Output {
        Self($(
           core::array::from_fn(|i| self.0.$l[i].mul(rhs.0.$l[i]))
        ,)*)
    }
}

impl<'a, const N: usize, $(const [< Q $l >]: u64,)*> AddAssign<&'a Self> for [< Pow2CyclotomicPolyRingNTT $L >]<N, $([< Q $l >],)*>
{
    fn add_assign(&mut self, rhs: &'a Self) {
       $(
        self.0.$l.as_mut_slice().zip(rhs.0.$l.iter()).for_each(|(a, b)| { a.add_assign(b); });
        )*
    }
}

impl<'a, const N: usize, $(const [< Q $l >]: u64,)*> SubAssign<&'a Self> for [< Pow2CyclotomicPolyRingNTT $L >]<N, $([< Q $l >],)*>
{
    fn sub_assign(&mut self, rhs: &'a Self) {
        $(
        self.0.$l.as_mut_slice().zip(rhs.0.$l.iter()).for_each(|(a, b)| { a.sub_assign(b); });
        )*
    }
}

impl<'a, const N: usize, $(const [< Q $l >]: u64,)*> MulAssign<&'a Self> for [< Pow2CyclotomicPolyRingNTT $L >]<N, $([< Q $l >],)*>
{
    fn mul_assign(&mut self, rhs: &'a Self) {
        $(
        self.0.$l.as_mut_slice().zip(rhs.0.$l.iter()).for_each(|(a, b)| { a.mul_assign(b); });
        )*
    }
}

// impl<'a, BaseRing: NttRing<N>, const N: usize> Add<&'a mut Self>
//     for Pow2CyclotomicPolyRingNTT<BaseRing, N>
// {
//     type Output = Self;
//
//     fn add(self, rhs: &'a mut Self) -> Self::Output {
//         self.0.add(rhs.0).into()
//     }
// }
//
// impl<'a, BaseRing: NttRing<N>, const N: usize> Sub<&'a mut Self>
//     for Pow2CyclotomicPolyRingNTT<BaseRing, N>
// {
//     type Output = Self;
//
//     fn sub(self, rhs: &'a mut Self) -> Self::Output {
//         self.0.sub(rhs.0).into()
//     }
// }

// impl<'a, BaseRing: NttRing<N>, const N: usize> Mul<&'a mut Self>
//     for Pow2CyclotomicPolyRingNTT<BaseRing, N>
// {
//     type Output = Self;
//
//     fn mul(self, rhs: &'a mut Self) -> Self::Output {
//         self.0.component_mul(&rhs.0).into()
//     }
// }

// impl<'a, BaseRing: NttRing<N>, const N: usize> AddAssign<&'a mut Self>
//     for Pow2CyclotomicPolyRingNTT<BaseRing, N>
// {
//     fn add_assign(&mut self, rhs: &'a mut Self) {
//         self.0.add_assign(rhs.0)
//     }
// }
//
// impl<'a, BaseRing: NttRing<N>, const N: usize> SubAssign<&'a mut Self>
//     for Pow2CyclotomicPolyRingNTT<BaseRing, N>
// {
//     fn sub_assign(&mut self, rhs: &'a mut Self) {
//         self.0.sub_assign(rhs.0)
//     }
// }
//
// impl<'a, BaseRing: NttRing<N>, const N: usize> MulAssign<&'a mut Self>
//     for Pow2CyclotomicPolyRingNTT<BaseRing, N>
// {
//     fn mul_assign(&mut self, rhs: &'a mut Self) {
//         self.0.component_mul_assign(&rhs.0)
//     }
// }

impl<const N: usize, $(const [< Q $l >]: u64,)*> From<Pow2CyclotomicPolyRing<Self::BaseRing, N>> for [< Pow2CyclotomicPolyRingNTT $L >]<N, $([< Q $l >],)*>
{
    fn from(value: Pow2CyclotomicPolyRing<Self::BaseRing, N>) -> Self {
        Self::from_coefficient_array(value.coefficient_array())
    }
}

impl<const N: usize, $(const [< Q $l >]: u64,)*> From<[< Pow2CyclotomicPolyRingNTT $L >]<N, $([< Q $l >],)*>> for Pow2CyclotomicPolyRing<BaseRing, N>
{
    fn from(val: [< Pow2CyclotomicPolyRingNTT $L >]<N, $([< Q $l >],)*>) -> Self {
        Self::from(val.coefficients())
    }
}

impl<const N: usize, $(const [< Q $l >]: u64,)*> Mul<Self::BaseRing> for [< Pow2CyclotomicPolyRingNTT $L >]<N, $([< Q $l >],)*>
{
    type Output = Self;

    fn mul(self, rhs: BaseRing) -> Self::Output {
        self.mul(Self::from_scalar(rhs))
    }
}

macro_rules! impl_try_from_primitive_type {
    ($primitive_type: ty) => {
        impl<const N: usize, $(const [< Q $l >]: u64,)*> TryFrom<$primitive_type> for [< Pow2CyclotomicPolyRingNTT $L >]<N, $([< Q $l >],)*>
        {
            type Error = <Self::BaseRing as TryFrom<$primitive_type>>::Error;

            fn try_from(value: $primitive_type) -> Result<Self, Self::Error> {
                Ok(Self::from_scalar(Self::BaseRing::try_from(value)?))
            }
        }
    }
}


macro_rules! impl_from_primitive_type {
    ($primitive_type: ty) => {
        impl<const N: usize, $(const [< Q $l >]: u64,)*> From<$primitive_type> for [< Pow2CyclotomicPolyRingNTT $L >]<N, $([< Q $l >],)*>
        {
            fn from(value: $primitive_type) -> Self {
                Self::from_scalar(BaseRing::from(value))
            }
        }
    }
}

impl_from_primitive_type!(bool);
impl_try_from_primitive_type!(u8);
impl_try_from_primitive_type!(u16);
impl_try_from_primitive_type!(u32);
impl_try_from_primitive_type!(u64);
impl_try_from_primitive_type!(u128);

impl<const N: usize, $(const [< Q $l >]: u64,)*> Sum<Self> for [< Pow2CyclotomicPolyRingNTT $L >]<N, $([< Q $l >],)*>
{
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::zero(), |acc, x| acc + x)
    }
}

impl<'a, const N: usize, $(const [< Q $l >]: u64,)*> Sum<&'a Self> for [< Pow2CyclotomicPolyRingNTT $L >]<N, $([< Q $l >],)*>
{
    fn sum<I: Iterator<Item = &'a Self>>(iter: I) -> Self {
        iter.fold(Self::zero(), |acc, x| acc + x)
    }
}

impl<const N: usize, $(const [< Q $l >]: u64,)*> Product<Self> for [< Pow2CyclotomicPolyRingNTT $L >]<N, $([< Q $l >],)*>
{
    fn product<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::one(), |acc, x| acc * x)
    }
}

impl<'a, const N: usize, $(const [< Q $l >]: u64,)*> Product<&'a Self> for [< Pow2CyclotomicPolyRingNTT $L >]<N, $([< Q $l >],)*>
{
    fn product<I: Iterator<Item = &'a Self>>(iter: I) -> Self {
        iter.fold(Self::one(), |acc, x| acc * x)
    }
}

impl<const N: usize, $(const [< Q $l >]: u64,)*> PolyRing for [< Pow2CyclotomicPolyRingNTT $L >]<N, $([< Q $l >],)*>
{
    type BaseRing = BaseRing;

    /// Return the coefficients of the polynomial in non-NTT form.
    fn coefficients(&self) -> Vec<Self::BaseRing> {
        let mut coeffs = self.0 .0.into();
        Self::intt_inplace(&mut coeffs);
        coeffs.to_vec()
    }

    fn try_from_coefficients(coeffs: &[Self::BaseRing]) -> Option<Self> {
        let arr: [BaseRing; N] = coeffs.try_into().ok()?;
        Some(Self::from_coefficient_array(arr))
    }

    fn dimension() -> usize {
        N
    }

    fn from_scalar(v: Self::BaseRing) -> Self {
        // NTT([v, 0, ..., 0]) = ([v, ..., v])
        Self::from_ntt_array([v; N])
    }
}

impl<const N: usize, $(const [< Q $l >]: u64,)*> From<Vec<Self::BaseRing>> for [< Pow2CyclotomicPolyRingNTT $L >]<N, $([< Q $l >],)*>
{
    /// Construct a polynomial from a vector of coefficients in non-NTT form.
    fn from(value: Vec<Self::BaseRing>) -> Self {
        let n = value.len();
        let array = TryInto::<[Self::BaseRing; N]>::try_into(value).unwrap_or_else(|_| {
            panic!(
                "Invalid vector length {} for polynomial ring dimension N={}",
                n, N
            )
        });
        Self::from_coefficient_array(array)
    }
}

impl<const N: usize, $(const [< Q $l >]: u64,)*> From<Self::BaseRing> for [< Pow2CyclotomicPolyRingNTT $L >]<N, $([< Q $l >],)*>
{
    fn from(value: BaseRing) -> Self {
        Self::from_scalar(value)
    }
}

impl<const N: usize, $(const [< Q $l >]: u64,)*> WithConjugationAutomorphism for [< Pow2CyclotomicPolyRingNTT $L >]<N, $([< Q $l >],)*>
{
    fn apply_automorphism(&self) -> Self {
        // TODO: can we implement the automorphism directly in NTT form?
        Into::<Pow2CyclotomicPolyRing<BaseRing, N>>::into(*self)
            .apply_automorphism()
            .into()
    }
}

impl<const N: usize, $(const [< Q $l >]: u64,)*> WithL2Norm for [< Pow2CyclotomicPolyRingNTT $L >]<N, $([< Q $l >],)*>
where
    Vec<Self::BaseRing>: WithL2Norm,
{
    fn l2_norm_squared(&self) -> BigUint {
        self.coefficients().l2_norm_squared()
    }
}

impl<const N: usize, $(const [< Q $l >]: u64,)*> WithLinfNorm for [< Pow2CyclotomicPolyRingNTT $L >]<N, $([< Q $l >],)*>
where
    Vec<Self::BaseRing>: WithLinfNorm,
{
    fn linf_norm(&self) -> BigUint {
        self.coefficients().linf_norm()
    }
}

}
}
}

pow2_cyclotomic_poly_ring_ntt_crt!(1, 1);
pow2_cyclotomic_poly_ring_ntt_crt!(2, 1, 2);
pow2_cyclotomic_poly_ring_ntt_crt!(3, 1, 2, 3);
pow2_cyclotomic_poly_ring_ntt_crt!(4, 1, 2, 3, 4);
pow2_cyclotomic_poly_ring_ntt_crt!(5, 1, 2, 3, 4, 5);
pow2_cyclotomic_poly_ring_ntt_crt!(6, 1, 2, 3, 4, 5, 6);
pow2_cyclotomic_poly_ring_ntt_crt!(7, 1, 2, 3, 4, 5, 6, 7);
pow2_cyclotomic_poly_ring_ntt_crt!(8, 1, 2, 3, 4, 5, 6, 7, 8);
pow2_cyclotomic_poly_ring_ntt_crt!(9, 1, 2, 3, 4, 5, 6, 7, 8, 9);
pow2_cyclotomic_poly_ring_ntt_crt!(10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10);
    
#[cfg(test)]
mod test {
    use crate::*;
    use super::*;

    const N: usize = 64;
    const Q: u64 = 65537;
    type PR = Pow2CyclotomicPolyRingNTT1<Q, N>;
    type BR = PR::BaseRing;
    const NUM_TEST_REPETITIONS: usize = 100;

    test_ring!(PR, NUM_TEST_REPETITIONS);

    test_polyring!(PR, NUM_TEST_REPETITIONS);

    // #[test]
    // fn test_from_into_pow2cyclotomicpolyring() {
    //     let rng = &mut ark_std::test_rng();
    //     for _ in 0..NUM_TEST_REPETITIONS {
    //         let p_ntt = <PR as UniformRand>::rand(rng);
    //         let p: Pow2CyclotomicPolyRing<BR, N> = p_ntt.into();
    //         let p_ntt_: Pow2CyclotomicPolyRingNTT<BR, N> = p.into();
    //         assert_eq!(p_ntt, p_ntt_);
    //
    //         let p = Pow2CyclotomicPolyRing::<BR, N>::rand(rng);
    //         let p_ntt: Pow2CyclotomicPolyRingNTT<BR, N> = p.into();
    //         let p_: Pow2CyclotomicPolyRing<BR, N> = p_ntt.into();
    //         assert_eq!(p, p_);
    //     }
    // }

    // #[test]
    // fn test_ntt_add() {
    //     let rng = &mut ark_std::test_rng();
    //     for _ in 0..NUM_TEST_REPETITIONS {
    //         let a = Pow2CyclotomicPolyRing::<BR, N>::rand(rng);
    //         let b = Pow2CyclotomicPolyRing::<BR, N>::rand(rng);
    //         let a_ntt: PR = a.into();
    //         let b_ntt: PR = b.into();
    //
    //         let a_plus_b = a + b;
    //         let a_plus_b_ntt = a_ntt + b_ntt;
    //         assert_eq!(a_plus_b_ntt.coefficients(), a_plus_b.coefficients());
    //
    //         let a_plus_b_: Pow2CyclotomicPolyRing<BR, N> = a_plus_b_ntt.into();
    //         assert_eq!(a_plus_b_, a_plus_b);
    //     }
    // }
    //
    // #[test]
    // fn test_ntt_mul() {
    //     use crate::ring::ntt::Ntt;
    //
    //     let rng = &mut ark_std::test_rng();
    //
    //     let mut a: [BR; N] = core::array::from_fn(|_| BR::rand(rng));
    //     let mut b: [BR; N] = core::array::from_fn(|_| BR::rand(rng));
    //
    //     let mut a_mul_b_naive: [BR; N] = core::array::from_fn(|_| BR::zero());
    //     for i in 0..N {
    //         for j in 0..N {
    //             if i + j < N {
    //                 a_mul_b_naive[i + j] += a[i] * b[j];
    //             } else {
    //                 a_mul_b_naive[i + j - N] -= a[i] * b[j];
    //             }
    //         }
    //     }
    //
    //     BR::ntt_inplace&mut a);
    //     BR::ntt_inplace&mut b);
    //     let mut a_mul_b_ntt: [BR; N] = core::array::from_fn(|i| a[i] * b[i]);
    //     BR::intt_inplace&mut a_mul_b_ntt);
    //
    //     assert_eq!(a_mul_b_ntt, a_mul_b_naive);
    //
    //     for _ in 0..NUM_TEST_REPETITIONS {
    //         let a_arr: [BR; N] = core::array::from_fn(|_| BR::rand(rng));
    //         let b_arr: [BR; N] = core::array::from_fn(|_| BR::rand(rng));
    //
    //         let a = Pow2CyclotomicPolyRing::from(a_arr.clone()); //Pow2CyclotomicPolyRing::<BR, N>::rand(rng);
    //         let b = Pow2CyclotomicPolyRing::from(b_arr.clone()); // Pow2CyclotomicPolyRing::<BR, N>::rand(rng);
    //         assert_eq!(a.coefficients(), a_arr.to_vec());
    //         assert_eq!(b.coefficients(), b_arr.to_vec());
    //
    //         //
    //         let a_ntt: PR = a.clone().into();
    //         let b_ntt: PR = b.clone().into();
    //         assert_eq!(
    //             a.coefficients(),
    //             a_ntt.coefficients(),
    //             "Mismatch in a → a_ntt"
    //         );
    //         assert_eq!(
    //             b.coefficients(),
    //             b_ntt.coefficients(),
    //             "Mismatch in b → b_ntt"
    //         );
    //
    //         let mut a_mul_b_arr_naive: [BR; N] = core::array::from_fn(|_| BR::zero());
    //         for i in 0..N {
    //             for j in 0..N {
    //                 if i + j < N {
    //                     a_mul_b_arr_naive[i + j] += a_arr[i] * b_arr[j];
    //                 } else {
    //                     a_mul_b_arr_naive[i + j - N] -= a_arr[i] * b_arr[j];
    //                 }
    //             }
    //         }
    //
    //         //
    //         let mut a_arr_ntt = a_arr.clone();
    //         let mut b_arr_ntt = b_arr.clone();
    //         BR::ntt_inplace&mut a_arr_ntt);
    //         BR::ntt_inplace&mut b_arr_ntt);
    //         assert_eq!(a_ntt.ntt_values(), a_arr_ntt.to_vec());
    //         assert_eq!(b_ntt.ntt_values(), b_arr_ntt.to_vec());
    //
    //         //
    //         let a_mul_b_arr_ntt: [BR; N] = core::array::from_fn(|i| a_arr_ntt[i] * b_arr_ntt[i]);
    //         let mut a_mul_b_arr = a_mul_b_arr_ntt.clone();
    //         BR::intt_inplace&mut a_mul_b_arr);
    //         assert_eq!(a_mul_b_arr, a_mul_b_arr_naive);
    //
    //         let a_mul_b = a * b;
    //         let a_mul_b_ntt = a_ntt * b_ntt;
    //
    //         assert_eq!(
    //             a_mul_b_ntt.ntt_values(),
    //             a_mul_b_arr_ntt.to_vec(),
    //             "PolyNTT and NTT should match"
    //         );
    //         assert_eq!(
    //             a_mul_b.coefficients(),
    //             a_mul_b_arr_naive.to_vec(),
    //             "Poly and naive should match"
    //         );
    //         assert_eq!(
    //             a_mul_b.coefficients(),
    //             a_mul_b_arr.to_vec(),
    //             "Poly and INTT(NTT) should match"
    //         );
    //         assert_eq!(
    //             a_mul_b_ntt.coefficients(),
    //             a_mul_b_arr.to_vec(),
    //             "INTT(PolyNTT) and INTT(NTT) should match"
    //         );
    //         assert_eq!(
    //             a_mul_b_ntt.coefficients(),
    //             a_mul_b.coefficients(),
    //             "INTT(PolyNTT) and Poly should match"
    //         );
    //
    //         //
    //         let a_mul_b_: Pow2CyclotomicPolyRing<BR, N> = a_mul_b_ntt.into();
    //         assert_eq!(a_mul_b_, a_mul_b);
    //     }
    // }
    //
    // test_conjugation_automorphism!(PR, NUM_TEST_REPETITIONS);
    //
    // #[test]
    // fn test_conjugation_automorphism_inner_product_debug() {
    //     use crate::linear_algebra::Vector;
    //
    //     let rng = &mut ark_std::test_rng();
    //     let a = <PR as UniformRand>::rand(rng);
    //     let b = <PR as UniformRand>::rand(rng);
    //
    //     let a_coeffs = Vector::<<PR as PolyRing>::BaseRing>::from_vec(a.coefficients());
    //     let b_coeffs = Vector::<<PR as PolyRing>::BaseRing>::from_vec(b.coefficients());
    //
    //     let a_non_ntt = Into::<Pow2CyclotomicPolyRing<BR, N>>::into(a.clone());
    //     let b_non_ntt = Into::<Pow2CyclotomicPolyRing<BR, N>>::into(b.clone());
    //
    //     // Check that <a_coeffs, b_coeffs> = ct(sigma_{-1}(a) * b)
    //     let a_sigma = a.apply_automorphism();
    //     let a_nonntt_sigma = a_non_ntt.apply_automorphism();
    //     let a_sigma_nonntt = Into::<Pow2CyclotomicPolyRing<BR, N>>::into(a_sigma.clone());
    //     assert_eq!(a_nonntt_sigma, a_sigma_nonntt);
    //
    //     assert_eq!(
    //         Into::<Pow2CyclotomicPolyRing<BR, N>>::into(a_sigma.clone()),
    //         a_nonntt_sigma
    //     );
    //     assert_eq!(
    //         Into::<Pow2CyclotomicPolyRingNTT<BR, N>>::into(a_nonntt_sigma.clone()),
    //         a_sigma
    //     );
    //
    //     let a_b = a_sigma * b;
    //     let a_b_nonntt = a_nonntt_sigma * b_non_ntt;
    //     assert_eq!(
    //         Into::<Pow2CyclotomicPolyRingNTT<BR, N>>::into(a_nonntt_sigma.clone()),
    //         a_sigma
    //     );
    //     assert_eq!(
    //         Into::<Pow2CyclotomicPolyRingNTT<BR, N>>::into(b_non_ntt.clone()),
    //         b
    //     );
    //     assert_eq!(
    //         Into::<Pow2CyclotomicPolyRingNTT<BR, N>>::into(a_b_nonntt.clone()),
    //         a_b
    //     );
    //     assert_eq!(
    //         Into::<Pow2CyclotomicPolyRing<BR, N>>::into(a_b.clone()),
    //         a_b_nonntt
    //     );
    //
    //     let inner_prod = a_coeffs.dot(&b_coeffs);
    //     assert_eq!(inner_prod, a_b_nonntt.coefficients()[0]);
    //     assert_eq!(inner_prod, a_b.coefficients()[0]);
    // }
}
