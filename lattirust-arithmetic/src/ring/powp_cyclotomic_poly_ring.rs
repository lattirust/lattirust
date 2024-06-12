use std::fmt::{Debug, Display, Formatter};
use std::ops::{Add, Mul};

use ark_std::rand::Rng;
use ark_std::UniformRand;
use num_bigint::BigUint;
use num_traits::{One, Zero};

use crate::linear_algebra::SVector;
use crate::ring::Ring;
use crate::ring::{ConvertibleRing, PolyRing};

use crate::traits::{FromRandomBytes, Modulus, WithConjugationAutomorphism, WithL2Norm, WithLinfNorm};

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct PowPCyclotomicPolyRing<BaseRing: ConvertibleRing, const P: usize, const N: usize>(
    SVector<BaseRing, N>,
);

impl<BaseRing: ConvertibleRing, const P: usize, const N: usize>
    PowPCyclotomicPolyRing<BaseRing, P, N>
{
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

impl<BaseRing: ConvertibleRing, const P: usize, const N: usize> From<[BaseRing; N]>
    for PowPCyclotomicPolyRing<BaseRing, P, N>
{
    fn from(value: [BaseRing; N]) -> Self {
        Self(Self::Inner::const_from_array(value))
    }
}

impl<BaseRing: ConvertibleRing, const P: usize, const N: usize> Modulus
    for PowPCyclotomicPolyRing<BaseRing, P, N>
{
    fn modulus() -> BigUint {
        BaseRing::modulus()
    }
}

macro_rules! impl_from_primitive_type {
    ($primitive_type: ty) => {
        impl<BaseRing: ConvertibleRing, const P: usize, const N: usize> From<$primitive_type>
            for PowPCyclotomicPolyRing<BaseRing, P, N>
        {
            fn from(value: $primitive_type) -> Self {
                Self::from_scalar(BaseRing::from(value))
            }
        }
    };
}

impl_from_primitive_type!(BaseRing);
impl_from_primitive_type!(u128);
impl_from_primitive_type!(u64);
impl_from_primitive_type!(u32);
impl_from_primitive_type!(u16);
impl_from_primitive_type!(u8);
impl_from_primitive_type!(bool);

// impl<BaseRing: ConvertibleRing, const N: usize> Default for PowPCyclotomicPolyRing<BaseRing, N> {
//     #[inline(always)]
//     fn default() -> Self {
//         Self::zero()
//     }
// }

// impl<BaseRing: ConvertibleRing, const P: usize, const N: usize> Zero
//     for PowPCyclotomicPolyRing<BaseRing, P, N>
// {
//     #[inline(always)]
//     fn zero() -> Self {
//         Self::ZERO
//     }
//
//     #[inline(always)]
//     fn is_zero(&self) -> bool {
//         self.eq(&Self::ZERO)
//     }
// }

// impl<BaseRing: ConvertibleRing, const N: usize> One for PowPCyclotomicPolyRing<BaseRing, N> {
//     #[inline(always)]
//     fn one() -> Self {
//         Self::ONE
//     }
// }

impl<BaseRing: ConvertibleRing, const P: usize, const N: usize> Mul<Self>
    for PowPCyclotomicPolyRing<BaseRing, P, N>
{
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        // This is for a power of two cyclotomic ring
        // let mut out = vec![BaseRing::zero(); N];
        // for i in 0..N {
        //     for j in 0..N - i {
        //         out[i + j] += self.0[i] * rhs.0[j];
        //     }
        //     for j in N - i..N {
        //         out[i + j - N] -= self.0[i] * rhs.0[j];
        //     }
        // }
        // Self::from(out)
        todo!()
    }
}

impl<BaseRing: ConvertibleRing, const P: usize, const N: usize> Ring
    for PowPCyclotomicPolyRing<BaseRing, P, N>
{
    const ZERO: Self = Self::const_from_element(BaseRing::ZERO);
    const ONE: Self = Self::const_from_element(BaseRing::ONE);
}

impl<BaseRing: ConvertibleRing, const P: usize, const N: usize> FromRandomBytes<Self>
    for PowPCyclotomicPolyRing<BaseRing, P, N>
{
    fn byte_size() -> usize {
        N * BaseRing::byte_size()
    }

    fn try_from_random_bytes(bytes: &[u8]) -> Option<Self> {
        assert_eq!(bytes.len(), Self::byte_size());
        let coeffs = core::array::from_fn(|i| {
            BaseRing::try_from_random_bytes(
                &bytes[i * BaseRing::byte_size()..(i + 1) * BaseRing::byte_size()],
            )
            .unwrap()
        });
        Some(Self::from(coeffs))
    }
}

// ============================================================================
// ============================================================================
// ============================================================================
// ============================================================================
// ============================================================================

impl<BaseRing: ConvertibleRing, const P: usize, const N: usize> Display
    for PowPCyclotomicPolyRing<BaseRing, P, N>
{
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        std::fmt::Display::fmt(&self.0, f)
    }
}

impl<BaseRing: ConvertibleRing, const P: usize, const N: usize> UniformRand
    for PowPCyclotomicPolyRing<BaseRing, P, N>
{
    fn rand<R: Rng + ?Sized>(rng: &mut R) -> Self {
        Self::from_fn(|_| BaseRing::rand(rng))
    }
}
// impl<'a, BaseRing: ConvertibleRing, const N: usize> Add<&'a Self>
//     for PowPCyclotomicPolyRing<BaseRing, N>
// {
//     type Output = Self;
//
//     fn add(self, rhs: &'a Self) -> Self::Output {
//         self.0.add(rhs.0).into()
//     }
// }

impl<BaseRing: ConvertibleRing, const P: usize, const N: usize> Mul<BaseRing>
    for PowPCyclotomicPolyRing<BaseRing, P, N>
{
    type Output = Self;

    fn mul(self, rhs: BaseRing) -> Self::Output {
        self.mul(Self::from_scalar(rhs))
    }
}
impl<BaseRing: ConvertibleRing, const P: usize, const N: usize> PolyRing
    for PowPCyclotomicPolyRing<BaseRing, P, N>
{
    type BaseRing = BaseRing;
    fn coeffs(&self) -> Vec<Self::BaseRing> {
        self.0.into_iter().map(|v_i| *v_i).collect()
    }
    fn dimension() -> usize {
        N
    }

    fn from_scalar(v: Self::BaseRing) -> Self {
        Self::from_fn(|i| if i == 0 { v } else { BaseRing::zero() })
    }
}

impl<BaseRing: ConvertibleRing, const P: usize, const N: usize> From<Vec<BaseRing>>
    for PowPCyclotomicPolyRing<BaseRing, P, N>
{
    fn from(value: Vec<BaseRing>) -> Self {
        Self::try_from(value).unwrap()
    }
}

impl<BaseRing: ConvertibleRing, const P: usize, const N: usize> WithConjugationAutomorphism
    for PowPCyclotomicPolyRing<BaseRing, P, N>
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

impl<BaseRing: ConvertibleRing, const P: usize, const N: usize> WithL2Norm
    for PowPCyclotomicPolyRing<BaseRing, P, N>
{
    fn l2_norm_squared(&self) -> u128 {
        self.coeffs().l2_norm_squared()
    }
}

impl<BaseRing: ConvertibleRing, const P: usize, const N: usize> WithLinfNorm
    for PowPCyclotomicPolyRing<BaseRing, P, N>
{
    fn linf_norm(&self) -> u128 {
        self.coeffs().linf_norm()
    }
}
