use std::fmt::{Debug, Display, Formatter};
use std::hash::Hash;
use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

use ark_ff::{Field, PrimeField};
use ark_std::rand::Rng;
use ark_std::UniformRand;
use derive_more::{Add, AddAssign, From, Into, Sub, SubAssign, Sum};
use nalgebra::{ArrayStorage, SVector};
use num_traits::{One, Zero};
use serde::{Deserialize, Deserializer, Serialize, Serializer};

use crate::ntt::NTT;
use crate::poly_ring::PolyRing;
use crate::pow2_cyclotomic_poly_ring::Pow2CyclotomicPolyRing;
use crate::ring::{Fq, Ring};
use crate::serde::{ark_de, ark_se};
use crate::traits::{
    FromRandomBytes, Modulus, WithConjugationAutomorphism, WithL2Norm, WithLinfNorm,
};

#[derive(
    Clone, Copy, Debug, Eq, PartialEq, Hash, Add, AddAssign, Sum, Sub, SubAssign, From, Into,
)]
pub struct Pow2CyclotomicPolyRingNTT<const Q: u64, const N: usize>(SVector<Fq<Q>, N>);

impl<const Q: u64, const N: usize> Pow2CyclotomicPolyRingNTT<Q, N> {
    /// Constructs a polynomial from a slice of coefficients in NTT form.
    fn from_slice(coeffs_ntt: &[Fq<Q>]) -> Self {
        // private since we can't enforce coeffs being in NTT form.
        assert_eq!(coeffs_ntt.len(), N);
        Self {
            0: SVector::<Fq<Q>, N>::from_row_slice(coeffs_ntt),
        }
    }

    /// Constructs a polynomial from a function specifying coefficients in non-NTT form.
    pub fn from_fn<F>(mut f: F) -> Self
    where
        F: FnMut(usize) -> Fq<Q>,
    {
        let mut coeffs = (0..N).map(|i| f(i)).collect::<Vec<Fq<Q>>>();
        Self::ntt(coeffs.as_mut_slice());
        Self::from_slice(&coeffs)
    }
}

impl<const Q: u64, const N: usize> const NTT<Q, N> for Pow2CyclotomicPolyRingNTT<Q, N> {
    fn ntt_coeffs(&self) -> Vec<Fq<Q>> {
        self.0.iter().cloned().collect()
    }
}

impl<const Q: u64, const N: usize> Modulus for Pow2CyclotomicPolyRingNTT<Q, N> {
    fn modulus() -> u64 {
        Fq::<Q>::MODULUS.0[0]
    }
}

const fn vec_from_element<const Q: u64, const N: usize>(elem: Fq<Q>) -> SVector<Fq<Q>, N> {
    let coeffs = [elem; N];
    SVector::<Fq<Q>, N>::from_array_storage(ArrayStorage::<Fq<Q>, { N }, 1> { 0: [coeffs; 1] })
}

impl<const Q: u64, const N: usize> Ring for Pow2CyclotomicPolyRingNTT<Q, N> {
    const ZERO: Self = Self {
        0: vec_from_element(Fq::<Q>::ZERO),
    };
    const ONE: Self = Self {
        0: vec_from_element(Fq::<Q>::ONE),
    };
}

impl<const Q: u64, const N: usize> FromRandomBytes<Self> for Pow2CyclotomicPolyRingNTT<Q, N> {
    fn byte_size() -> usize {
        N * Fq::<Q>::byte_size()
    }

    fn try_from_random_bytes(bytes: &[u8]) -> Option<Self> {
        assert_eq!(bytes.len(), Self::byte_size());
        let coeffs = SVector::<Fq<Q>, N>::from_fn(|i, _| {
            Fq::<Q>::try_from_random_bytes(
                &bytes[i * Fq::<Q>::byte_size()..(i + 1) * Fq::<Q>::byte_size()],
            )
            .unwrap()
        });
        Some(Self::from(coeffs))
    }
}

impl<const Q: u64, const N: usize> Serialize for Pow2CyclotomicPolyRingNTT<Q, N> {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        ark_se(&self.coeffs(), serializer)
    }
}

impl<'a, const Q: u64, const N: usize> Deserialize<'a> for Pow2CyclotomicPolyRingNTT<Q, N> {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'a>,
    {
        ark_de(deserializer).map(|v: Vec<Fq<Q>>| Self::from(v))
    }
}

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
        Self::from_fn(|_| Fq::<Q>::rand(rng))
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
        self.0.add(&rhs.0).into()
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
        self.0.add_assign(&rhs.0)
    }
}

impl<'a, const Q: u64, const N: usize> SubAssign<&'a Self> for Pow2CyclotomicPolyRingNTT<Q, N> {
    fn sub_assign(&mut self, rhs: &'a Self) {
        self.0.sub_assign(&rhs.0)
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
        self.0.add(&rhs.0).into()
    }
}

impl<'a, const Q: u64, const N: usize> Sub<&'a mut Self> for Pow2CyclotomicPolyRingNTT<Q, N> {
    type Output = Self;

    fn sub(self, rhs: &'a mut Self) -> Self::Output {
        self.0.sub(&rhs.0).into()
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
        self.0.add_assign(&rhs.0).into()
    }
}

impl<'a, const Q: u64, const N: usize> SubAssign<&'a mut Self> for Pow2CyclotomicPolyRingNTT<Q, N> {
    fn sub_assign(&mut self, rhs: &'a mut Self) {
        self.0.sub_assign(&rhs.0).into()
    }
}

impl<'a, const Q: u64, const N: usize> MulAssign<&'a mut Self> for Pow2CyclotomicPolyRingNTT<Q, N> {
    fn mul_assign(&mut self, rhs: &'a mut Self) {
        self.0.component_mul_assign(&rhs.0)
    }
}

impl<const Q: u64, const N: usize> From<Pow2CyclotomicPolyRing<Fq<Q>, N>>
    for Pow2CyclotomicPolyRingNTT<Q, N>
{
    fn from(value: Pow2CyclotomicPolyRing<Fq<Q>, N>) -> Self {
        let mut coeffs = value.coeffs();
        Self::ntt(coeffs.as_mut_slice());
        Self::from_slice(coeffs.as_slice())
    }
}

impl<const Q: u64, const N: usize> Into<Pow2CyclotomicPolyRing<Fq<Q>, N>>
    for Pow2CyclotomicPolyRingNTT<Q, N>
{
    fn into(self) -> Pow2CyclotomicPolyRing<Fq<Q>, N> {
        Pow2CyclotomicPolyRing::<Fq<Q>, N>::from(self.coeffs())
    }
}

impl<const Q: u64, const N: usize> Mul<Fq<Q>> for Pow2CyclotomicPolyRingNTT<Q, N> {
    type Output = Self;

    fn mul(self, rhs: Fq<Q>) -> Self::Output {
        self.mul(Self::from_scalar(rhs))
    }
}

impl<const Q: u64, const N: usize> PolyRing for Pow2CyclotomicPolyRingNTT<Q, N> {
    type BaseRing = Fq<Q>;
    fn coeffs(&self) -> Vec<Fq<Q>> {
        let mut coeffs = self.ntt_coeffs();
        Self::intt(&mut coeffs);
        coeffs
    }
    fn dimension() -> usize {
        N
    }

    fn from_scalar(v: Self::BaseRing) -> Self {
        // NTT([v, 0, ..., 0]) = ([v, ..., v])
        Self {
            0: SVector::<Fq<Q>, N>::from_fn(|_i, _| v),
        }
    }
}

impl<const Q: u64, const N: usize> From<Vec<Fq<Q>>> for Pow2CyclotomicPolyRingNTT<Q, N> {
    fn from(mut value: Vec<Fq<Q>>) -> Self {
        Self::ntt(value.as_mut_slice());
        Self {
            0: SVector::<Fq<Q>, N>::from_vec(value),
        }
    }
}

impl<const Q: u64, const N: usize> From<Fq<Q>> for Pow2CyclotomicPolyRingNTT<Q, N> {
    fn from(value: Fq<Q>) -> Self {
        Self::from_scalar(value)
    }
}

impl<const Q: u64, const N: usize> From<u128> for Pow2CyclotomicPolyRingNTT<Q, N> {
    fn from(value: u128) -> Self {
        Self::from_scalar(value.into())
    }
}

impl<const Q: u64, const N: usize> WithConjugationAutomorphism for Pow2CyclotomicPolyRingNTT<Q, N> {
    fn sigma(&self) -> Self {
        // TODO: can we implement the automorphism directly in NTT form?
        Into::<Pow2CyclotomicPolyRing<Fq<Q>, N>>::into(*self)
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
