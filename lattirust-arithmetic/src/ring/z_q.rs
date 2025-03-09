#![allow(warnings)] // TODO: remove
#![allow(long_running_const_eval)]

use std::array;
use std::convert::Into;
use std::fmt::{Debug};
use std::hash::Hash;
use std::io::{Read, Write};
use std::iter::Iterator;
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

use ark_ff::{
    AdditiveGroup, BigInt, BigInteger, FftField, Field, Fp, Fp64, FpConfig, MontBackend,
    MontConfig, PrimeField,
};
use ark_serialize::{
    CanonicalDeserialize, CanonicalDeserializeWithFlags, CanonicalSerialize,
    CanonicalSerializeWithFlags, Compress, EmptyFlags, Flags, SerializationError, Valid, Validate,
};
use ark_std::rand::distributions::Distribution;
use ark_std::rand::Rng;
use ark_std::UniformRand;
use derivative::Derivative;
use derive_more::Display;
use num_bigint::BigUint;
use num_traits::{One, Signed, Zero};
use rounded_div::RoundedDiv;
use zeroize::Zeroize;

use crate::balanced_decomposition::DecompositionFriendlySignedRepresentative;
use crate::impl_try_from_primitive_type;
use crate::ring::{NttRing, Ring};
use crate::ring::ntt::Ntt;
use crate::ring::representatives::WithSignedRepresentative;
use crate::ring::z_q_signed_representative::ZqSignedRepresentative;
use crate::traits::{FromRandomBytes, Modulus, WithLinfNorm};

/// Returns an array containing the prime factors of `n`.
/// The length of the array is fixed to 64, with remaining slots filled with 0s.
/// This is necessary since we cannot return a dynamically-sized array in a `const fn`.
pub const fn prime_factors(mut n: u64) -> [u64; 64] {
    let mut factors = [0u64; 64];
    let mut index = 0;
    let mut divisor = 2;

    while n > 1 {
        if n % divisor == 0 {
            factors[index] = divisor;
            index += 1;
            n /= divisor;
        } else {
            divisor = const_primes::next_prime(divisor + 1).unwrap();
        }
    }

    factors
}

/*
const fn generator<const Q: u64>() -> u64 {
    assert!(const_primes::is_prime(Q));

    let mut quotients = prime_factors(Q - 1);
    let mut j: usize = 1;
    while j < 64 && quotients[0] != 0 {
        quotients[j] = (Q - 1) / quotients[j];
    }

    let mut next_pow_2 = 1;
    // while next_pow_2 * 2 > Q {
    //     next_pow_2 *= 2;
    // }

    let mut i = 1;
    while i < Q {
        let mut idx = 0;
        let mut is_generator = true;
        while quotients[idx] != 0 {
            if const_pow_mod::<Q>(i, quotients[idx]) == 1 {
                is_generator = false;
                break;
            }
            idx += 1;
        }
        if is_generator {
            // Check that generator generates primitive roots of unity for 2^k, 0 < k < ceil(log2(Q))
            for k in 0..next_pow_2 {
                if !is_primitive_root_of_unity::<Q>(i, k) {
                    break;
                }
            }
            return i;
        }
        i += 1;
    }
    panic!("No generator found");
}
*/

pub struct FqConfig<const Q: u64> {}

impl<const Q: u64> FqConfig<Q> {
    /// If `_IS_ODD_PRIME` is evaluated somewhere, this will cause a compilation error if FqConfig is instantiated with `Q` not an odd prime.
    const _IS_ODD_PRIME: () = assert!(
        Q > 2 && const_primes::is_prime(Q),
        "You tried to instantiate an FqConfig<Q> with either Q = 2 or Q not a prime"
    );
}

/// Hacky way to ensure that `_IS_ODD_PRIME` is evaluated at compile time, and that the assertion is actually called.
const fn to_bigint_assert_odd_prime<const Q: u64>() -> BigInt<1> {
    let _ = FqConfig::<Q>::_IS_ODD_PRIME;
    BigInt::<1>([Q])
}

impl<const Q: u64> MontConfig<1> for FqConfig<Q> {
    const MODULUS: BigInt<1> = to_bigint_assert_odd_prime::<Q>();

    // TODO: As far as I can tell, `GENERATOR`and `TWO_ADIC_ROOT_OF_UNITY` are only used for FftField. For our small Fft sizes, we might be able to significantly speed up the finding of a generator.
    const GENERATOR: Fp<MontBackend<Self, 1>, 1> = Fp::new(BigInt::<1>([0u64])); //generator::<Q>()]));
    const TWO_ADIC_ROOT_OF_UNITY: Fp<MontBackend<Self, 1>, 1> = Fp::new(BigInt::<1>([0u64]));
}

/// `Fq<Q>` is a prime field with modulus `Q`, where `Q` is less than 64 bits.
pub type Fq<const Q: u64> = Fp64<MontBackend<FqConfig<Q>, 1>>;

pub const fn fq_zero<const Q: u64>() -> Fq<Q> {
    MontBackend::<FqConfig<Q>, 1>::ZERO
}

pub trait ZqConfig<const L: usize>: Send + Sync + 'static + Sized {
    const MODULI: [u64; L];
    type Limbs: Sync
        + Send
        + Sized
        + Clone
        + Copy
        + Debug
        + Hash
        + PartialEq
        + Eq
        + Zeroize
        + CanonicalSerialize
        + CanonicalDeserialize
        + FromRandomBytes<Self::Limbs>;

    const ZERO: Zq<Self, L>;
    const ONE: Zq<Self, L>;

    /// Set a += b.
    fn add_assign(a: &mut Zq<Self, L>, b: &Zq<Self, L>);

    /// Set a -= b.
    fn sub_assign(a: &mut Zq<Self, L>, b: &Zq<Self, L>);

    /// Set a = a + a.
    fn double_in_place(a: &mut Zq<Self, L>);

    /// Set a = -a;
    fn neg_in_place(a: &mut Zq<Self, L>);

    /// Set a *= b.
    fn mul_assign(a: &mut Zq<Self, L>, b: &Zq<Self, L>);

    /// Compute the inner product `<a, b>`.
    fn sum_of_products<const T: usize>(a: &[Zq<Self, L>; T], b: &[Zq<Self, L>; T]) -> Zq<Self, L>;

    /// Set a *= b.
    fn square_in_place(a: &mut Zq<Self, L>);

    /// Compute a^{-1} if `a` is not zero.
    fn inverse(a: &Zq<Self, L>) -> Option<Zq<Self, L>>;

    /// Construct a field element from an integer in the range `0..(Self::MODULUS - 1)`. Returns `None` if the integer is outside this range.
    fn from_bigint(other: BigInt<L>) -> Option<Zq<Self, L>>;

    /// Convert a field element to an integer in the range `0..(Self::MODULUS - 1)`.
    fn into_bigint(other: Zq<Self, L>) -> BigInt<L>;

    fn rand<R: Rng + ?Sized>(rng: &mut R) -> Zq<Self, L>;
}

#[derive(Derivative, PartialOrd, Ord, Zeroize, Display)]
#[derivative(
    Debug(bound = ""),
    Hash(bound = ""),
    Clone(bound = ""),
    Copy(bound = ""),
    PartialEq(bound = ""),
    Eq(bound = "")
)]
#[display("{}", C::into_bigint(*self))]
pub struct Zq<C: ZqConfig<L>, const L: usize>(C::Limbs);

impl<C: ZqConfig<L>, const L: usize> Zq<C, L> {}

impl<C: ZqConfig<L>, const L: usize> Default for Zq<C, L> {
    fn default() -> Self {
        C::ZERO
    }
}

impl<C: ZqConfig<L>, const L: usize> Zero for Zq<C, L> {
    fn zero() -> Self {
        C::ZERO
    }

    fn is_zero(&self) -> bool {
        *self == C::ZERO
    }
}

impl<C: ZqConfig<L>, const L: usize> One for Zq<C, L> {
    fn one() -> Self {
        C::ONE
    }
}

impl<C: ZqConfig<L>, const L: usize> Ring for Zq<C, L> {
    const ZERO: Self = C::ZERO;
    const ONE: Self = C::ONE;

    fn inverse(&self) -> Option<Self> {
        C::inverse(self)
    }
}

impl<C: ZqConfig<L>, const L: usize> From<bool> for Zq<C, L> {
    fn from(b: bool) -> Self {
        if b {
            Self::one()
        } else {
            Self::zero()
        }
    }
}

impl<C: ZqConfig<L>, const L: usize> TryFrom<BigUint> for Zq<C, L> {
    type Error = ();

    fn try_from(value: BigUint) -> Result<Self, Self::Error> {
        match C::from_bigint(BigInt::<L>::try_from(value)?) {
            None => Err(()),
            Some(elem) => Ok(elem),
        }
    }
}

#[macro_export]
macro_rules! impl_try_from_primitive_type {
    ($primitive_type: ty) => {
        impl<C: ZqConfig<L>, const L: usize> TryFrom<$primitive_type> for Zq<C, L> {
            type Error = <Zq<C, L> as TryFrom<BigUint>>::Error;

            fn try_from(value: $primitive_type) -> Result<Self, Self::Error> {
                Self::try_from(BigUint::from(value))
            }
        }
    };
}

impl_try_from_primitive_type!(u8);
impl_try_from_primitive_type!(u16);
impl_try_from_primitive_type!(u32);
impl_try_from_primitive_type!(u64);
impl_try_from_primitive_type!(u128);

impl<C: ZqConfig<L>, const L: usize> Neg for Zq<C, L> {
    type Output = Self;
    #[inline]
    #[must_use]
    fn neg(mut self) -> Self {
        C::neg_in_place(&mut self);
        self
    }
}

impl<'a, C: ZqConfig<L>, const L: usize> Add<&'a Zq<C, L>> for Zq<C, L> {
    type Output = Self;

    #[inline]
    fn add(mut self, other: &Self) -> Self {
        self.add_assign(other);
        self
    }
}

impl<'a, C: ZqConfig<L>, const L: usize> Sub<&'a Zq<C, L>> for Zq<C, L> {
    type Output = Self;

    #[inline]
    fn sub(mut self, other: &Self) -> Self {
        self.sub_assign(other);
        self
    }
}

impl<'a, C: ZqConfig<L>, const L: usize> Mul<&'a Zq<C, L>> for Zq<C, L> {
    type Output = Self;

    #[inline]
    fn mul(mut self, other: &Self) -> Self {
        self.mul_assign(other);
        self
    }
}

impl<'a, C: ZqConfig<L>, const L: usize> Div<&'a Zq<C, L>> for Zq<C, L> {
    type Output = Self;

    /// Returns `self * other.inverse()` if `other.inverse()` is `Some`, and
    /// panics otherwise.
    #[inline]
    fn div(mut self, other: &Self) -> Self {
        self.mul_assign(&other.inverse().unwrap());
        self
    }
}

impl<'a, 'b, C: ZqConfig<L>, const L: usize> Add<&'b Zq<C, L>> for &'a Zq<C, L> {
    type Output = Zq<C, L>;

    #[inline]
    fn add(self, other: &'b Zq<C, L>) -> Zq<C, L> {
        let mut result = *self;
        result.add_assign(other);
        result
    }
}

impl<'a, 'b, C: ZqConfig<L>, const L: usize> Sub<&'b Zq<C, L>> for &'a Zq<C, L> {
    type Output = Zq<C, L>;

    #[inline]
    fn sub(self, other: &Zq<C, L>) -> Zq<C, L> {
        let mut result = *self;
        result.sub_assign(other);
        result
    }
}

impl<'a, 'b, C: ZqConfig<L>, const L: usize> Mul<&'b Zq<C, L>> for &'a Zq<C, L> {
    type Output = Zq<C, L>;

    #[inline]
    fn mul(self, other: &Zq<C, L>) -> Zq<C, L> {
        let mut result = *self;
        result.mul_assign(other);
        result
    }
}

impl<'a, 'b, C: ZqConfig<L>, const L: usize> Div<&'b Zq<C, L>> for &'a Zq<C, L> {
    type Output = Zq<C, L>;

    #[inline]
    fn div(self, other: &Zq<C, L>) -> Zq<C, L> {
        let mut result = *self;
        result.div_assign(other);
        result
    }
}

impl<'a, C: ZqConfig<L>, const L: usize> AddAssign<&'a Self> for Zq<C, L> {
    #[inline]
    fn add_assign(&mut self, other: &Self) {
        C::add_assign(self, other)
    }
}

impl<'a, C: ZqConfig<L>, const L: usize> SubAssign<&'a Self> for Zq<C, L> {
    #[inline]
    fn sub_assign(&mut self, other: &Self) {
        C::sub_assign(self, other);
    }
}

#[allow(unused_qualifications)]
impl<C: ZqConfig<L>, const L: usize> core::ops::Add<Self> for Zq<C, L> {
    type Output = Self;

    #[inline]
    fn add(mut self, other: Self) -> Self {
        self.add_assign(&other);
        self
    }
}

#[allow(unused_qualifications)]
impl<'a, C: ZqConfig<L>, const L: usize> core::ops::Add<&'a mut Self> for Zq<C, L> {
    type Output = Self;

    #[inline]
    fn add(mut self, other: &'a mut Self) -> Self {
        self.add_assign(&*other);
        self
    }
}

#[allow(unused_qualifications)]
impl<C: ZqConfig<L>, const L: usize> core::ops::Sub<Self> for Zq<C, L> {
    type Output = Self;

    #[inline]
    fn sub(mut self, other: Self) -> Self {
        self.sub_assign(&other);
        self
    }
}

#[allow(unused_qualifications)]
impl<'a, C: ZqConfig<L>, const L: usize> core::ops::Sub<&'a mut Self> for Zq<C, L> {
    type Output = Self;

    #[inline]
    fn sub(mut self, other: &'a mut Self) -> Self {
        self.sub_assign(&*other);
        self
    }
}

#[allow(unused_qualifications)]
impl<C: ZqConfig<L>, const L: usize> core::iter::Sum<Self> for Zq<C, L> {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::zero(), core::ops::Add::add)
    }
}

#[allow(unused_qualifications)]
impl<'a, C: ZqConfig<L>, const L: usize> core::iter::Sum<&'a Self> for Zq<C, L> {
    fn sum<I: Iterator<Item = &'a Self>>(iter: I) -> Self {
        iter.fold(Self::zero(), core::ops::Add::add)
    }
}

#[allow(unused_qualifications)]
impl<C: ZqConfig<L>, const L: usize> core::ops::AddAssign<Self> for Zq<C, L> {
    #[inline(always)]
    fn add_assign(&mut self, other: Self) {
        self.add_assign(&other)
    }
}

#[allow(unused_qualifications)]
impl<C: ZqConfig<L>, const L: usize> core::ops::SubAssign<Self> for Zq<C, L> {
    #[inline(always)]
    fn sub_assign(&mut self, other: Self) {
        self.sub_assign(&other)
    }
}

#[allow(unused_qualifications)]
impl<'a, C: ZqConfig<L>, const L: usize> core::ops::AddAssign<&'a mut Self> for Zq<C, L> {
    #[inline(always)]
    fn add_assign(&mut self, other: &'a mut Self) {
        self.add_assign(&*other)
    }
}

#[allow(unused_qualifications)]
impl<'a, C: ZqConfig<L>, const L: usize> core::ops::SubAssign<&'a mut Self> for Zq<C, L> {
    #[inline(always)]
    fn sub_assign(&mut self, other: &'a mut Self) {
        self.sub_assign(&*other)
    }
}

impl<'a, C: ZqConfig<L>, const L: usize> MulAssign<&'a Self> for Zq<C, L> {
    fn mul_assign(&mut self, other: &Self) {
        C::mul_assign(self, other)
    }
}

/// Computes `self *= other.inverse()` if `other.inverse()` is `Some`, and
/// panics otherwise.
impl<'a, C: ZqConfig<L>, const L: usize> DivAssign<&'a Self> for Zq<C, L> {
    #[inline(always)]
    fn div_assign(&mut self, other: &Self) {
        self.mul_assign(&other.inverse().unwrap());
    }
}

#[allow(unused_qualifications)]
impl<C: ZqConfig<L>, const L: usize> core::ops::Mul<Self> for Zq<C, L> {
    type Output = Self;

    #[inline(always)]
    fn mul(mut self, other: Self) -> Self {
        self.mul_assign(&other);
        self
    }
}

#[allow(unused_qualifications)]
impl<C: ZqConfig<L>, const L: usize> core::ops::Div<Self> for Zq<C, L> {
    type Output = Self;

    #[inline(always)]
    fn div(mut self, other: Self) -> Self {
        self.div_assign(&other);
        self
    }
}

#[allow(unused_qualifications)]
impl<'a, C: ZqConfig<L>, const L: usize> core::ops::Mul<&'a mut Self> for Zq<C, L> {
    type Output = Self;

    #[inline(always)]
    fn mul(mut self, other: &'a mut Self) -> Self {
        self.mul_assign(&*other);
        self
    }
}

#[allow(unused_qualifications)]
impl<'a, C: ZqConfig<L>, const L: usize> core::ops::Div<&'a mut Self> for Zq<C, L> {
    type Output = Self;

    #[inline(always)]
    fn div(mut self, other: &'a mut Self) -> Self {
        self.div_assign(&*other);
        self
    }
}

#[allow(unused_qualifications)]
impl<C: ZqConfig<L>, const L: usize> core::iter::Product<Self> for Zq<C, L> {
    fn product<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::one(), core::ops::Mul::mul)
    }
}

#[allow(unused_qualifications)]
impl<'a, C: ZqConfig<L>, const L: usize> core::iter::Product<&'a Self> for Zq<C, L> {
    fn product<I: Iterator<Item = &'a Self>>(iter: I) -> Self {
        iter.fold(Self::one(), Mul::mul)
    }
}

#[allow(unused_qualifications)]
impl<C: ZqConfig<L>, const L: usize> core::ops::MulAssign<Self> for Zq<C, L> {
    #[inline(always)]
    fn mul_assign(&mut self, other: Self) {
        self.mul_assign(&other)
    }
}

#[allow(unused_qualifications)]
impl<'a, C: ZqConfig<L>, const L: usize> core::ops::DivAssign<&'a mut Self> for Zq<C, L> {
    #[inline(always)]
    fn div_assign(&mut self, other: &'a mut Self) {
        self.div_assign(&*other)
    }
}

#[allow(unused_qualifications)]
impl<'a, C: ZqConfig<L>, const L: usize> core::ops::MulAssign<&'a mut Self> for Zq<C, L> {
    #[inline(always)]
    fn mul_assign(&mut self, other: &'a mut Self) {
        self.mul_assign(&*other)
    }
}

#[allow(unused_qualifications)]
impl<C: ZqConfig<L>, const L: usize> core::ops::DivAssign<Self> for Zq<C, L> {
    #[inline(always)]
    fn div_assign(&mut self, other: Self) {
        self.div_assign(&other)
    }
}

impl<C: ZqConfig<L>, const L: usize> Modulus for Zq<C, L> {
    fn modulus() -> BigUint {
        C::MODULI.iter().map(|q| BigUint::from(*q)).product()
    }
}

impl<C: ZqConfig<L>, const L: usize> FromRandomBytes<Zq<C, L>> for Zq<C, L> 
{
    fn needs_bytes() -> usize {
        C::Limbs::needs_bytes()
    }

    fn try_from_random_bytes_inner(bytes: &[u8]) -> Option<Zq<C, L>> {
        C::Limbs::try_from_random_bytes_inner(bytes).map(Self)
    }
}

impl<C: ZqConfig<L>, const L: usize> CanonicalSerialize for Zq<C, L> {
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

impl<C: ZqConfig<L>, const L: usize> Valid for Zq<C, L> {
    fn check(&self) -> Result<(), SerializationError> {
        self.0.check()
    }
}

impl<C: ZqConfig<L>, const L: usize> CanonicalDeserialize for Zq<C, L> {
    fn deserialize_with_mode<R: Read>(
        reader: R,
        compress: Compress,
        validate: Validate,
    ) -> Result<Self, SerializationError> {
       C::Limbs::deserialize_with_mode(reader, compress, validate).map(Self)
    }
}

impl<C: ZqConfig<L>, const L: usize> Distribution<Zq<C, L>> for Zq<C, L> {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> Zq<C, L> {
        todo!()
    }
}

impl<C: ZqConfig<L>, const L: usize> UniformRand for Zq<C, L> {
    fn rand<R: Rng + ?Sized>(rng: &mut R) -> Self {
        C::rand(rng)
    }
}

/// Given a set of co-prime RNS moduli `qs = [q1, ..., qL]`, compute the RNS coefficients `mod_inv(C, qi) * C` for each modulus, where `C = (q1 * ... * qL) / qi`.
///
/// This is using BigUint and does not implement any fancy tricks.  
fn rns_coeffs<const L: usize>(qs: [u64; L]) -> [BigUint; L] {
    let mut res: [BigUint; L] = array::from_fn(|i| BigUint::one());
    let qs_biguint: [BigUint; L] = array::from_fn(|i| BigUint::from(qs[i]));
    for i in 0..L {
        for j in 0..L {
            if j != i {
                res[i] *= qs_biguint[j].clone();
            }
        }
    }
    for i in 0..L {
        res[i] *= res[i].modinv(&qs_biguint[i]).unwrap();
    }
    res
}

/// Map (-h, h] to [0, q) using [0, h] -> [0, h] and (-h, 0) -> (h, q), where h = floor((q-1)/2)
impl<C: ZqConfig<L>, const L: usize> From<ZqSignedRepresentative> for Zq<C, L> {
    fn from(value: ZqSignedRepresentative) -> Self {
        if value.0.is_negative() {
            C::from_bigint(
                BigInt::<L>::try_from(Zq::<C, L>::modulus() - BigUint::try_from(-value.0).unwrap())
                    .unwrap(),
            )
            .unwrap()
        } else {
            C::from_bigint(BigInt::<L>::try_from(BigUint::try_from(value.0).unwrap()).unwrap())
                .unwrap()
        }
    }
}

/// Map [0, q) to (-h, h] using [0, h] -> [0, h] and (h, q) -> (-h, 0), where h = floor((q-1)/2)
impl<C: ZqConfig<L>, const L: usize> From<Zq<C, L>> for ZqSignedRepresentative {
    fn from(value: Zq<C, L>) -> Self {
        let q_half = (Zq::<C, L>::modulus() - BigUint::one()) / BigUint::from(2u8);
        let v = num_bigint::BigInt::from(BigUint::from(C::into_bigint(value)));
        if v > num_bigint::BigInt::from(q_half) {
            ZqSignedRepresentative(v - num_bigint::BigInt::from(Zq::<C, L>::modulus()))
        } else {
            ZqSignedRepresentative(v)
        }
    }
}

impl<C: ZqConfig<L>, const L: usize> WithSignedRepresentative for Zq<C, L> {
    type SignedRepresentative = ZqSignedRepresentative;
}

/// To be called as `zq_config_impl!(L, 0, 1, ..., L-1);` for some non-negative number of limbs `L`
macro_rules! zq_config_impl {
    ($L:literal $(,$i:literal)+) => {
        paste::expr! {
            pub struct [< Zq $L ConfigImpl >]<$(const [< Q $i >]: u64,)*>;

            impl<$(const [< Q $i >]: u64,)*> [< Zq $L ConfigImpl >]<$([< Q $i >],)*> {
                $(
                type [< F $i >] = Fq<[< Q $i >]>;
                type [< FConfig $i >] = FqConfig<[< Q $i >]>;
                )*
            }

            impl<$(const [< Q $i >]: u64,)*> ZqConfig<$L> for [< Zq $L ConfigImpl >]<$([< Q $i >],)*> {
                const MODULI: [u64; $L] = [$([< Q $i >],)*];
                type Limbs = ($(Fq<[< Q $i >]>,)*);
                const ZERO: Zq<Self, $L> = Zq(($(<Fq::<[< Q $i >]> as AdditiveGroup>::ZERO,)*));
                const ONE: Zq<Self, $L> = Zq(($(<Fq::<[< Q $i >]> as Field>::ONE,)*));

                fn add_assign(a: &mut Zq<Self, $L>, b: &Zq<Self, $L>) {
                    $(Self::[< FConfig $i >]::add_assign(&mut a.0.$i, &b.0.$i);)*
                }

                fn sub_assign(a: &mut Zq<Self, $L>, b: &Zq<Self, $L>) {
                    $(Self::[< FConfig $i >]::sub_assign(&mut a.0.$i, &b.0.$i);)*
                }

                fn double_in_place(a: &mut Zq<Self, $L>) {
                    $(Self::[< FConfig $i >]::double_in_place(&mut a.0.$i);)*
                }

                fn neg_in_place(a: &mut Zq<Self, $L>) {
                    $(Self::[< FConfig $i >]::neg_in_place(&mut a.0.$i);)*
                }

                fn mul_assign(a: &mut Zq<Self, $L>, b: &Zq<Self, $L>) {
                    $(Self::[< FConfig $i >]::mul_assign(&mut a.0.$i, &b.0.$i);)*
                }

                fn sum_of_products<const T: usize>(a: &[Zq<Self, $L>; T], b: &[Zq<Self, $L>; T]) -> Zq<Self, $L> {
                    todo!()
                }

                fn square_in_place(a: &mut Zq<Self, $L>) {
                    $(Self::[< FConfig $i >]::square_in_place(&mut a.0.$i);)*
                }

                fn inverse(a: &Zq<Self, $L>) -> Option<Zq<Self, $L>> {
                    $(
                    let [< inv_ $i >] = Self::[< FConfig $i >]::inverse(& a.0.$i)?;
                    )*
                    Some(Zq::<Self, $L>(($( [< inv_ $i >] ,)+)))
                }


                fn from_bigint(other: BigInt<$L>) -> Option<Zq<Self, $L>> {
                    let biguint = BigUint::from(other);
                    if biguint >= Zq::<Self, $L>::modulus() {
                        return None;
                    }
                    // Use BigUint to perform modular reductions for each modulus
                    Some(Zq::<Self, $L>((
                        $( Self::[< F $i >]::new(
                            BigInt::try_from(&biguint % &BigUint::from([< Q $i >])).unwrap()
                        ) ,)+
                    )))
                }

                fn into_bigint(other: Zq<Self, $L>) -> BigInt<$L> {
                    // Recomputed at runtime because it's cumbersome to generate at compile time; may revisit this at a future point if this turns out to be a performance bottleneck.
                    let rns_coeffs = rns_coeffs::<$L>(Self::MODULI);

                    let prod = $(
                        &rns_coeffs[$i] * &BigUint::from(Self::[< FConfig $i >]::into_bigint(other.0.$i)) +
                    )+ &BigUint::zero();
                    BigInt::<$L>::try_from(
                        prod % Zq::<Self, $L>::modulus()
                    ).unwrap()
                }

                fn rand<R: Rng + ?Sized>(rng: &mut R) -> Zq<Self, $L> {
                    Zq::<Self, $L>(($(Self::[< F $i >]::rand(rng),)*))
                }
            }

            impl<$(const [< Q $i >]: u64,)* const N: usize> Ntt<N> for Zq<[< Zq $L ConfigImpl >]<$([< Q $i >],)*>, $L>
                where $(Fq<[< Q $i >]>: Ntt<N>,)*
            {
                fn ntt(coeffs: &mut [Self; N]) {
                    $(Fq::<[< Q $i >]>::ntt(&mut coeffs.map(|x| x.0 .$i));)*
                }

                fn intt(evals: &mut [Self; N]) {
                    $(Fq::<[< Q $i >]>::intt(&mut evals.map(|x| x.0 .$i));)*
                }
            }
        }
    }
}

zq_config_impl!(1, 0);
zq_config_impl!(2, 0, 1);
zq_config_impl!(3, 0, 1, 2);
zq_config_impl!(4, 0, 1, 2, 3);
zq_config_impl!(5, 0, 1, 2, 3, 4);
// zq_config_impl!(6, 0, 1, 2, 3, 4, 5);
// zq_config_impl!(7, 0, 1, 2, 3, 4, 5, 6);
// zq_config_impl!(8, 0, 1, 2, 3, 4, 5, 6, 7);
// zq_config_impl!(9, 0, 1, 2, 3, 4, 5, 6, 7, 8);
// zq_config_impl!(10, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9);

pub type Zq1<const Q: u64> = Zq<Zq1ConfigImpl<Q>, 1>;
pub type Zq2<const Q1: u64, const Q2: u64> = Zq<Zq2ConfigImpl<Q1, Q2>, 2>;
pub type Zq3<const Q1: u64, const Q2: u64, const Q3: u64> = Zq<Zq3ConfigImpl<Q1, Q2, Q3>, 3>;
pub type Zq4<const Q1: u64, const Q2: u64, const Q3: u64, const Q4: u64> =
    Zq<Zq4ConfigImpl<Q1, Q2, Q3, Q4>, 4>;
pub type Zq5<const Q1: u64, const Q2: u64, const Q3: u64, const Q4: u64, const Q5: u64> =
    Zq<Zq5ConfigImpl<Q1, Q2, Q3, Q4, Q5>, 5>;

macro_rules! test_zq_config_impl {
    ($L:literal $(,$Q:ident)+) => {
        paste::expr! {
            #[test]
            fn [< test_zq_config_impl_from_into_bigint $L >]() {
                let modulus: BigUint = Zq::<[< Zq $L ConfigImpl >]<$($Q,)*>, $L>::modulus();
                for x in [0, 1, $($Q/2, $Q-1, $Q, $Q+1,)*] {
                    let x_bigint = BigInt::<$L>::from(x);

                    let x_zq = [< Zq $L ConfigImpl >]::<$($Q,)*>::from_bigint(x_bigint);

                    if BigUint::from(x) >= modulus {
                        assert!(x_zq.is_none());
                    } else {
                        assert!(x_zq.is_some());
                        let x_out = [< Zq $L ConfigImpl >]::<$($Q,)*>::into_bigint(x_zq.unwrap());
                        if (x_bigint != x_out) {
                            dbg!(x);
                            dbg!(x_bigint.clone());
                            dbg!(x_zq.clone());
                            dbg!(x_out.clone());
                        }
                        assert_eq!(x_bigint, x_out);
                    }
                }
            }
            
            #[test]
            fn [< test_zq_config_impl_canonical_serialize_deserialize_compressed $L >]() {
                let modulus: BigUint = Zq::<[< Zq $L ConfigImpl >]<$($Q,)*>, $L>::modulus();
                for x in [0, 1, $($Q/2, $Q-1,)*] {
                    let x_in = Zq::<[< Zq $L ConfigImpl >]<$($Q,)*>, $L>::try_from(x).unwrap();
                    let mut bytes = Vec::new();
                    x_in.serialize_with_mode(&mut bytes, Compress::No).unwrap();
                    let mut x_out = Zq::<[< Zq $L ConfigImpl >]<$($Q,)*>, $L>::deserialize_with_mode(&*bytes, Compress::No, Validate::Yes).unwrap();
                    assert_eq!(x_in, x_out);
                }
            }

            #[test]
            fn [< test_zq_config_impl_canonical_serialize_deserialize_uncompressed $L >]() {
                let modulus: BigUint = Zq::<[< Zq $L ConfigImpl >]<$($Q,)*>, $L>::modulus();
                for x in [0, 1, $($Q/2, $Q-1,)*] {
                    let x_in = Zq::<[< Zq $L ConfigImpl >]<$($Q,)*>, $L>::try_from(x).unwrap();
                    let mut bytes = Vec::new();
                    x_in.serialize_with_mode(&mut bytes, Compress::Yes).unwrap();
                    let x_out = Zq::<[< Zq $L ConfigImpl >]<$($Q,)*>, $L>::deserialize_with_mode(&*bytes, Compress::Yes, Validate::Yes).unwrap();
                    assert_eq!(x_in, x_out);
                }
            }
        }
    }
}

macro_rules! test_fq {
    ($($q:expr,)*) => {
        paste::expr! {
            $(
                #[test]
                fn [< test_f $q >]() {
                    let biguint = BigUint::from(1u64);
                    let x = BigInt::<1>::try_from(&biguint % &BigUint::from($q as u64)).unwrap();
                    let f = Fq::<$q>::new(x);
                    let y = f.into_bigint();
                    assert_eq!(x, y);
                    assert!(f.is_one());
                    assert!(!f.is_zero());
                }
            )*
        }
    };
}

#[cfg(test)]
mod test {
    use super::*;

    // Some primes, big and small. In particular, this tests that the implementation does not rely on any special structure of the prime, nor on the primes being specified in any particular order.
    const Q1: u64 = 3;
    const Q2: u64 = (1u64 << 31) - 1;
    const Q3: u64 = ((1u128 << 64) - 59) as u64;
    const Q4: u64 = ((1u128 << 61) - 1) as u64;
    const Q5: u64 = 7;
    const Q6: u64 = (1u64 << 19) - 1;
    const Q7: u64 = (1u64 << 13) - 1;
    const Q8: u64 = 27644437;
    const Q9: u64 = 200560490131;
    const Q10: u64 = 87178291199;

    test_zq_config_impl!(1, Q1);
    test_zq_config_impl!(2, Q1, Q2);
    test_zq_config_impl!(3, Q1, Q2, Q3);
    test_zq_config_impl!(4, Q1, Q2, Q3, Q4);
    test_zq_config_impl!(5, Q1, Q2, Q3, Q4, Q5);
    // test_zq_config_impl!(6, Q1, Q2, Q3, Q4, Q5, Q6);
    // test_zq_config_impl!(7, Q1, Q2, Q3, Q4, Q5, Q6, Q7);
    // test_zq_config_impl!(8, Q1, Q2, Q3, Q4, Q5, Q6, Q7, Q8);
    // test_zq_config_impl!(9, Q1, Q2, Q3, Q4, Q5, Q6, Q7, Q8, Q9);
    // test_zq_config_impl!(10, Q1, Q2, Q3, Q4, Q5, Q6, Q7, Q8, Q9, Q10);
}
