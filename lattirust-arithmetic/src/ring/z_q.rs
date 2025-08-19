#![allow(warnings)] // TODO: remove
#![allow(long_running_const_eval)]

use std::array;
use std::convert::Into;
use std::fmt::Debug;
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
use num_traits::{One, Signed, ToPrimitive, Zero};
use rounded_div::RoundedDiv;
use zeroize::Zeroize;

use crate::decomposition::DecompositionFriendlySignedRepresentative;
use crate::impl_try_from_primitive_type;
use crate::ring::ntt::{
    const_fq_from, const_pow_mod, generator, is_primitive_root_of_unity, two_adic_root_of_unity,
    Ntt,
};
use crate::ring::representatives::{SignedRepresentative, WithSignedRepresentative};
use crate::ring::{NttRing, Ring};
use crate::traits::{FromRandomBytes, Modulus, WithLinfNorm};

use crate::ring::f_p::{fq_zero, Fq, FqConfig};

const fn to_bigint_assert_odd_prime<const Q: u64>() -> BigInt<1> {
    assert!(
        Q > 2 && const_primes::is_prime(Q),
        "You tried to instantiate an FqConfig<Q> with either Q = 2 or Q not a prime"
    );
    BigInt::<1>([Q])
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

impl<C: ZqConfig<L>, const L: usize> FromRandomBytes<Zq<C, L>> for Zq<C, L> {
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

impl<C: ZqConfig<L>, const L: usize> From<Zq<C, L>> for BigUint {
    fn from(value: Zq<C, L>) -> Self {
        C::into_bigint(value).into()
    }
}

impl<M: Modulus, C: ZqConfig<L>, const L: usize> From<SignedRepresentative<M>> for Zq<C, L> {
    fn from(value: SignedRepresentative<M>) -> Self {
        if Self::modulus() < M::modulus() {
            panic!("Cannot convert signed representative to {}, as the signed representative modulus {} is larger than the modulus {}.", M::modulus(), std::any::type_name::<Self>(),  Self::modulus());
        }
        let mut bigint = value.0;
        if bigint.is_negative() {
            let modulus: num_bigint::BigInt = Self::modulus().into();
            bigint = bigint + modulus;
        }
        let biguint: BigUint = bigint.to_biguint().unwrap();
        Self::try_from(biguint).unwrap()
    }
}



// TODO: more efficient implementation based on CRT decomposition and signed Fp implementation?
impl<C: ZqConfig<L>, const L: usize> WithSignedRepresentative for Zq<C, L> {
    type SignedRepresentative = SignedRepresentative<Zq<C, L>>;

    fn as_signed_representative(&self) -> Self::SignedRepresentative {
        (*self).clone().into()
    }

    fn signed_representative_to_bigint(repr: &Self::SignedRepresentative) ->num_bigint::BigInt {
        repr.0.clone()  
    }
    
    fn signed_representative_from_bigint(value: num_bigint::BigInt) -> Option<Self::SignedRepresentative> {
        Some(SignedRepresentative::new(value))
    }
}

const fn all_distinct<const N: usize>(xs: [u64; N]) -> bool {
    let mut i = 0;
    while i < N {
        let mut j = i + 1;
        while j < N {
            if xs[i] == xs[j] {
                return false;
            }
            j += 1;
        }
        i += 1;
    }
    return true;
}

const fn id_assert_all_distinct<const L: usize>(xs: [u64; L]) -> [u64; L] {
    assert!(
        all_distinct(xs),
        "You tried to instantiate an ZqConfigLImpl<Q1, ..., QL> with Qi not distinct."
    );
    xs
}

/// To be called as `zq_config_impl!(L, 0, 1, ..., L-1);` for some non-zero number of limbs `L`
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
                const MODULI: [u64; $L] = id_assert_all_distinct([$([< Q $i >],)*]); // Fails at compile-time if the moduli are not distinct
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
                    Zq::<Self, $L>((
                        $(<Fq::<[< Q $i >]> as Field>::sum_of_products(
                            &a.map(|x| x.0.$i),
                            &b.map(|x| x.0.$i)
                        ),)*
                    ))
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

            // TODO: this might not be the more efficient implementation, we're using an array-of-structs, and not doing NTTs/INTTs in-place.
            impl<$(const [< Q $i >]: u64,)* const N: usize> Ntt<N> for Zq<[< Zq $L ConfigImpl >]<$([< Q $i >],)*>, $L>
                where $(Fq<[< Q $i >]>: Ntt<N>,)*
            {
                fn ntt_inplace(coeffs: &mut [Self; N]) {
                    $(
                    let mut [< arr $i >] = coeffs.map(|x| x.0 .$i);
                    Fq::<[< Q $i >]>::ntt_inplace(&mut [< arr $i >]);
                    )*
                    *coeffs = core::array::from_fn(|j| Self(($([< arr $i >][j],)*)))
                }

                fn intt_inplace(evals: &mut [Self; N]) {
                    $(
                    let mut [< arr $i >] = evals.map(|x| x.0 .$i);
                    Fq::<[< Q $i >]>::intt_inplace(&mut [< arr $i >]);
                    )*
                    *evals = core::array::from_fn(|j| Self(($([< arr $i >][j],)*)))
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
// TODO: arkworks only automatically derives CanonicalSerialize/CanonicalDeserialize for up to 5-tuples.
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

#[cfg(test)]
mod test {
    use super::*;
    use crate::*;

    // Some primes, big and small. In particular, this tests that the implementation does not rely on any special structure of the prime, nor on the primes being specified in any particular order.
    const Q1: u64 = (1 << 31) - (1 << 27) + 1; // BabyBear prime, NTT-friendly up to N=2^26
    const Q2: u64 = 274177; // LaBRADOR modulus factor 1, NTT-friendly up to N=2^7
    const Q3: u64 = 67280421310721; // LaBRADOR modulus factor 2, NTT-friendly up to N=2^7
    const Q4: u64 = ((1u128 << 64) - (1u128 << 32) + 1) as u64; // Goldilocks prime, NTT-friendly up to N=2^31
    const Q5: u64 = 3; // Not NTT-friendly
    const Q6: u64 = (1 << 31) - 1; // Mersenne prime, not NTT-friendly
    const Q7: u64 = (1 << 13) - 1; // Not NTT-friendly
    const Q8: u64 = 27644437; // Not NTT-friendly
    const Q9: u64 = 200560490131; // Not NTT-friendly
    const Q10: u64 = 7; // Not NTT-friendly

    #[cfg(test)]
    mod test_z1 {
        use super::*;
        type Z1 = Zq1<Q1>;

        test_ring!(Z1, 100);
        test_zq_config_impl!(1, Q1);
        test_ntt_intt!(Z1, Z1, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192);
        test_ntt_add!(Z1, Z1, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192);
        test_ntt_mul!(Z1, Z1, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192);
    }

    #[cfg(test)]
    mod test_z2 {
        use super::*;
        type Z2 = Zq2<Q1, Q2>;

        test_ring!(Z2, 100);
        test_zq_config_impl!(2, Q1, Q2);
        test_ntt_intt!(Z2, Z2, 32, 64, 128);
        test_ntt_add!(Z2, Z2, 32, 64, 128);
        test_ntt_mul!(Z2, Z2, 32, 64, 128);
    }

    #[cfg(test)]
    mod test_z3 {
        use super::*;
        type Z3 = Zq3<Q1, Q2, Q3>;

        test_ring!(Z3, 100);
        test_zq_config_impl!(3, Q1, Q2, Q3);
        test_ntt_intt!(Z3, Z3, 32, 64, 128);
        test_ntt_add!(Z3, Z3, 32, 64, 128);
        test_ntt_mul!(Z3, Z3, 32, 64, 128);
    }

    #[cfg(test)]
    mod test_z4 {
        use super::*;
        type Z4 = Zq4<Q1, Q2, Q3, Q4>;

        test_ring!(Z4, 100);
        test_zq_config_impl!(4, Q1, Q2, Q3, Q4);
        test_ntt_intt!(Z4, Z4, 32, 64, 128);
        test_ntt_add!(Z4, Z4, 32, 64, 128);
        test_ntt_mul!(Z4, Z4, 32, 64, 128);
    }

    #[cfg(test)]
    mod test_z5 {
        use super::*;
        type Z5 = Zq5<Q1, Q2, Q3, Q4, Q5>;

        test_ring!(Z5, 100);
        test_zq_config_impl!(5, Q1, Q2, Q3, Q4, Q5);
    }
}
