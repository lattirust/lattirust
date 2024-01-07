use std::cmp::Ordering;
use std::fmt::{Debug, Display, Formatter};
use std::hash::{Hash, Hasher};
use std::io::{Read, Write};
use std::iter::{Product, Sum};
use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

use ark_ff::{BigInt, BitIteratorBE, BitIteratorLE, Field, Fp, Fp64, MontBackend, MontConfig, PrimeField};
use ark_serialize::{CanonicalDeserialize, CanonicalDeserializeWithFlags, CanonicalSerialize, CanonicalSerializeHashExt, CanonicalSerializeWithFlags, Compress, Flags, SerializationError, Valid, Validate};
use ark_std::UniformRand;
use delegate_attr::delegate;
use derive_more::{Deref, DerefMut, Mul, MulAssign, Neg};
use derive_more::Into;
use num_traits::{One, Zero};
use rand::Rng;
use serde::{self, Deserialize, Deserializer, Serialize, Serializer};
use zeroize::Zeroize;

use crate::lattice_arithmetic::traits::{FromRandomBytes, IntegerDiv, Modulus, WithLog2};
use crate::nimue::serialization::FromBytes;

pub trait Ring:
'static
+ Copy
+ Clone
+ Debug
+ Display
+ Default
+ Send
+ Sync
+ Eq
+ Zero
+ One
+ Ord
+ Neg<Output=Self>
+ UniformRand
//+ Zeroize
+ Sized
+ Hash
+ CanonicalSerialize
+ CanonicalSerializeWithFlags
+ CanonicalDeserialize
+ CanonicalDeserializeWithFlags
+ Add<Self, Output=Self>
+ Sub<Self, Output=Self>
+ Mul<Self, Output=Self>
+ AddAssign<Self>
+ SubAssign<Self>
+ MulAssign<Self>
+ for<'a> Add<&'a Self, Output=Self>
+ for<'a> Sub<&'a Self, Output=Self>
+ for<'a> Mul<&'a Self, Output=Self>
+ for<'a> AddAssign<&'a Self>
+ for<'a> SubAssign<&'a Self>
+ for<'a> MulAssign<&'a Self>
+ for<'a> Add<&'a mut Self, Output=Self>
+ for<'a> Sub<&'a mut Self, Output=Self>
+ for<'a> Mul<&'a mut Self, Output=Self>
+ for<'a> AddAssign<&'a mut Self>
+ for<'a> SubAssign<&'a mut Self>
+ for<'a> MulAssign<&'a mut Self>
+ core::iter::Sum<Self>
//+ for<'a> core::iter::Sum<&'a Self>
+ core::iter::Product<Self>
+ for<'a> core::iter::Product<&'a Self>
+ From<u128>
+ From<u64>
+ From<u32>
+ From<u16>
+ From<u8>
+ From<bool>
// Differs from arkworks
+ Serialize
+ for<'a> Deserialize<'a>
+ FromRandomBytes<Self>
+ FromBytes
{
    /// The additive identity of the ring.
    const ZERO: Self;
    /// The multiplicative identity of the ring.
    const ONE: Self;

    /// Returns `sum([a_i * b_i])`.
    #[inline]
    fn sum_of_products<const T: usize>(a: &[Self; T], b: &[Self; T]) -> Self {
        let mut sum = Self::zero();
        for i in 0..a.len() {
            sum += a[i] * b[i];
        }
        sum
    }

    fn square_in_place(&mut self) -> &mut Self {
        *self *= *self;
        self
    }


    /// Returns `self^exp`, where `exp` is an integer represented with `u64` limbs,
    /// least significant limb first.
    #[must_use]
    fn pow<S: AsRef<[u64]>>(&self, exp: S) -> Self {
        let mut res = Self::one();

        for i in BitIteratorBE::without_leading_zeros(exp) {
            res.square_in_place();

            if i {
                res *= self;
            }
        }
        res
    }


    /// Exponentiates a field element `f` by a number represented with `u64`
    /// limbs, using a precomputed table containing as many powers of 2 of
    /// `f` as the 1 + the floor of log2 of the exponent `exp`, starting
    /// from the 1st power. That is, `powers_of_2` should equal `&[p, p^2,
    /// p^4, ..., p^(2^n)]` when `exp` has at most `n` bits.
    ///
    /// This returns `None` when a power is missing from the table.
    #[inline]
    fn pow_with_table<S: AsRef<[u64]>>(powers_of_2: &[Self], exp: S) -> Option<Self> {
        let mut res = Self::one();
        for (pow, bit) in BitIteratorLE::without_trailing_zeros(exp).enumerate() {
            if bit {
                res *= powers_of_2.get(pow)?;
            }
        }
        Some(res)
    }
}

pub struct FqConfig<const Q: u64> {}

impl<const Q: u64> MontConfig<1> for FqConfig<Q> {
    const MODULUS: BigInt<1> = BigInt::<1> { 0: [Q] };
    const GENERATOR: Fp<MontBackend<Self, 1>, 1> = Fp::new(BigInt::<1> { 0: [2u64] });
    // TODO: check if this needed/makes sense
    const TWO_ADIC_ROOT_OF_UNITY: Fp<MontBackend<Self, 1>, 1> = Fp::new(BigInt::zero()); // TODO
}

type Fq<const Q: u64> = (Fp64<MontBackend<FqConfig<Q>, 1>>);

#[derive(Deref, DerefMut, Clone, Debug, Eq, PartialEq, Into, Mul, MulAssign, Neg, Zeroize)]
#[derive(Serialize, Deserialize)]
pub struct Zq<const Q: u64>(
    #[serde(with = "crate::nimue::serialization::ser")]
    Fq<Q>
);

impl<const Q: u64> Zq<Q> {
    pub fn Q() -> Self { Self::from(Q) }
}

impl<const Q: u64> From<Fq<Q>> for Zq<Q> {
    fn from(value: Fq<Q>) -> Self {
        Self(value)
    }
}

impl<const Q: u64> From<Zq<Q>> for u64 {
    fn from(value: Zq<Q>) -> Self {
        let res = value.0.into_bigint().0[0];
        assert!(res < Q);
        res
    }
}

impl<const Q: u64> From<Zq<Q>> for i64 {
    fn from(value: Zq<Q>) -> Self {
        let res_u64 = u64::from(value);
        if res_u64 > (Q / 2) {
            -(Q as i64) + (res_u64 as i64)
        } else {
            res_u64 as i64
        }
    }
}

impl<const Q: u64> From<Zq<Q>> for i128 {
    fn from(value: Zq<Q>) -> Self {
        i64::from(value) as i128
    }
}


#[delegate(self.0)]
impl<const Q: u64> Zq<Q> {}

impl<const Q: u64> Copy for Zq<Q> {}

impl<const Q: u64> Display for Zq<Q> {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        // std::fmt::Display::fmt(&u64::from(*self), f)
        std::fmt::Display::fmt(&self.0, f)
    }
}

impl<const Q: u64> Default for Zq<Q> {
    fn default() -> Self {
        Fq::default().into()
    }
}


impl<const Q: u64> Zero for Zq<Q> {
    fn zero() -> Self { Self::ZERO }
    fn is_zero(&self) -> bool { self.0.is_zero() }
}

impl<const Q: u64> Add<Self> for Zq<Q> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        self.0.add(rhs.0).into()
    }
}

impl<const Q: u64> One for Zq<Q> {
    fn one() -> Self {
        Self::ONE
    }
}

impl<const Q: u64> Mul<Self> for Zq<Q> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        self.0.mul(rhs.0).into()
    }
}

impl<const Q: u64> Ord for Zq<Q> {
    fn cmp(&self, other: &Self) -> Ordering {
        self.0.cmp(other)
    }
}

impl<const Q: u64> PartialOrd<Self> for Zq<Q> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        self.0.partial_cmp(&other.0)
    }
}

impl<const Q: u64> UniformRand for Zq<Q> {
    fn rand<R: Rng + ?Sized>(rng: &mut R) -> Self {
        Fq::rand(rng).into()
    }
}

impl<const Q: u64> Hash for Zq<Q> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        Hash::hash(&self.0, state)
    }
}

impl<const Q: u64> CanonicalSerialize for Zq<Q> {
    fn serialize_with_mode<W: Write>(&self, writer: W, compress: Compress) -> Result<(), SerializationError> {
        self.0.serialize_with_mode(writer, compress)
    }

    fn serialized_size(&self, compress: Compress) -> usize {
        self.0.serialized_size(compress)
    }
}

impl<const Q: u64> CanonicalSerializeWithFlags for Zq<Q> {
    fn serialize_with_flags<W: Write, F: Flags>(&self, writer: W, flags: F) -> Result<(), SerializationError> {
        self.0.serialize_with_flags(writer, flags)
    }

    fn serialized_size_with_flags<F: Flags>(&self) -> usize {
        self.0.serialized_size_with_flags::<F>()
    }
}

impl<const Q: u64> CanonicalDeserialize for Zq<Q> {
    fn deserialize_with_mode<R: Read>(reader: R, compress: Compress, validate: Validate) -> Result<Self, SerializationError> {
        Fq::deserialize_with_mode(reader, compress, validate).map(|fq| fq.into())
    }
}

impl<const Q: u64> Valid for Zq<Q> {
    fn check(&self) -> Result<(), SerializationError> {
        self.0.check()
    }
}

impl<const Q: u64> CanonicalDeserializeWithFlags for Zq<Q> {
    fn deserialize_with_flags<R: Read, F: Flags>(reader: R) -> Result<(Self, F), SerializationError> {
        Fq::deserialize_with_flags(reader).map(|(fq, flags)| (fq.into(), flags))
    }
}

impl<const Q: u64> Sub<Self> for Zq<Q> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        self.0.sub(rhs.0).into()
    }
}

impl<const Q: u64> AddAssign<Self> for Zq<Q> {
    fn add_assign(&mut self, rhs: Self) {
        self.0.add_assign(rhs.0)
    }
}

impl<const Q: u64> SubAssign<Self> for Zq<Q> {
    fn sub_assign(&mut self, rhs: Self) { self.0.sub_assign(rhs.0) }
}

impl<const Q: u64> MulAssign<Self> for Zq<Q> {
    fn mul_assign(&mut self, rhs: Self) { self.0.mul_assign(rhs.0) }
}

impl<'a, const Q: u64> Add<&'a Self> for Zq<Q> {
    type Output = Self;

    fn add(self, rhs: &'a Self) -> Self::Output {
        self.0.add(rhs.0).into()
    }
}

impl<'a, const Q: u64> Sub<&'a Self> for Zq<Q> {
    type Output = Self;

    fn sub(self, rhs: &'a Self) -> Self::Output {
        self.0.sub(rhs.0).into()
    }
}

impl<'a, const Q: u64> Mul<&'a Self> for Zq<Q> {
    type Output = Self;

    fn mul(self, rhs: &'a Self) -> Self::Output {
        self.0.mul(rhs.0).into()
    }
}

impl<'a, const Q: u64> AddAssign<&'a Self> for Zq<Q> {
    fn add_assign(&mut self, rhs: &'a Self) {
        self.0.add_assign(rhs.0)
    }
}

impl<'a, const Q: u64> SubAssign<&'a Self> for Zq<Q> {
    fn sub_assign(&mut self, rhs: &'a Self) { self.0.sub_assign(rhs.0) }
}

impl<'a, const Q: u64> MulAssign<&'a Self> for Zq<Q> {
    fn mul_assign(&mut self, rhs: &'a Self) {
        self.0.mul_assign(rhs.0)
    }
}

impl<'a, const Q: u64> Add<&'a mut Self> for Zq<Q> {
    type Output = Self;

    fn add(self, rhs: &'a mut Self) -> Self::Output {
        self.0.add(rhs.0).into()
    }
}

impl<'a, const Q: u64> Sub<&'a mut Self> for Zq<Q> {
    type Output = Self;

    fn sub(self, rhs: &'a mut Self) -> Self::Output {
        self.0.sub(rhs.0).into()
    }
}

impl<'a, const Q: u64> Mul<&'a mut Self> for Zq<Q> {
    type Output = Self;

    fn mul(self, rhs: &'a mut Self) -> Self::Output {
        self.0.mul(rhs.0).into()
    }
}

impl<'a, const Q: u64> AddAssign<&'a mut Self> for Zq<Q> {
    fn add_assign(&mut self, rhs: &'a mut Self) {
        self.0.add_assign(rhs.0)
    }
}

impl<'a, const Q: u64> SubAssign<&'a mut Self> for Zq<Q> {
    fn sub_assign(&mut self, rhs: &'a mut Self) {
        self.0.sub_assign(rhs.0)
    }
}

impl<'a, const Q: u64> MulAssign<&'a mut Self> for Zq<Q> {
    fn mul_assign(&mut self, rhs: &'a mut Self) {
        self.0.mul_assign(rhs.0)
    }
}

impl<'a, const Q: u64> Sum<Self> for Zq<Q> {
    fn sum<I: Iterator<Item=Self>>(iter: I) -> Self { Self(Sum::sum(iter.map(|x| x.0))) }
}

impl<'a, const Q: u64> Sum<&'a Self> for Zq<Q> {
    fn sum<I: Iterator<Item=&'a Self>>(iter: I) -> Self {
        todo!()
    }
}

impl<const Q: u64> Product<Self> for Zq<Q> {
    fn product<I: Iterator<Item=Self>>(iter: I) -> Self {
        todo!()
    }
}

impl<'a, const Q: u64> Product<&'a Self> for Zq<Q> {
    fn product<I: Iterator<Item=&'a Self>>(iter: I) -> Self {
        todo!()
    }
}

impl<const Q: u64> From<u128> for Zq<Q> {
    fn from(value: u128) -> Self { Fq::from(value).into() }
}

impl<const Q: u64> From<u64> for Zq<Q> {
    fn from(value: u64) -> Self { Fq::from(value).into() }
}

impl<const Q: u64> From<u32> for Zq<Q> {
    fn from(value: u32) -> Self { Fq::from(value).into() }
}

impl<const Q: u64> From<u16> for Zq<Q> {
    fn from(value: u16) -> Self { Fq::from(value).into() }
}

impl<const Q: u64> From<u8> for Zq<Q> {
    fn from(value: u8) -> Self { Fq::from(value).into() }
}

impl<const Q: u64> From<bool> for Zq<Q> {
    fn from(value: bool) -> Self { Fq::from(value).into() }
}

impl<const Q: u64> Ring for Zq<Q> {
    const ZERO: Self = Zq(Fq::new(BigInt::<1> { 0: [0u64] }));
    const ONE: Self = Zq(Fq::new(BigInt::<1> { 0: [1u64] }));
}

impl<const Q: u64> FromRandomBytes<Self> for Zq<Q> {
    fn byte_size() -> usize {
        Fq::<Q>::MODULUS_BIT_SIZE as usize / 8
    }

    fn from_random_bytes(bytes: &[u8]) -> Option<Self> {
        Fq::<Q>::from_random_bytes(bytes).map(|fq| fq.into())
    }
}

impl<const Q: u64> IntegerDiv for Zq<Q> {
    fn integer_div(&self, rhs: Self) -> Self {
        Zq::from((u64::from(*self)).div_euclid(u64::from(rhs)))
    }

    fn div_round(&self, rhs: Self) -> Self {
        Zq::from(((u64::from(*self) as f64) / (u64::from(rhs) as f64)).round() as u64)
    }
}

impl<const Q: u64> WithLog2 for Zq<Q> {
    fn log2(&self) -> f64 {
        f64::log2(u64::from(*self) as f64)
    }

    fn log2_q() -> f64 {
        f64::log2(Q as f64)
    }
}

impl<const Q: u64> Modulus for Zq<Q> {
    fn modulus() -> u64 {
        Q
    }
}


mod tests {
    use num_traits::{One, Zero};
    use serde::Serialize;

    use crate::lattice_arithmetic::ring::Ring;

    const Q: u64 = 11;

    type R = super::Zq<Q>;

    #[test]
    fn test_zero() {
        let r = R::from(0u64);
        assert_eq!(r, R::ZERO);
        assert_eq!(r, R::zero());
        assert!(r.is_zero());
    }

    #[test]
    fn test_one() {
        let r = R::from(1u64);
        assert_eq!(r, R::ONE);
        assert_eq!(r, R::one());
        assert!(r.is_one());
    }

    #[test]
    fn test_from() {
        for v in 0..Q {
            assert_eq!(v, u64::from(R::from(v)));
        }
    }

    #[test]
    fn test_ser() {
        for v in 0..Q {
            let r = R::from(v as u64);
            let ser = bincode::serialize(&r);
            assert!(ser.is_ok());
        }
    }

    #[test]
    fn test_des() {
        for v in 0..Q {
            let r = R::from(v as u64);
            let s = bincode::serialize(&r);
            let der = bincode::deserialize(&s.unwrap());
            assert!(der.is_ok());
            let r_ = der.unwrap();
            assert_eq!(r, r_, "r = {} != des(ser(r)) = {}", r, r_);
        }
    }

    #[test]
    fn test_from_random_bytes() {
        todo!("ensure from_random_bytes() returns a valid result when called on byte arrays of size byte_size()");
    }
}