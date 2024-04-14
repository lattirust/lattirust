#![allow(non_snake_case)]

use std::fmt::{Debug, Display};
use std::hash::Hash;
use std::iter::Sum;
use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

use ark_ff::{BigInt, BitIteratorBE, BitIteratorLE, Fp, Fp64, MontBackend, MontConfig};
use ark_std::UniformRand;
use num_traits::{One, Zero};
use serde::{self, Deserialize, Serialize};

use crate::nimue::serialization::{FromBytes, ToBytes};
use crate::traits::{FromRandomBytes, Modulus};

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
// + Ord
+ Neg<Output=Self>
+ UniformRand
//+ Zeroize
+ Sized
+ Hash
// + CanonicalSerialize
// + CanonicalSerializeWithFlags
// + CanonicalDeserialize
// + CanonicalDeserializeWithFlags
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
+ Sum<Self>
// Differs from arkworks
+ Serialize
+ for<'a> Deserialize<'a>
+ FromRandomBytes<Self>
+ FromBytes
+ ToBytes
+ Modulus
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
    const TWO_ADIC_ROOT_OF_UNITY: Fp<MontBackend<Self, 1>, 1> = todo!(); // TODO
}

pub type Fq<const Q: u64> = Fp64<MontBackend<FqConfig<Q>, 1>>;
pub type Fq2<const Q: u64> = Fp64<MontBackend<FqConfig<Q>, 2>>;

impl<const Q: u64> const Modulus for Fq<Q> {
    fn modulus() -> u64 {
        Q
    }
}

pub const fn const_fq_from<const Q: u64>(val: u64) -> Fq<Q> {
    Fq::new(BigInt::<1> { 0: [val] })
}