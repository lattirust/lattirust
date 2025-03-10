#![allow(non_snake_case)]

use std::fmt::{Debug, Display};
use std::hash::Hash;
use std::iter::{Product, Sum};
use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

use ark_ff::{AdditiveGroup, BitIteratorBE, BitIteratorLE, Field};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::UniformRand;
use num_traits::{One, Zero};

pub use ntt::NttRing;
// Exports
pub use poly_ring::PolyRing;
pub use pow2_cyclotomic_poly_ring::Pow2CyclotomicPolyRing;
pub use pow2_cyclotomic_poly_ring_ntt::Pow2CyclotomicPolyRingNTT;
pub use z_2::*;
pub use z_2_128::*;
pub use z_2_64::*;
pub use z_q::*;

use crate::nimue::serialization::{FromBytes, ToBytes};
use crate::traits::{FromRandomBytes, Modulus, WithL2Norm, WithLinfNorm};

pub mod ntt;
mod poly_ring;
pub(crate) mod pow2_cyclotomic_poly_ring;
pub(crate) mod pow2_cyclotomic_poly_ring_ntt;
pub mod representatives;
pub mod util;
mod z_2;
mod z_2_128;
mod z_2_64;
mod z_q;
pub(crate) mod z_q_signed_representative;

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
+ CanonicalSerialize
// + CanonicalSerializeWithFlags
+ CanonicalDeserialize
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
+ for<'a> Sum<&'a Self>
+ Product<Self>
+ for<'a> Product<&'a Self>
+ Sum<Self>
+ TryFrom<u128, Error: Debug>
+ TryFrom<u64, Error: Debug>
+ TryFrom<u32, Error: Debug>
+ TryFrom<u16, Error: Debug>
+ TryFrom<u8, Error: Debug>
+ From<bool>
// Differs from arkworks
+ FromRandomBytes<Self>
+ FromBytes
+ ToBytes
+ Modulus
+ WithL2Norm
+ WithLinfNorm
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

     fn inverse(&self) -> Option<Self>;
}

impl<T> Ring for T
where
    T: Field + AdditiveGroup + Modulus + FromRandomBytes<T> + WithLinfNorm + WithL2Norm,
{
    const ZERO: Self = <Self as AdditiveGroup>::ZERO;
    const ONE: Self = <Self as Field>::ONE;

    fn inverse(&self) -> Option<Self> {
        <Self as Field>::inverse(self)
    }
}

#[macro_export]
macro_rules! test_ring {
    ($T:ty, $N:expr) => {
        test_associative_addition!($T, $N);
        test_associative_multiplication!($T, $N);
        test_distributive!($T, $N);
        test_identity_addition!($T, $N);
        test_identity_multiplication!($T, $N);
        test_inverse_addition!($T, $N);
        test_inverse_multiplication!($T, $N);
        test_canonical_serialize_deserialize_uncompressed!($T, $N);
        test_canonical_serialize_deserialize_compressed!($T, $N);
    };
}

#[macro_export]
macro_rules! test_distributive {
    ($T:ty, $N:expr) => {
        #[test]
        fn test_distributive() {
            let rng = &mut ark_std::test_rng();
            for _ in 0..$N {
                let a = <$T as UniformRand>::rand(rng);
                let b = <$T as UniformRand>::rand(rng);
                let c = <$T as UniformRand>::rand(rng);

                assert_eq!(a * (b + c), a * b + a * c);
                assert_eq!((a + b) * c, a * c + b * c);
            }
        }
    };
}

#[macro_export]
macro_rules! test_associative_addition {
    ($T:ty, $N:expr) => {
        #[test]
        fn test_associative_addition() {
            let rng = &mut ark_std::test_rng();
            for _ in 0..$N {
                let a = <$T as UniformRand>::rand(rng);
                let b = <$T as UniformRand>::rand(rng);
                let c = <$T as UniformRand>::rand(rng);

                assert_eq!(a + (b + c), (a + b) + c);
            }
        }
    };
}

#[macro_export]
macro_rules! test_associative_multiplication {
    ($T:ty, $N:expr) => {
        #[test]
        fn test_associative_multiplication() {
            let rng = &mut ark_std::test_rng();
            for _ in 0..$N {
                let a = <$T as UniformRand>::rand(rng);
                let b = <$T as UniformRand>::rand(rng);
                let c = <$T as UniformRand>::rand(rng);

                assert_eq!(a * (b * c), (a * b) * c);
            }
        }
    };
}

#[macro_export]
macro_rules! test_identity_addition {
    ($T:ty, $N:expr) => {
        #[test]
        fn test_identity_addition() {
            let rng = &mut ark_std::test_rng();
            for _ in 0..$N {
                let a = <$T as UniformRand>::rand(rng);

                assert_eq!(a + <$T as Ring>::ZERO, a);
                assert_eq!(<$T as Ring>::ZERO + a, a);
            }
        }
    };
}

#[macro_export]
macro_rules! test_identity_multiplication {
    ($T:ty, $N:expr) => {
        #[test]
        fn test_identity_multiplication() {
            let rng = &mut ark_std::test_rng();
            for _ in 0..$N {
                let a = <$T as UniformRand>::rand(rng);

                assert_eq!(a * <$T as Ring>::ONE, a);
                assert_eq!(<$T as Ring>::ONE * a, a);
            }
        }
    };
}

#[macro_export]
macro_rules! test_inverse_addition {
    ($T:ty, $N:expr) => {
        #[test]
        fn test_inverse_addition() {
            let rng = &mut ark_std::test_rng();
            for _ in 0..$N {
                let a = <$T as UniformRand>::rand(rng);

                assert_eq!(a + -a, <$T as Ring>::ZERO);
                assert_eq!(-a + a, <$T as Ring>::ZERO);
            }
        }
    };
}

#[macro_export]
macro_rules! test_inverse_multiplication {
    ($T:ty, $N:expr) => {
        #[test]
        fn test_inverse_multiplication() {
            let rng = &mut ark_std::test_rng();
            assert_eq!(
                <$T as Ring>::inverse(&<$T as Ring>::ONE).unwrap(),
                <$T as Ring>::ONE
            );
            assert_eq!(<$T as Ring>::inverse(&<$T as Ring>::ZERO), None);
            for _ in 0..$N {
                let a = <$T as UniformRand>::rand(rng);
                let inv = <$T as Ring>::inverse(&a);
                if inv.is_some() {
                    assert_eq!(a * inv.unwrap(), <$T as Ring>::ONE);
                    assert_eq!(inv.unwrap() * a, <$T as Ring>::ONE);
                }
            }
        }
    };
}

#[macro_export]
macro_rules! test_canonical_serialize_deserialize_compressed {
    ($T:ty, $N:expr) => {
        #[test]
        fn test_canonical_serialize_deserialize_compressed() {
            let rng = &mut ark_std::test_rng();

            for _ in 0..$N {
                let a = <$T as UniformRand>::rand(rng);
                let mut bytes = Vec::new();
                a.serialize_with_mode(&mut bytes, Compress::Yes).unwrap();
                let a2 = <$T as CanonicalDeserialize>::deserialize_with_mode(
                    &*bytes,
                    Compress::Yes,
                    Validate::Yes,
                )
                .unwrap();
                assert_eq!(a, a2);
            }
        }
    };
}

#[macro_export]
macro_rules! test_canonical_serialize_deserialize_uncompressed {
    ($T:ty, $N:expr) => {
        #[test]
        fn test_canonical_serialize_deserialize_uncompressed() {
            let rng = &mut ark_std::test_rng();

            for _ in 0..$N {
                let a = <$T as UniformRand>::rand(rng);
                let mut bytes = Vec::new();
                a.serialize_with_mode(&mut bytes, Compress::No).unwrap();
                let a2 = <$T as CanonicalDeserialize>::deserialize_with_mode(
                    &*bytes,
                    Compress::No,
                    Validate::Yes,
                )
                .unwrap();
                assert_eq!(a, a2);
            }
        }
    };
}
