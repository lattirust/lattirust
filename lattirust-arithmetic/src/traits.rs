use std::array;

use num_bigint::BigUint;
use num_traits::{Signed, ToPrimitive};

use crate::linear_algebra::Vector;
use crate::ring::representatives::WithSignedRepresentative;
use crate::ring::PolyRing;

pub trait WithL2Norm {
    fn l2_norm(&self) -> f64 {
        self.l2_norm_squared().to_f64().unwrap().sqrt()
    }
    fn l2_norm_squared(&self) -> BigUint;
}

impl<R: WithSignedRepresentative> WithL2Norm for R
where
    R::SignedRepresentative: Signed + Into<num_bigint::BigInt>,
{
    fn l2_norm_squared(&self) -> BigUint {
        let signed: num_bigint::BigInt = self.as_signed_representative().into();
        BigUint::try_from(&signed * &signed).unwrap()
    }
}

impl<R: WithL2Norm> WithL2Norm for [R] {
    fn l2_norm_squared(&self) -> BigUint {
        self.iter().map(|x| x.l2_norm_squared()).sum()
    }
}

impl<R: WithL2Norm> WithL2Norm for Vec<R> {
    fn l2_norm_squared(&self) -> BigUint {
        self.as_slice().l2_norm_squared()
    }
}

pub trait WithLinfNorm {
    fn linf_norm(&self) -> BigUint;
}

impl<R: WithSignedRepresentative> WithLinfNorm for R
where
    R::SignedRepresentative: Signed + Into<num_bigint::BigInt>,
{
    fn linf_norm(&self) -> BigUint {
        let abs: num_bigint::BigInt = self.as_signed_representative().abs().into();
        abs.try_into().unwrap()
    }
}

impl<R: WithLinfNorm> WithLinfNorm for [R] {
    fn linf_norm(&self) -> BigUint {
        self.iter().map(|x| x.linf_norm()).max().unwrap()
    }
}

impl<R: WithLinfNorm> WithLinfNorm for Vec<R> {
    fn linf_norm(&self) -> BigUint {
        self.as_slice().linf_norm()
    }
}

pub trait IntegerDiv<Rhs = Self> {
    /// Divides `self` by `rhs`, returning the quotient.
    fn integer_div(&self, rhs: &Rhs) -> Self;

    fn div_round(&self, rhs: &Rhs) -> Self;
}

pub trait Modulus {
    fn modulus() -> BigUint;
}

/// If `C: FromRandomBytes<T>`, then `C` defines a distribution over `T`, and defines how values `T` can be created from random bytes.
/// This is mostly used to generate verifier challenges from the output of a hash function.
pub trait FromRandomBytes<T> {
    const SECURITY_PARAMETER: usize = 128;

    /// `has_no_biases()` should return false for distributions `C` that are not typically bit-aligned (e.g., for generating a uniformly random element in $Z_p$, where generating a random integer and outputting it $\mod p$ would introduce some bias).
    /// If `has_no_biases() == false` (and rather than using rejection sampling) this implementation requires `SECURITY_PARAMETER / 8` additional (leading) random bytes and throws them away, and allows an implementation to use modular reduction on the remaining bytes.
    /// For a suitably large `SECURITY_PARAMETER`, the bias introduced by this method is negligible.
    /// Do *not* override this to return `true` unless you are absolutely confident that `try_from_random_bytes_inner` implements the intended distribution without bias.
    fn has_no_bias() -> bool {
        false
    }

    /// Returns the minimum number of bytes required to generate an element (including `SECURITY_PARAMETER / 8` bytes if `has_no_bias() == false`)
    fn byte_size() -> usize {
        if Self::has_no_bias() {
            Self::needs_bytes()
        } else {
            Self::needs_bytes() + Self::SECURITY_PARAMETER / 8
        }
    }

    /// Returns the minimum number of bytes required to generate an element
    fn needs_bytes() -> usize;

    /// Returns `Some(t)` for a `t` created from `bytes`, or `None` if there was an error.
    fn try_from_random_bytes(bytes: &[u8]) -> Option<T> {
        if bytes.len() < Self::byte_size() {
            return None;
        }
        if Self::has_no_bias() {
            Self::try_from_random_bytes_inner(&bytes)
        } else {
            Self::try_from_random_bytes_inner(&bytes[(Self::SECURITY_PARAMETER / 8)..])
        }
    }

    /// Returns `Some(t)` for a `t` created from `bytes`, or `None` if there was an error.
    fn try_from_random_bytes_inner(bytes: &[u8]) -> Option<T>;
}

impl<T, const N: usize> FromRandomBytes<[T; N]> for [T; N]
where
    T: FromRandomBytes<T>,
{
    fn needs_bytes() -> usize {
        N * T::needs_bytes() // Use T::needs_byte() rather than T::byte_size() for efficiency: this returns `SECPARAM + N * T::needs_bytes()` instead of `N * (SECPARAM + T::needs_bytes())` in the biased case.
    }

    fn try_from_random_bytes_inner(bytes: &[u8]) -> Option<[T; N]> {
        let chunks: Vec<&[u8]> = bytes.chunks_exact(T::needs_bytes()).collect();
        Some(array::from_fn(|i| {
            T::try_from_random_bytes_inner(chunks[i]).unwrap()
        }))
    }
}

impl<T1> FromRandomBytes<(T1,)> for (T1,)
where
    T1: FromRandomBytes<T1>,
{
    fn needs_bytes() -> usize {
        T1::needs_bytes()
    }

    fn try_from_random_bytes_inner(bytes: &[u8]) -> Option<(T1,)> {
        Some((T1::try_from_random_bytes_inner(bytes)?,))
    }
}

impl<T1, T2> FromRandomBytes<(T1, T2)> for (T1, T2)
where
    T1: FromRandomBytes<T1>,
    T2: FromRandomBytes<T2>,
{
    fn needs_bytes() -> usize {
        T1::needs_bytes() + T2::needs_bytes()
    }

    fn try_from_random_bytes_inner(bytes: &[u8]) -> Option<(T1, T2)> {
        let t1_len = T1::needs_bytes();
        let t1 = T1::try_from_random_bytes_inner(&bytes[..t1_len])?;
        let t2 = T2::try_from_random_bytes_inner(&bytes[t1_len..])?;
        Some((t1, t2))
    }
}

impl<T1, T2, T3> FromRandomBytes<(T1, T2, T3)> for (T1, T2, T3)
where
    T1: FromRandomBytes<T1>,
    T2: FromRandomBytes<T2>,
    T3: FromRandomBytes<T3>,
{
    fn needs_bytes() -> usize {
        T1::needs_bytes() + T2::needs_bytes() + T3::needs_bytes()
    }

    fn try_from_random_bytes_inner(bytes: &[u8]) -> Option<(T1, T2, T3)> {
        let t1_len = T1::needs_bytes();
        let t2_len = T2::needs_bytes();
        let t1 = T1::try_from_random_bytes_inner(&bytes[..t1_len])?;
        let t2 = T2::try_from_random_bytes_inner(&bytes[t1_len..t1_len + t2_len])?;
        let t3 = T3::try_from_random_bytes_inner(&bytes[t1_len + t2_len..])?;
        Some((t1, t2, t3))
    }
}

impl<T1, T2, T3, T4> FromRandomBytes<(T1, T2, T3, T4)> for (T1, T2, T3, T4)
where
    T1: FromRandomBytes<T1>,
    T2: FromRandomBytes<T2>,
    T3: FromRandomBytes<T3>,
    T4: FromRandomBytes<T4>,
{
    fn needs_bytes() -> usize {
        T1::needs_bytes() + T2::needs_bytes() + T3::needs_bytes() + T4::needs_bytes()
    }

    fn try_from_random_bytes_inner(bytes: &[u8]) -> Option<(T1, T2, T3, T4)> {
        let t1_len = T1::needs_bytes();
        let t2_len = T2::needs_bytes();
        let t3_len = T3::needs_bytes();
        let t1 = T1::try_from_random_bytes_inner(&bytes[..t1_len])?;
        let t2 = T2::try_from_random_bytes_inner(&bytes[t1_len..t1_len + t2_len])?;
        let t3 =
            T3::try_from_random_bytes_inner(&bytes[t1_len + t2_len..t1_len + t2_len + t3_len])?;
        let t4 = T4::try_from_random_bytes_inner(&bytes[t1_len + t2_len + t3_len..])?;
        Some((t1, t2, t3, t4))
    }
}

impl<T1, T2, T3, T4, T5> FromRandomBytes<(T1, T2, T3, T4, T5)> for (T1, T2, T3, T4, T5)
where
    T1: FromRandomBytes<T1>,
    T2: FromRandomBytes<T2>,
    T3: FromRandomBytes<T3>,
    T4: FromRandomBytes<T4>,
    T5: FromRandomBytes<T5>,
{
    fn needs_bytes() -> usize {
        T1::needs_bytes()
            + T2::needs_bytes()
            + T3::needs_bytes()
            + T4::needs_bytes()
            + T5::needs_bytes()
    }

    fn try_from_random_bytes_inner(bytes: &[u8]) -> Option<(T1, T2, T3, T4, T5)> {
        let t1_len = T1::needs_bytes();
        let t2_len = T2::needs_bytes();
        let t3_len = T3::needs_bytes();
        let t4_len = T4::needs_bytes();
        let t1 = T1::try_from_random_bytes_inner(&bytes[..t1_len])?;
        let t2 = T2::try_from_random_bytes_inner(&bytes[t1_len..t1_len + t2_len])?;
        let t3 =
            T3::try_from_random_bytes_inner(&bytes[t1_len + t2_len..t1_len + t2_len + t3_len])?;
        let t4 = T4::try_from_random_bytes_inner(
            &bytes[t1_len + t2_len + t3_len..t1_len + t2_len + t3_len + t4_len],
        )?;
        let t5 = T5::try_from_random_bytes_inner(&bytes[t1_len + t2_len + t3_len + t4_len..])?;
        Some((t1, t2, t3, t4, t5))
    }
}

impl<T1, T2, T3, T4, T5, T6> FromRandomBytes<(T1, T2, T3, T4, T5, T6)> for (T1, T2, T3, T4, T5, T6)
where
    T1: FromRandomBytes<T1>,
    T2: FromRandomBytes<T2>,
    T3: FromRandomBytes<T3>,
    T4: FromRandomBytes<T4>,
    T5: FromRandomBytes<T5>,
    T6: FromRandomBytes<T6>,
{
    fn needs_bytes() -> usize {
        T1::needs_bytes()
            + T2::needs_bytes()
            + T3::needs_bytes()
            + T4::needs_bytes()
            + T5::needs_bytes()
            + T6::needs_bytes()
    }

    fn try_from_random_bytes_inner(bytes: &[u8]) -> Option<(T1, T2, T3, T4, T5, T6)> {
        let t1_len = T1::needs_bytes();
        let t2_len = T2::needs_bytes();
        let t3_len = T3::needs_bytes();
        let t4_len = T4::needs_bytes();
        let t5_len = T5::needs_bytes();
        let t1 = T1::try_from_random_bytes_inner(&bytes[..t1_len])?;
        let t2 = T2::try_from_random_bytes_inner(&bytes[t1_len..t1_len + t2_len])?;
        let t3 =
            T3::try_from_random_bytes_inner(&bytes[t1_len + t2_len..t1_len + t2_len + t3_len])?;
        let t4 = T4::try_from_random_bytes_inner(
            &bytes[t1_len + t2_len + t3_len..t1_len + t2_len + t3_len + t4_len],
        )?;
        let t5 = T5::try_from_random_bytes_inner(
            &bytes[t1_len + t2_len + t3_len + t4_len..t1_len + t2_len + t3_len + t4_len + t5_len],
        )?;
        let t6 =
            T6::try_from_random_bytes_inner(&bytes[t1_len + t2_len + t3_len + t4_len + t5_len..])?;
        Some((t1, t2, t3, t4, t5, t6))
    }
}

pub trait WithConjugationAutomorphism {
    fn sigma(&self) -> Self;
    fn sigma_vec(vec: &Vector<Self>) -> Vector<Self>
    where
        Self: PolyRing,
    {
        Vector::<Self>::from(vec.iter().map(|x| x.sigma()).collect::<Vec<Self>>())
    }
}
