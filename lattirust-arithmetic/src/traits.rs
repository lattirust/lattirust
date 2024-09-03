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

pub trait WithConjugationAutomorphism {
    fn sigma(&self) -> Self;
    fn sigma_vec(vec: &Vector<Self>) -> Vector<Self>
    where
        Self: PolyRing,
    {
        Vector::<Self>::from(vec.iter().map(|x| x.sigma()).collect::<Vec<Self>>())
    }
}
