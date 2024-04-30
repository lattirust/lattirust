use num_bigint::BigUint;

use crate::linear_algebra::Vector;
use crate::ring::PolyRing;

pub trait WithL2Norm {
    fn l2_norm(&self) -> f64 {
        (self.l2_norm_squared() as f64).sqrt()
    }
    fn l2_norm_squared(&self) -> u128;
}

impl<R: WithL2Norm> WithL2Norm for [R] {
    fn l2_norm_squared(&self) -> u128 {
        self.iter().map(|x| x.l2_norm_squared()).sum()
    }
}

impl<R: WithL2Norm> WithL2Norm for Vec<R> {
    fn l2_norm_squared(&self) -> u128 {
        self.as_slice().l2_norm_squared()
    }
}

pub trait WithLinfNorm {
    fn linf_norm(&self) -> u128;
}

impl<R: WithLinfNorm> WithLinfNorm for [R] {
    fn linf_norm(&self) -> u128 {
        self.iter().map(|x| x.linf_norm()).max().unwrap()
    }
}

impl<R: WithLinfNorm> WithLinfNorm for Vec<R> {
    fn linf_norm(&self) -> u128 {
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

pub trait FromRandomBytes<T> {
    fn byte_size() -> usize;
    fn try_from_random_bytes(bytes: &[u8]) -> Option<T>;
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
