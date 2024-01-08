use crate::lattice_arithmetic::matrix::Vector;
use crate::lattice_arithmetic::poly_ring::PolyRing;

pub trait Normed<T> {
    /// Returns the L2 norm.
    fn norm(&self) -> f64;
    /// Returns the squared L2 norm.
    fn norm_squared(&self) -> T;
}

pub trait IntegerDiv<Rhs = Self> {
    /// Divides `self` by `rhs`, returning the quotient.
    fn integer_div(&self, rhs: Rhs) -> Self;

    fn div_round(&self, rhs: Rhs) -> Self;
}

pub trait WithLog2 {
    fn log2(&self) -> f64;
    fn log2_q() -> f64;
}

pub trait Modulus {
    fn modulus() -> u64;
}

pub trait FromRandomBytes<T> {
    fn byte_size() -> usize;
    fn from_random_bytes(bytes: &[u8]) -> Option<T>;
}

pub trait WithConjugationAutomorphism {
    fn sigma(&self) -> Self;
    fn sigma_vec(vec: &Vector<Self>) -> Vector<Self> where Self: PolyRing {
        Vector::<Self>::from(vec.iter().map(|x| x.sigma()).collect::<Vec<Self>>())
    }
}