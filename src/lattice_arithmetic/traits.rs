use crate::lattice_arithmetic::matrix::Vector;
use crate::lattice_arithmetic::poly_ring::PolyRing;

pub trait WithL2Norm {
    fn l2_norm(&self) -> f64 {
        (self.l2_norm_squared() as f64).sqrt()
    }
    fn l2_norm_squared(&self) -> u64;
}

pub trait WithLinfNorm {
    fn linf_norm(&self) -> u128;
}

pub trait IntegerDiv<Rhs = Self> {
    /// Divides `self` by `rhs`, returning the quotient.
    fn integer_div(&self, rhs: &Rhs) -> Self;

    fn div_round(&self, rhs: &Rhs) -> Self;
}

pub trait Modulus {
    fn modulus() -> u64;
}

pub trait FromRandomBytes<T> {
    fn byte_size() -> usize;
    fn try_from_random_bytes(bytes: &[u8]) -> Option<T>;
}

pub trait WithConjugationAutomorphism {
    fn sigma(&self) -> Self;
    fn sigma_vec(vec: &Vector<Self>) -> Vector<Self> where Self: PolyRing {
        Vector::<Self>::from(vec.iter().map(|x| x.sigma()).collect::<Vec<Self>>())
    }
}