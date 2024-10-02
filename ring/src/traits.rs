use lattirust_linear_algebra::{Matrix, Vector};
use num_bigint::BigUint;

use crate::{PolyRing, Ring};

pub trait WithL2Norm {
    fn l2_norm_squared(&self) -> BigUint;
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

pub trait FromRandomBytes<T> {
    fn byte_size() -> usize;
    fn try_from_random_bytes(bytes: &[u8]) -> Option<T>;
}

pub trait Cyclotomic: PolyRing {
    #[inline(always)]
    fn degree() -> usize {
        Self::dimension() - 1
    }

    #[inline]
    fn x() -> Self {
        Self::from(vec![Self::BaseRing::ZERO, Self::BaseRing::ONE])
    }

    fn rot(&mut self);

    fn into_rot_iter(self) -> Rotation<Self> {
        Rotation { curr: self }
    }
}

pub fn rot_matrix<C: Cyclotomic>(mut f: C) -> Matrix<C::BaseRing> {
    let degree = C::degree();
    let mut columns = Vec::with_capacity(degree);

    (0..degree).for_each(|_| {
        f.rot();

        columns.push(Vector::from_slice(f.coeffs()))
    });

    Matrix::from_columns(&columns)
}

pub struct Rotation<R: Cyclotomic> {
    curr: R,
}

impl<R: Cyclotomic> Iterator for Rotation<R> {
    type Item = R;

    fn next(&mut self) -> Option<R> {
        let curr = self.curr;

        self.curr.rot();

        Some(curr)
    }
}
