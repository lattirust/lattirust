use std::ops::Mul;

use ark_ff::{BigInt, BigInteger, Field, Fp, Fp64, FpConfig, PrimeField};
use ark_serialize::CanonicalSerialize;
use bincode::Options;
use num_traits::Zero;
use serde::{Deserialize, Serialize};

use crate::lattice_arithmetic::matrix::Vector;
use crate::lattice_arithmetic::ring::{Fq, Fq2, Ring};
use crate::lattice_arithmetic::traits::{FromRandomBytes, IntegerDiv, Modulus, WithConjugationAutomorphism, WithL2Norm, WithLinfNorm};

pub trait ConvertibleField: Field + Modulus + Into<UnsignedRepresentative> + Into<SignedRepresentative> + From<SignedRepresentative> + FromRandomBytes<Self> {}

pub trait PolyRing:
Ring
+ Mul<Self::BaseRing, Output=Self>
+ From<Vec<Self::BaseRing>>
+ WithConjugationAutomorphism
+ WithL2Norm
+ WithLinfNorm
+ FromRandomBytes<Self>
+ From<u128>
+ From<Self::BaseRing>
+ Serialize
+ for<'a> Deserialize<'a>
{
    type BaseRing: ConvertibleField;

    fn coeffs(&self) -> Vec<Self::BaseRing>;
    fn flattened(vec: &Vector<Self>) -> Vector<Self::BaseRing> {
        Self::flattened_coeffs(vec).into()
    }
    fn flattened_coeffs(vec: &Vector<Self>) -> Vec<Self::BaseRing> {
        vec.as_slice().into_iter().flat_map(|x| x.coeffs()).collect::<Vec<Self::BaseRing>>()
    }
    fn dimension() -> usize;

    fn from_scalar(scalar: Self::BaseRing) -> Self;
}

impl<C: FpConfig<N>, const N: usize> FromRandomBytes<Fp<C, N>> for Fp<C, N> {
    fn byte_size() -> usize {
        Self::zero().uncompressed_size() + 9 // TODO: check if this is correct; this is inferred from Fp<C, N>::from_random_bytes()
    }

    fn try_from_random_bytes(bytes: &[u8]) -> Option<Self> {
        Self::from_random_bytes(&mut bytes.as_ref())
    }
}

impl<const Q: u64> ConvertibleField for Fq<Q> {}

// Norms
impl<F: ConvertibleField> WithL2Norm for F
{
    fn l2_norm_squared(&self) -> u64 {
        Into::<SignedRepresentative>::into(*self).0.pow(2) as u64
    }
}

impl<R: WithL2Norm> WithL2Norm for [R] {
    fn l2_norm_squared(&self) -> u64 {
        self.iter().map(|x| x.l2_norm_squared()).sum()
    }
}

impl<R: WithL2Norm> WithL2Norm for Vec<R> {
    fn l2_norm_squared(&self) -> u64 {
        self.as_slice().l2_norm_squared()
    }
}

impl<R: WithL2Norm> WithL2Norm for Vector<R> {
    fn l2_norm_squared(&self) -> u64 {
        self.as_slice().l2_norm_squared()
    }
}


impl<F: ConvertibleField> WithLinfNorm for F
{
    fn linf_norm(&self) -> u128 {
        Into::<SignedRepresentative>::into(*self).0.abs() as u128
    }
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

impl<R: WithLinfNorm> WithLinfNorm for Vector<R> {
    fn linf_norm(&self) -> u128 {
        self.as_slice().linf_norm()
    }
}


// Work-around to allow us implementing From traits
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct UnsignedRepresentative(pub u128);

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct SignedRepresentative(pub i128);


impl From<BigInt<1>> for UnsignedRepresentative {
    fn from(value: BigInt<1>) -> Self {
        UnsignedRepresentative(value.0[0] as u128)
    }
}

impl<C: ark_ff::FpConfig<1>> From<Fp64<C>> for UnsignedRepresentative {
    fn from(value: Fp64<C>) -> Self {
        UnsignedRepresentative::from(value.into_bigint())
    }
}

/// Map [0, q[ to [-m, m] using [0, m] -> [0, m] and ]m, q[ -> [-m, 0[, where m = (q-1)/2, assuming q is odd
impl<C: ark_ff::FpConfig<1>> From<Fp64<C>> for SignedRepresentative
{
    fn from(value: Fp64<C>) -> Self {
        debug_assert!(Fp64::<C>::MODULUS.is_odd());
        let unsigned = UnsignedRepresentative::from(value).0 as i128;
        let v: BigInt<1> = value.into();
        let q_half: BigInt<1> = Fp64::<C>::MODULUS_MINUS_ONE_DIV_TWO;
        let q = UnsignedRepresentative::from(Fp64::<C>::MODULUS).0 as i128;
        if v > q_half {
            SignedRepresentative(unsigned - q)
        } else {
            SignedRepresentative(unsigned)
        }
    }
}

/// Map [-m, m] to [0, q[ using [0, m] -> [0, m] and [-m, 0[ -> [m, q[, where m = (q-1)/2, assuming q is odd
impl<C: ark_ff::FpConfig<1>> From<SignedRepresentative> for Fp64<C> {
    fn from(value: SignedRepresentative) -> Self {
        debug_assert!(Fp64::<C>::MODULUS.is_odd());
        let q: i128 = UnsignedRepresentative::from(Fp64::<C>::MODULUS).0 as i128;
        if value.0 < 0 {
            Fp64::<C>::from((value.0 + q) as u128)
        } else {
            Fp64::<C>::from(value.0 as u128)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lattice_arithmetic::ring::Fq;

    const Q: u128 = 2u128.pow(61) - 1;
    const Q_HALF: i128 = (Q as i128 - 1) / 2;
    const TEST_STEP_SIZE: usize = (Q / 10) as usize;

    type F = Fq<{ Q as u64 }>;


    #[test]
    fn test_unsigned_representative() {
        for i in (0..Q).step_by(TEST_STEP_SIZE) {
            let f1 = F::from(i);
            let v1 = UnsignedRepresentative::from(f1);
            assert_eq!(i, v1.0);
        }
    }

    #[test]
    fn test_signed_representative() {
        assert_eq!(Q_HALF, F::MODULUS_MINUS_ONE_DIV_TWO.0[0] as i128);
        println!("Q_HALF = {Q_HALF}");
        println!("Q = {Q}");
        for i in (-Q_HALF..=Q_HALF).step_by(TEST_STEP_SIZE) {
            let v1 = SignedRepresentative(i);
            let f2 = F::from(v1.clone());
            let v2 = SignedRepresentative::from(f2);
            assert_eq!(v1.0, v2.0);
        }
    }
}