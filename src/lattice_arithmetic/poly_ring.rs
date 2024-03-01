use std::ops::Mul;

use ark_ff::{BigInt, Field, Fp64, PrimeField};

use crate::lattice_arithmetic::matrix::Vector;
use crate::lattice_arithmetic::ring::{Fq, Ring};
use crate::lattice_arithmetic::traits::{FromRandomBytes, IntegerDiv, WithConjugationAutomorphism, WithL2Norm, WithLinfNorm};

pub trait ConvertibleField: Field + Into<UnsignedRepresentative> + Into<SignedRepresentative> + FromRandomBytes<Self> + IntegerDiv {}

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

    fn signed_repr(x: &Self::BaseRing) -> i128;

    fn unsigned_repr(x: &Self::BaseRing) -> u128;
}

impl<const Q: u64> FromRandomBytes<Fq<Q>> for Fq<Q> {
    fn byte_size() -> usize {
        todo!()
    }

    fn try_from_random_bytes(bytes: &[u8]) -> Option<Self> {
        todo!()
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
pub struct UnsignedRepresentative(pub u128);

pub struct SignedRepresentative(pub i128);


impl From<BigInt<1>> for UnsignedRepresentative {
    fn from(value: BigInt<1>) -> Self {
        UnsignedRepresentative(value.0[0] as u128)
    }
}

impl<C: ark_ff::FpConfig<1>> From<Fp64<C>> for UnsignedRepresentative {
    fn from(value: Fp64<C>) -> Self {
        UnsignedRepresentative::from(value.0)
    }
}

impl<C: ark_ff::FpConfig<1>> From<Fp64<C>> for SignedRepresentative
{
    fn from(value: Fp64<C>) -> Self {
        let signed = UnsignedRepresentative::from(value).0 as i128;
        let v: BigInt<1> = value.into();
        let q_half: BigInt<1> = Fp64::<C>::MODULUS_MINUS_ONE_DIV_TWO;
        let q = UnsignedRepresentative::from(Fp64::<C>::MODULUS).0 as i128;
        if v >= q_half {
            SignedRepresentative(signed - q)
        } else {
            SignedRepresentative(signed)
        }
    }
}