use ark_crypto_primitives::sponge::Absorb;
use ark_ff::Field;
use ark_std::ops::Mul;

use crate::{traits::FromRandomBytes, Ring};

pub trait PolyRing: Ring + From<Vec<Self::BaseRing>> + FromRandomBytes<Self> + From<u128> {
    type BaseRing: Ring;

    fn coeffs(&self) -> &[Self::BaseRing];
    fn into_coeffs(self) -> Vec<Self::BaseRing>;

    fn dimension() -> usize;

    fn from_scalar(scalar: Self::BaseRing) -> Self;
}

pub trait OverField:
    PolyRing<BaseRing: Field<BasePrimeField: Absorb>>
    + Mul<Self::BaseRing, Output = Self>
    + From<Self::BaseRing>
{
}
