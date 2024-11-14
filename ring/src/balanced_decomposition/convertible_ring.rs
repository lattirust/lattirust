use ark_std::ops::{Add, BitXor, Div, DivAssign, Rem, Sub};
use num_bigint::{BigInt, BigUint};
use num_integer::Integer;
use num_traits::sign::Signed;

use crate::{
    balanced_decomposition::Decompose,
    traits::{WithL2Norm, WithLinfNorm},
    Ring,
};

use super::decompose_balanced_in_place;

pub trait ConvertibleRing:
    Ring
    + Into<Self::UnsignedInt>
    + Into<Self::SignedInt>
    + From<Self::UnsignedInt>
    + From<Self::SignedInt>
{
    type UnsignedInt: Integer
        + Into<BigUint>
        + DivAssign<Self::UnsignedInt>
        + From<u128>
        + Clone
        + BitXor<Output = Self::UnsignedInt>
        + Sized;
    type SignedInt: Integer
        + Signed
        + Into<BigInt>
        + DivAssign<Self::SignedInt>
        + DivAssign<i128>
        + From<i128>
        + Add<i128, Output = Self::SignedInt>
        + Sub<i128, Output = Self::SignedInt>
        + Div<i128, Output = Self::SignedInt>
        + Rem<i128, Output = Self::SignedInt>
        + Clone
        + BitXor<Output = Self::SignedInt>
        + Sized;
}

impl<R: ConvertibleRing> Decompose for R {
    fn decompose_in_place(&self, b: u128, out: &mut [Self]) {
        decompose_balanced_in_place(self, b, out);
    }
}

// Norms
impl<F: ConvertibleRing> WithL2Norm for F {
    fn l2_norm_squared(&self) -> BigUint {
        let x: <F as ConvertibleRing>::SignedInt = (*self).into();
        let x: BigInt = x.into();

        x.pow(2).try_into().unwrap()
    }
}

impl<F: ConvertibleRing> WithLinfNorm for F {
    fn linf_norm(&self) -> BigUint {
        let x: <F as ConvertibleRing>::SignedInt = (*self).into();
        let x: BigInt = x.into();

        x.abs().to_biguint().unwrap()
    }
}
