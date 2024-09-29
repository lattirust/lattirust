use ark_std::ops::{BitXor, DivAssign};
use num_bigint::{BigInt, BigUint};
use num_integer::Integer;
use num_traits::sign::Signed;

use crate::{
    balanced_decomposition::{decompose_balanced, Decompose},
    traits::{WithL2Norm, WithLinfNorm},
    Ring,
};

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
        + Into<BigInt>
        + DivAssign<Self::SignedInt>
        + From<i128>
        + Clone
        + BitXor<Output = Self::SignedInt>
        + Sized;
}

impl<R: ConvertibleRing> Decompose for R {
    fn decompose(&self, b: u128, padding_size: Option<usize>) -> Vec<Self> {
        decompose_balanced(self, b, padding_size)
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
