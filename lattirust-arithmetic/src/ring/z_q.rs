use ark_ff::{BigInt, Field, Fp, Fp64, FpConfig, MontBackend, MontConfig, PrimeField};
use num_bigint::{BigUint, ToBigUint};

use crate::ring::{ConvertibleRing, Ring, SignedRepresentative, UnsignedRepresentative};
use crate::traits::Modulus;

pub struct FqConfig<const Q: u64> {}

impl<const Q: u64> MontConfig<1> for FqConfig<Q> {
    const MODULUS: BigInt<1> = BigInt::<1>([Q]);
    const GENERATOR: Fp<MontBackend<Self, 1>, 1> = Fp::new(BigInt::<1>([2u64]));
    // TODO: check if this needed/makes sense
    const TWO_ADIC_ROOT_OF_UNITY: Fp<MontBackend<Self, 1>, 1> = Fp::new(BigInt::<1>([0u64])); // TODO: implement this? Not using todo!() here to generate the docs
}

pub type Zq<const Q: u64> = Fp64<MontBackend<FqConfig<Q>, 1>>;

impl<const Q: u64> const Modulus for Zq<Q> {
    fn modulus() -> BigUint {
        Q.to_biguint().unwrap()
    }
}

pub const fn const_fq_from<const Q: u64>(val: u64) -> Zq<Q> {
    Zq::new(BigInt::<1>([val]))
}

impl<const Q: u64> Ring for Zq<Q> {
    const ZERO: Self = <Zq<Q> as Field>::ZERO;
    const ONE: Self = <Zq<Q> as Field>::ONE;

    fn inverse(&self) -> Option<Self> {
        ark_ff::Field::inverse(&self)
    }
}

impl From<BigInt<1>> for UnsignedRepresentative {
    fn from(value: BigInt<1>) -> Self {
        UnsignedRepresentative(value.0[0] as u128)
    }
}

impl<C: FpConfig<1>> From<Fp64<C>> for UnsignedRepresentative {
    fn from(value: Fp64<C>) -> Self {
        UnsignedRepresentative::from(value.into_bigint())
    }
}

/// Map [0, q) to (-h, h] using [0, h] -> [0, h] and (h, q) -> (-h, 0), where h = floor((q-1)/2)
impl<C: FpConfig<1>> From<Fp64<C>> for SignedRepresentative {
    fn from(value: Fp64<C>) -> Self {
        let unsigned = UnsignedRepresentative::from(value).0 as i128;
        let v: BigInt<1> = value.into();
        let q = UnsignedRepresentative::from(Fp64::<C>::MODULUS).0 as i128;
        if v > Fp64::<C>::MODULUS_MINUS_ONE_DIV_TWO {
            SignedRepresentative(unsigned - q)
        } else {
            SignedRepresentative(unsigned)
        }
    }
}

/// Map (-h, h] to [0, q) using [0, h] -> [0, h] and (-h, 0) -> (h, q), where h = floor((q-1)/2)
impl<C: FpConfig<1>> From<SignedRepresentative> for Fp64<C> {
    fn from(value: SignedRepresentative) -> Self {
        let q: i128 = UnsignedRepresentative::from(Fp64::<C>::MODULUS).0 as i128;
        if value.0 < 0 {
            Fp64::<C>::from((value.0 + q) as u128)
        } else {
            Fp64::<C>::from(value.0 as u128)
        }
    }
}

impl<const Q: u64> ConvertibleRing for Zq<Q> {}
