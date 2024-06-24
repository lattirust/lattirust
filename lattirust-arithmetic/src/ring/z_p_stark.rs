use ark_ff::{ BigInt, BigInteger, Field, Fp, Fp64, FpConfig, MontBackend, MontConfig, PrimeField };
use num_bigint::{ BigUint, ToBigUint };

use crate::ring::{ ConvertibleRing, Ring, SignedRepresentative, UnsignedRepresentative };
use crate::traits::Modulus;

// pSTARK = 2^251 + 17 · 2^192 + 1
pub struct FPStarkConfig {}

impl MontConfig<4> for FPStarkConfig {
    const MODULUS: BigInt<4> = BigInt::<4> {
        0: [0x800000000000011, 0x0000000000000000, 0x0000000000000000, 0x0000000000000001],
    };
    const GENERATOR: Fp<MontBackend<Self, 4>, 4> = Fp::new(BigInt::<4> { 0: [0, 0, 0, 2] });
    // TODO: check if this needed/makes sense
    const TWO_ADIC_ROOT_OF_UNITY: Fp<MontBackend<Self, 4>, 4> = Fp::new(BigInt::<4> {
        0: [0, 0, 0, 0],
    }); // TODO: implement this? Not using todo!() here to generate the docs
}

pub type ZPStark = Fp<MontBackend<FPStarkConfig, 4>, 4>;

// Implementing Modulus trait for ZPStark
impl Modulus for ZPStark {
    fn modulus() -> BigUint {
        FPStarkConfig::MODULUS.into()
    }
}
// Function to create a constant element in ZPStark
pub const fn const_zp_stark_from(val1: u64, val2: u64, val3: u64, val4: u64) -> ZPStark {
    ZPStark::new(BigInt::<4> { 0: [val1, val2, val3, val4] })
}
