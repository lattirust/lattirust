use ark_ff::{ BigInt, Fp, MontBackend, MontConfig };
use num_bigint::BigUint;

use crate::ring::Ring;
use crate::traits::Modulus;

// pgold = 2^64 − 2^32 + 1.
pub struct FPGoldConfig {}

impl MontConfig<1> for FPGoldConfig {
    const MODULUS: BigInt<1> = BigInt::<1> {
        0: [0xffffffff00000001],
    };
    const GENERATOR: Fp<MontBackend<Self, 1>, 1> = Fp::new(BigInt::<1> { 0: [2] });
    // TODO: check if this needed/makes sense
    const TWO_ADIC_ROOT_OF_UNITY: Fp<MontBackend<Self, 1>, 1> = Fp::new(BigInt::<1> {
        0: [0],
    }); // TODO: implement this? Not using todo!() here to generate the docs
}

pub type ZPGold = Fp<MontBackend<FPGoldConfig, 1>, 1>;

// Implementing Modulus trait for ZPStark
impl Modulus for ZPGold {
    fn modulus() -> BigUint {
        FPGoldConfig::MODULUS.into()
    }
}
// Function to create a constant element in ZPStark
pub const fn const_zp_stark_from(val1: u64) -> ZPGold {
    ZPGold::new(BigInt::<1> { 0: [val1] })
}

impl Ring for ZPGold {
    const ZERO: Self = ZPGold::new(BigInt::<1> { 0: [0] });
    const ONE: Self = ZPGold::new(BigInt::<1> { 0: [1] });
}
