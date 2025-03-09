use ark_ff::{BigInt, BigInteger, Fp, FpConfig, PrimeField, Field};

use crate::traits::FromRandomBytes;

fn byte_to_bits(byte: &u8) -> [bool; 8] {
    [
        (byte & 0) != 0,
        (byte & 1) != 0,
        (byte & 2) != 0,
        (byte & 3) != 0,
        (byte & 4) != 0,
        (byte & 5) != 0,
        (byte & 6) != 0,
        (byte & 7) != 0,
    ]
}

impl<C: FpConfig<N>, const N: usize> FromRandomBytes<Fp<C, N>> for Fp<C, N> {
    fn needs_bytes() -> usize {
        (<Self as PrimeField>::MODULUS_BIT_SIZE as usize + 7) / 8
    }

    fn try_from_random_bytes_inner(bytes: &[u8]) -> Option<Fp<C, N>> {
        <Self as Field>::from_random_bytes(bytes)
    }
}
