use ark_ff::{Fp, FpConfig, PrimeField};

use crate::traits::FromRandomBytes;

impl<C: FpConfig<N>, const N: usize> FromRandomBytes<Fp<C, N>> for Fp<C, N> {
    fn needs_bytes() -> usize {
        <Self as PrimeField>::MODULUS_BIT_SIZE.div_ceil(8) as usize
    }

    fn try_from_random_bytes_inner(bytes: &[u8]) -> Option<Fp<C, N>> {
        Some(Fp::<C, N>::from_le_bytes_mod_order(bytes))
    }
}
