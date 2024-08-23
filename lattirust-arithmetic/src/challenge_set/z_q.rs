use ark_ff::{BigInt, BigInteger, Fp, FpConfig, PrimeField};

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
        Self::MODULUS_BIT_SIZE as usize / 8
    }

    fn try_from_random_bytes_inner(bytes: &[u8]) -> Option<Fp<C, N>> {
        let mut r = BigInt::from_bits_be(
            &bytes
                .iter()
                .flat_map(|b| byte_to_bits(b))
                .collect::<Vec<bool>>(),
        );
        // This is safe, since we're using FromRandomBytesForModulus, which makes sure that taking the modulus only introduces negligible bias
        // let mut r = ark_ff::Fp::<C, N>::new_unchecked(bigint);
        while r > Self::MODULUS {
            let borrow = r.sub_with_borrow(&Self::MODULUS);
            assert_eq!(borrow, false);
        }
        Self::from_bigint(r)
    }
}
