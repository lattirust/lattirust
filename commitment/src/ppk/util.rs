use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};

use crate::bfv::public_key::PublicKey;
// use lattirust_arithmetic::nimue::serialization::{FromBytes, ToBytes}

#[derive(Clone, Copy, Debug, CanonicalSerialize, CanonicalDeserialize, )]
pub struct PublicParameters<const Q: u64, const P: u64, const N: usize> {
    pub pk: PublicKey<Q, P, N>,
    pub q: u64, // ciphertext modulo
    pub p: u64, // plaintext modulo
    pub n: usize, // degree of the polynomial
}

impl<const Q: u64, const P: u64, const N: usize> Default for PublicParameters<Q, P, N> {
    fn default() -> Self {
        Self {
            pk: PublicKey::default(),
            q: Q,
            p: P,
            n: N,
        }
    }
}


