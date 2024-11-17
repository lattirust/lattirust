use lattirust_arithmetic::ring::ConvertibleRing;
use crate::bfv::{public_key::PublicKey, Ciphertext};

// TODO: PublicKeu is not a trait so I cannot use it as a generic, do I need to fix it?
pub struct PublicParameters<const Q: u64, const N: usize> {
    pub pk: Ciphertext<Q, N>,
}
