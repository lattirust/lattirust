use crate::lattice_arithmetic::ring::Ring;
use crate::lattice_arithmetic::traits::FromRandomBytes;

pub trait ChallengeSet<R>: FromRandomBytes<R> {}

