use num_traits::{One, Zero};
use crate::lattice_arithmetic::poly_ring::PolyRing;

use crate::lattice_arithmetic::pow2_cyclotomic_poly_ring_ntt::Pow2CyclotomicPolyRingNTT;
use crate::lattice_arithmetic::ring::{Ring, Zq};
use crate::lattice_arithmetic::traits::FromRandomBytes;

pub trait ChallengeSet<R: Ring>: FromRandomBytes<R> {}

