//! An integer-only sampler for the degree-64 LaBRADOR challenge shape.
//! The certified spectral filter is conservative at the boundary.
//! This is not the legacy sampler's wire format or a security estimate.

use crate::{
    ring::{PolyRing, Pow2CyclotomicPolyRing, Ring},
    traits::FromRandomBytes,
};
use ark_ff::{BigInteger, PrimeField};
use std::marker::PhantomData;

pub struct LabradorExactChallengeSet<R>(PhantomData<R>);

const DEGREE: usize = 64;
const NORM: i64 = 15;
const NORM_SCALE: i64 = 1 << 20;
const L1: i64 = 51; // 31 one-valued + 10 two-valued entries.
const CANDIDATE_BYTES: usize = 6 + (DEGREE - 1) * 8;
const CANDIDATES: usize = 128;

// cos(j*pi/64) * 2^20, absolute error < 1 at each entry. Reproducible via
// pcs/tools/generate_labrador_cosines.py, which certifies the intervals.
const COS: [i64; 128] = [
    1048576, 1047313, 1043527, 1037227, 1028428, 1017151, 1003425, 987281, 968758, 947901, 924761,
    899394, 871859, 842224, 810560, 776944, 741455, 704181, 665210, 624636, 582558, 539076, 494295,
    448324, 401273, 353255, 304386, 254783, 204567, 153858, 102778, 51451, 0, -51451, -102778,
    -153858, -204567, -254783, -304386, -353255, -401273, -448324, -494295, -539076, -582558,
    -624636, -665210, -704181, -741455, -776944, -810560, -842224, -871859, -899394, -924761,
    -947901, -968758, -987281, -1003425, -1017151, -1028428, -1037227, -1043527, -1047313,
    -1048576, -1047313, -1043527, -1037227, -1028428, -1017151, -1003425, -987281, -968758,
    -947901, -924761, -899394, -871859, -842224, -810560, -776944, -741455, -704181, -665210,
    -624636, -582558, -539076, -494295, -448324, -401273, -353255, -304386, -254783, -204567,
    -153858, -102778, -51451, 0, 51451, 102778, 153858, 204567, 254783, 304386, 353255, 401273,
    448324, 494295, 539076, 582558, 624636, 665210, 704181, 741455, 776944, 810560, 842224, 871859,
    899394, 924761, 947901, 968758, 987281, 1003425, 1017151, 1028428, 1037227, 1043527, 1047313,
];

impl<R> LabradorExactChallengeSet<R> {
    pub const OPERATOR_NORM_THRESHOLD: f64 = NORM as f64;
    pub const LINF_NORM: u128 = 2;
    pub const VARIANCE_SUM_COEFFS: f64 = 51.0;
    pub const L2_NORM_SQUARED: f64 = 71.0;
    pub const MAX_CANDIDATES: usize = CANDIDATES;
    pub const BYTES_PER_CANDIDATE: usize = CANDIDATE_BYTES;

    fn candidate(bytes: &[u8]) -> Option<[i8; DEGREE]> {
        let mut permutation: [usize; DEGREE] = std::array::from_fn(|i| i);
        for (i, word) in (1..DEGREE).rev().zip(bytes[6..].chunks_exact(8)) {
            let x = u64::from_le_bytes(word.try_into().ok()?);
            let n = (i + 1) as u64;
            if x < n.wrapping_neg() % n {
                return None;
            }
            permutation.swap(i, (x % n) as usize);
        }
        let mut out = [0; DEGREE];
        let mut next_sign = 0;
        for (o, i) in out.iter_mut().zip(permutation) {
            let v = match i {
                0..23 => continue,
                23..54 => 1,
                _ => 2,
            };
            *o = if bytes[next_sign / 8] >> (next_sign % 8) & 1 == 1 {
                -v
            } else {
                v
            };
            next_sign += 1;
        }
        Some(out)
    }

    // The negacyclic convolution matrix is diagonalized by the unitary
    // Fourier basis at the roots of X^64 + 1. Real c gives conjugate samples.
    // L1 bounds the absolute error of each real and imaginary dot product.
    fn bounded(c: &[i8; DEGREE]) -> bool {
        for j in 0..32 {
            let (mut re, mut im) = (0i64, 0i64);
            for (i, &x) in c.iter().enumerate() {
                let a = ((2 * j + 1) * i) & 127;
                re += i64::from(x) * COS[a];
                im += i64::from(x) * COS[(160 - a) & 127];
            }
            let r = re.abs() + L1;
            let s = im.abs() + L1;
            if r * r + s * s >= (NORM * NORM_SCALE).pow(2) {
                return false;
            }
        }
        true
    }
}

impl<B: Ring, const N: usize> FromRandomBytes<Pow2CyclotomicPolyRing<B, N>>
    for LabradorExactChallengeSet<Pow2CyclotomicPolyRing<B, N>>
{
    fn needs_bytes() -> usize {
        CANDIDATE_BYTES * CANDIDATES
    }

    fn has_no_bias() -> bool {
        // Conditioned on acceptance, each permutation is uniform. Rejection
        // of the whole candidate does not favor any accepted permutation.
        true
    }

    fn try_from_random_bytes_inner(bytes: &[u8]) -> Option<Pow2CyclotomicPolyRing<B, N>> {
        if N != DEGREE {
            return None;
        }
        let c = bytes
            .chunks_exact(CANDIDATE_BYTES)
            .take(CANDIDATES)
            .filter_map(Self::candidate)
            .find(Self::bounded)?;
        let converted = c.try_map(|x| {
            B::try_from(u64::from(x.unsigned_abs()))
                .ok()
                .map(|r| if x < 0 { -r } else { r })
        })?;
        Pow2CyclotomicPolyRing::try_from_coefficients(&converted)
    }
}

/// Uniform prime-field challenge, with a fixed six-candidate transcript.
/// Aborts (rather than introducing a bias) if all masked values exceed p.
pub struct LabradorUniformPrime<P>(PhantomData<P>);
/// Independently sample each coefficient with the same fixed candidate count.
pub struct LabradorUniformPolynomial<R>(PhantomData<R>);

impl<P: PrimeField> FromRandomBytes<P> for LabradorUniformPrime<P> {
    fn needs_bytes() -> usize {
        6 * (P::MODULUS_BIT_SIZE as usize).div_ceil(8)
    }
    fn has_no_bias() -> bool {
        true
    }
    fn try_from_random_bytes_inner(bytes: &[u8]) -> Option<P> {
        let width = (P::MODULUS_BIT_SIZE as usize).div_ceil(8);
        let high_mask = 255u8 >> ((8 - P::MODULUS_BIT_SIZE % 8) % 8);
        bytes.chunks_exact(width).take(6).find_map(|sample| {
            let mut candidate = sample.to_vec();
            candidate[width - 1] &= high_mask;
            let value = P::from_le_bytes_mod_order(&candidate);
            value
                .into_bigint()
                .to_bytes_le()
                .starts_with(&candidate)
                .then_some(value)
        })
    }
}

impl<R: PolyRing> FromRandomBytes<R> for LabradorUniformPolynomial<R>
where
    R::BaseRing: PrimeField,
{
    fn needs_bytes() -> usize {
        R::dimension() * LabradorUniformPrime::<R::BaseRing>::needs_bytes()
    }
    fn has_no_bias() -> bool {
        true
    }
    fn try_from_random_bytes_inner(bytes: &[u8]) -> Option<R> {
        let width = LabradorUniformPrime::<R::BaseRing>::needs_bytes();
        let view = bytes.get(..Self::needs_bytes())?;
        let coefficients = view
            .chunks_exact(width)
            .map(|lane| LabradorUniformPrime::<R::BaseRing>::try_from_random_bytes(lane))
            .collect::<Option<Vec<_>>>()?;
        R::try_from_coefficients(&coefficients)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ring::f_p::Fq;

    type R = Pow2CyclotomicPolyRing<Fq<4294967197>, 64>;
    type Sampler = LabradorExactChallengeSet<R>;

    #[test]
    fn bounded_sampler_has_exact_weights_and_no_panics() {
        use ark_std::rand::{RngCore, SeedableRng};
        let mut rng = ark_std::rand::rngs::StdRng::seed_from_u64(41);
        let mut buf = vec![0; Sampler::byte_size()];
        for _ in 0..8 {
            rng.fill_bytes(&mut buf);
            let c = Sampler::try_from_random_bytes(&buf).unwrap().coefficients();
            assert_eq!(c.iter().filter(|&&x| x == Fq::from(0u64)).count(), 23);
            assert_eq!(
                c.iter()
                    .filter(|&&x| x == Fq::from(1u64) || x == -Fq::from(1u64))
                    .count(),
                31
            );
            assert_eq!(
                c.iter()
                    .filter(|&&x| x == Fq::from(2u64) || x == -Fq::from(2u64))
                    .count(),
                10
            );
        }
        assert!(Sampler::try_from_random_bytes(&[]).is_none());
    }
}
