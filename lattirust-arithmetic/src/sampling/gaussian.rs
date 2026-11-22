//! Integer Gaussian masks with nominal standard deviation 1.55 * 2^k.
//!
//! Rust port of Biasioli et al., lazer-crypto/labrador c1ae0e6, gaussian.c,
//! Copyright (C) 2026, IBM Corp., Apache-2.0. The centered sampler uses the
//! same 128-bit CDF and floating point Bernoulli procedure as the reference.
//! The four test vectors were generated independently by that C routine.
//! This reference sampler is variable time. See NOTICE in this directory.
use ark_std::rand::{CryptoRng, RngCore};

const BLOCK: usize = 512;
const CDF: &[(u64, u64)] = &[
    (10894764499197476522, 10804844707381617341),
    (4761708367981796450, 6209732027000382074),
    (1476784279527800432, 14108379346150813303),
    (316388870594767345, 17827298407763885637),
    (46043503515468600, 18385899657892021654),
    (4503729779335039, 3860889375818664979),
    (294122444862326, 13947176349836216550),
    (12769598070895, 14321894682751135119),
    (367552986472, 10286761328368440884),
    (7001273393, 2287787188970898528),
    (88153536, 17843977990435837663),
    (733119, 12174894787802692461),
    (4024, 18067426722645776197),
    (14, 10764017821655913055),
    (0, 643125733022530080),
    (0, 1014291014134832),
    (0, 1055215183460),
    (0, 724109373),
    (0, 327744),
    (0, 98),
    (0, 0),
];

struct Buffered<'a, G> {
    rng: &'a mut G,
    bytes: [u8; BLOCK],
    used: usize,
    bits: u64,
    available: u32,
}

impl<G: RngCore> Buffered<'_, G> {
    fn take<const N: usize>(&mut self) -> [u8; N] {
        if self.used + N > BLOCK {
            self.rng.fill_bytes(&mut self.bytes);
            self.used = 0;
        }
        let out = self.bytes[self.used..self.used + N].try_into().unwrap();
        self.used += N;
        out
    }
    fn pair(&mut self) -> (u64, u64) {
        let b = self.take::<16>();
        (
            u64::from_le_bytes(b[..8].try_into().unwrap()),
            u64::from_le_bytes(b[8..].try_into().unwrap()),
        )
    }
    fn bit(&mut self) -> i32 {
        if self.available == 0 {
            self.bits = u64::from_le_bytes(self.take());
            self.available = 64;
        }
        let bit = (self.bits & 1) as i32;
        self.bits >>= 1;
        self.available -= 1;
        bit
    }
    /// Bernoulli(exp(-x)); the small exponential matches the fdlibm
    /// approximation used in the upstream sampler.
    fn bernoulli(&mut self, x: f64) -> bool {
        let (first, second) = self.pair();
        let exponent = (x * std::f64::consts::LOG2_E) as u32;
        let residue = -(x - std::f64::consts::LN_2 * f64::from(exponent));
        let t = residue * residue;
        let correction = residue
            - t * (1.6666666666666602e-1
                + t * (-2.7777777777015593e-3
                    + t * (6.613756321437934e-5
                        + t * (-1.6533902205465252e-6 + t * 4.1381367970572385e-8))));
        let probability = 1.0 - ((residue * correction) / (correction - 2.0) - residue);
        let mask = (1u64 << exponent.min(63)) - 1;
        (first & mask) == 0
            && (second & ((1u64 << 53) - 1)) < (probability * ((1u64 << 53) as f64)) as u64
    }
    fn centered(&mut self, center: f64) -> i64 {
        loop {
            let sign = self.bit();
            let random = self.pair();
            let absolute = CDF
                .iter()
                .take_while(|&&threshold| random < threshold)
                .count() as i32;
            let k = if sign == 1 { absolute + 1 } else { -absolute };
            let delta = f64::from(k) - center;
            let base = f64::from(k - sign);
            if self.bernoulli((delta * delta - base * base) / (2.0 * 1.55 * 1.55)) {
                return i64::from(k);
            }
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SamplingError {
    InvalidScale,
    SizeOverflow,
}

/// Samples the integer masks used by the LNP reduction. The wider result
/// avoids overflowing a signed i32 before the protocol's public tail check.
pub fn sample_masks<G: RngCore + CryptoRng>(
    out: &mut [i64],
    k: u32,
    rng: &mut G,
) -> Result<(), SamplingError> {
    if k > 32 {
        return Err(SamplingError::InvalidScale);
    }
    if out.is_empty() {
        return Ok(());
    }
    let bits = out
        .len()
        .checked_mul(k as usize)
        .ok_or(SamplingError::SizeOverflow)?;
    let size = bits
        .div_ceil(BLOCK * 8)
        .checked_add(1)
        .and_then(|n| n.checked_mul(BLOCK))
        .ok_or(SamplingError::SizeOverflow)?;
    let mut residues = vec![0; size];
    rng.fill_bytes(&mut residues);
    let mut rand = Buffered {
        rng,
        bytes: [0; BLOCK],
        used: BLOCK,
        bits: 0,
        available: 0,
    };
    let radix = 1i64 << k;
    for (i, sample) in out.iter_mut().enumerate() {
        let bit = i * k as usize;
        let window = u64::from_le_bytes(residues[bit / 8..bit / 8 + 8].try_into().unwrap());
        let residue = ((window >> (bit % 8)) & ((radix - 1) as u64)) as i64;
        let center = residue as f64 / radix as f64;
        *sample = rand.centered(center) * radix - residue;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_std::rand::Error;
    // Fixed reference entropy for cross-language tests only.
    struct Fixture(u64);
    impl RngCore for Fixture {
        fn next_u64(&mut self) -> u64 {
            self.0 = self.0.wrapping_add(0x9e3779b97f4a7c15);
            let z = (self.0 ^ (self.0 >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
            let z = (z ^ (z >> 27)).wrapping_mul(0x94d049bb133111eb);
            z ^ (z >> 31)
        }
        fn next_u32(&mut self) -> u32 {
            self.next_u64() as u32
        }
        fn fill_bytes(&mut self, d: &mut [u8]) {
            for chunk in d.chunks_mut(8) {
                chunk.copy_from_slice(&self.next_u64().to_le_bytes()[..chunk.len()]);
            }
        }
        fn try_fill_bytes(&mut self, d: &mut [u8]) -> Result<(), Error> {
            self.fill_bytes(d);
            Ok(())
        }
    }
    impl CryptoRng for Fixture {}
    #[test]
    fn matches_ccs_reference_mask_vectors() {
        for (k, expected) in [
            (0, [-1, 2, 3, 1, -1, 0, 1, -1, 1, 2, 1, -1, 1, 0, 0, -2]),
            (
                3,
                [3, 1, -6, 24, 2, 20, -25, -11, -3, 3, 2, 20, 9, 4, -11, -10],
            ),
            (
                9,
                [
                    67, 208, -217, 1227, 281, 1382, -1999, -918, -242, 175, 172, 1471, 594, 183,
                    -859, -605,
                ],
            ),
            (
                15,
                [
                    7747, 10548, -14822, 90470, 21321, 76706, -99156, -40305, -22249, 25507, 28508,
                    73242, 49334, 5906, -60866, -39234,
                ],
            ),
        ] {
            let mut values = [0; 16];
            sample_masks(&mut values, k, &mut Fixture(177)).unwrap();
            assert_eq!(values, expected);
        }
    }
}
