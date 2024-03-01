use bitter::BitReader;

use crate::lattice_arithmetic::challenge_set::challenge::ChallengeSet;
use crate::lattice_arithmetic::matrix::Matrix;
use crate::lattice_arithmetic::poly_ring::PolyRing;
use crate::lattice_arithmetic::pow2_cyclotomic_poly_ring::Pow2CyclotomicPolyRing;
use crate::lattice_arithmetic::pow2_cyclotomic_poly_ring_ntt::Pow2CyclotomicPolyRingNTT;
use crate::lattice_arithmetic::ring::Fq;
use crate::lattice_arithmetic::traits::FromRandomBytes;

pub struct LabradorChallengeSet<R: PolyRing> {
    _marker: std::marker::PhantomData<R>,
}


/// Challenge set for Fq[X]/(X^64+1) where entries have 23 zero coefficients, 31 coefficients with value ±1, and 10 coefficients with value ±2
/// There are more than 2^128 such elements, and they all have l2-norm 71.
/// In addition, rejection sampling is used to restrict to challenges with operator norm at most 15.
/// On average, 6 elements need to be sampled to get some c with ||c||_op < 15.
/// Differences of distinct challenges are invertible.
impl<R: PolyRing> LabradorChallengeSet<R> {
    pub type Field = R;
    pub type BaseRing = R::BaseRing;

    #[allow(dead_code)]
    const NUM_ZEROS: usize = 23;
    const NUM_PM_ONES: usize = 31;
    const NUM_PM_TWOS: usize = 10;
    #[allow(dead_code)]
    const NUM_COEFFS: usize = Self::NUM_ZEROS + Self::NUM_PM_ONES + Self::NUM_PM_TWOS; // 64

    pub const OPERATOR_NORM_THRESHOLD: f64 = 15.;
    pub const L2_NORM_SQUARED: f64 = (1 * Self::NUM_PM_ONES + 2 * 2 * Self::NUM_PM_TWOS) as f64; // 71

    // Return the variance of a the sum of the coefficients of a challenge polynomial
    pub const VARIANCE_SUM_COEFFS: f64 = (Self::NUM_PM_ONES + 2 * Self::NUM_PM_TWOS) as f64; // 51
}

impl<const Q: u64, const N: usize> LabradorChallengeSet<Pow2CyclotomicPolyRing<Fq<Q>, N>> {
    const CUTOFF_OPERATOR_NORM_REJECTION_SAMPLES: usize = 64;
    // TODO: find a value with a solid theoretical justification

    /// Returns an integer uniformly sampled from [0, p[ using rejection sampling, as well as the unused random bytes
    fn u8_from_random_bytes(p: u8, bytes: &[u8]) -> Option<(u8, &[u8])> {
        let mut i = 0;
        loop {
            let val = bytes[i] & (p.next_multiple_of(2) - 1);
            if val < p {
                return Some((val, &bytes[i + 1..]));
            }
            i += 1;
            if i >= bytes.len() {
                return None;
            }
        }
    }

    fn unchecked_coeffs_from_random_bytes(bytes: &[u8]) -> Vec<i8> {
        assert!(bytes.len() >= Self::sample_byte_size());
        assert_eq!(N, 64); // The current implementation can easily be generalized to powers of 2, but this is not implemented yet

        // Sample 31 + 10 sign bits
        let num_sign_bits: usize = Self::NUM_PM_ONES + Self::NUM_PM_TWOS;
        let sign_bits_bytesize: usize = num_sign_bits.div_ceil(8);
        let mut sign_bits = bitter::LittleEndianReader::new(&bytes[0..sign_bits_bytesize]);
        let mut bytes = &bytes[sign_bits_bytesize..];

        // Sample a permutation of [0, N=64[ to decide the values, using the Fisher-Yates shuffle
        let mut indices = (0..N).collect::<Vec<usize>>();
        let mut j: u8;
        for i in (1..N - 1).rev() {
            assert!(i > 0);
            (j, bytes) = Self::u8_from_random_bytes((i + 1) as u8, &bytes).unwrap();
            indices.swap(i, j as usize);
        }
        assert_eq!(indices.clone().len(), N);
        for i in 0..N {
            assert!(indices.clone().into_iter().find(|x| *x == i).is_some());
        }

        // Set coefficients according to permutation and sign bits
        let mut coeffs = Vec::with_capacity(N);
        for i in 0..N {
            let val = match indices[i] {
                0..23 => 0,  // 23 times 0
                23..54 => 1, // 31 times ±1
                54..64 => 2, // 10 times ±2
                _ => unreachable!(),
            };
            if val == 0 {
                coeffs.push(0);
            } else {
                let sign_bit = sign_bits.read_bit().unwrap();
                coeffs.push(if sign_bit { -val } else { val });
            }
        }
        assert_eq!(bytes[0..sign_bits_bytesize].len(), sign_bits_bytesize);
        let b_ = sign_bits.bits_remaining().unwrap();
        debug_assert_eq!(b_, sign_bits_bytesize * 8 - 8 * (num_sign_bits / 8) - 1);
        coeffs
    }

    fn checked_coeffs_from_random_bytes(bytes: &[u8]) -> Vec<i8> {
        let mut i = 0;
        let sample_bytesize = Self::sample_byte_size();
        loop {
            let coeffs = Self::unchecked_coeffs_from_random_bytes(&bytes[i * sample_bytesize..(i + 1) * sample_bytesize]);
            if Self::operator_norm(&coeffs) < Self::OPERATOR_NORM_THRESHOLD {
                return coeffs;
            }
            i += 1;
            if i == Self::CUTOFF_OPERATOR_NORM_REJECTION_SAMPLES {
                panic!("Could not sample a challenge with operator norm < {} after {} attempts", Self::OPERATOR_NORM_THRESHOLD, Self::CUTOFF_OPERATOR_NORM_REJECTION_SAMPLES);
            }
        }
    }

    fn sign_bits_bytesize() -> usize {
        let num_sign_bits: usize = 31 + 10;
        num_sign_bits.div_ceil(8)
    }

    fn sample_u8_bytesize(p: u8) -> usize {
        let log2_threshold = -40;
        // Pr[failure] = (1-p/P)^k, where P is the next power of 2 above p, and k is the number of repetitions
        // Pr[failure] <= t <=> k >= log(t) / log(1-p/P)
        let p_pow2 = p.next_power_of_two();
        let k = (log2_threshold as f64) / (1. - p as f64 / p_pow2 as f64).log2();
        k.ceil() as usize // We only need one byte per sample since p < 2^8
    }

    fn permutation_byte_size() -> usize {
        let mut bs = 0usize;
        for i in (1..N - 1).rev() {
            bs += Self::sample_u8_bytesize((i + 1) as u8);
        }
        bs
    }

    fn sample_byte_size() -> usize {
        Self::permutation_byte_size() + Self::sign_bits_bytesize()
    }
}

impl<const Q: u64, const N: usize> FromRandomBytes<Pow2CyclotomicPolyRing<Fq<Q>, N>> for LabradorChallengeSet<Pow2CyclotomicPolyRing<Fq<Q>, N>> {
    fn byte_size() -> usize {
        Self::CUTOFF_OPERATOR_NORM_REJECTION_SAMPLES * Self::sample_byte_size()
    }

    fn try_from_random_bytes(bytes: &[u8]) -> Option<Pow2CyclotomicPolyRing<Fq<Q>, N>> {
        assert_eq!(bytes.len(), Self::byte_size());
        Some(Self::Field::from(Self::checked_coeffs_from_random_bytes(bytes).iter().cloned().map(|c|
            if c >= 0 { Self::BaseRing::from(c as u32) } else { -Self::BaseRing::from(-c as u32) }
        ).collect::<Vec<Self::BaseRing>>()))
    }
}

impl<const Q: u64, const N: usize> LabradorChallengeSet<Pow2CyclotomicPolyRing<Fq<Q>, N>> {
    fn challenge_to_matrix(c: &Vec<i8>) -> Matrix<f64> {
        assert_eq!(c.len(), N);

        let mut c_mat = Matrix::<f64>::zeros(N, N);
        for i in 0..N {
            for j in 0..N - i {
                c_mat[(i + j, j)] = c[i] as f64;
            }
            for j in N - i..N {
                c_mat[(i + j - N, j)] = -(c[i] as f64);
            }
        }
        c_mat
    }

    pub fn operator_norm(c: &Vec<i8>) -> f64 {
        let c_mat = Self::challenge_to_matrix(c);
        let eig = (c_mat.transpose() * c_mat).symmetric_eigenvalues().iter().map(|e| e.abs()).into_iter().max_by(|a, b| a.partial_cmp(b).unwrap()).unwrap();
        eig.sqrt()
    }
}

impl<const Q: u64, const N: usize> ChallengeSet<Pow2CyclotomicPolyRing<Fq<Q>, N>> for LabradorChallengeSet<Pow2CyclotomicPolyRing<Fq<Q>, N>> {}

impl<const Q: u64, const N: usize> FromRandomBytes<Pow2CyclotomicPolyRingNTT<Q, N>> for LabradorChallengeSet<Pow2CyclotomicPolyRingNTT<Q, N>> {
    fn byte_size() -> usize {
        LabradorChallengeSet::<Pow2CyclotomicPolyRing<Fq::<Q>, N>>::byte_size()
    }

    fn try_from_random_bytes(bytes: &[u8]) -> Option<Pow2CyclotomicPolyRingNTT<Q, N>> {
        LabradorChallengeSet::<Pow2CyclotomicPolyRing<Fq::<Q>, N>>::try_from_random_bytes(bytes).map(|x| x.into())
    }
}

impl<const Q: u64, const N: usize> ChallengeSet<Pow2CyclotomicPolyRingNTT<Q, N>> for LabradorChallengeSet<Pow2CyclotomicPolyRingNTT<Q, N>> {}

#[cfg(test)]
mod tests {
    use ark_std::UniformRand;
    use rand::{Rng, thread_rng};

    use crate::lattice_arithmetic::matrix::{Matrix, Vector};
    use crate::lattice_arithmetic::ntt::ntt_modulus;
    use crate::lattice_arithmetic::poly_ring::PolyRing;
    use crate::lattice_arithmetic::pow2_cyclotomic_poly_ring::Pow2CyclotomicPolyRing;
    use crate::lattice_arithmetic::ring::Fq;
    use crate::lattice_arithmetic::traits::{FromRandomBytes, WithL2Norm};

    use super::LabradorChallengeSet;

    const Q: u64 = ntt_modulus::<64>(32);
    const D: usize = 64;

    type R = Fq<Q>;

    type PolyR = Pow2CyclotomicPolyRing<R, 64>;
    type LabCS = LabradorChallengeSet<PolyR>;

    const NUM_REPETITIONS: usize = 1000;
    const TOLERANCE: f64 = 1e-60;

    #[test]
    fn test_challenge_to_matrix() {
        let mut bytes = vec![0u8; LabCS::byte_size()];
        thread_rng().fill(bytes.as_mut_slice());

        let c = LabCS::checked_coeffs_from_random_bytes(&bytes.as_slice());
        let op_norm = LabCS::operator_norm(&c);
        assert!(op_norm <= LabCS::OPERATOR_NORM_THRESHOLD, "||c||_op = {} is not < {}", op_norm, LabCS::OPERATOR_NORM_THRESHOLD);

        let c_poly = PolyR::from_fn(|i| if c[i] >= 0 { R::from(c[i] as u32) } else { -R::from(-c[i] as u32) });

        let c_mat_int = LabCS::challenge_to_matrix(&c);
        let c_mat = Matrix::<R>::from_fn(D, D,
                                         |i, j| if c_mat_int[(i, j)] >= 0. { R::from(c_mat_int[(i, j)] as u32) } else { -R::from(-c_mat_int[(i, j)] as u32) },
        );

        for _ in 0..NUM_REPETITIONS {
            let r = PolyR::rand(&mut thread_rng());
            let r_vec = Vector::<R>::from(r.coeffs());
            let cr_vec = &c_mat * r_vec;
            let cr = PolyR::from(cr_vec.data.as_vec().to_vec());
            assert_eq!(c_poly * r, cr);
        }
    }

    #[test]
    fn test_operator_norm() {
        let norm: f64;
        let z = vec![0i8; D];
        norm = LabCS::operator_norm(&z);
        assert_eq!(norm, 0.);

        let mut bytes = vec![0u8; LabCS::byte_size()];
        thread_rng().fill(bytes.as_mut_slice());

        let c = LabCS::checked_coeffs_from_random_bytes(&bytes.as_slice());

        // Check that the operator norm is less than the threshold
        let op_norm = LabCS::operator_norm(&c);
        assert!(op_norm <= LabCS::OPERATOR_NORM_THRESHOLD, "||c||_op = {} is not < {}", op_norm, LabCS::OPERATOR_NORM_THRESHOLD);

        // Check that from_random_bytes() is consistent with check_coeffs_from_random_bytes()
        let c_poly = LabradorChallengeSet::<PolyR>::try_from_random_bytes(&bytes.as_slice()).unwrap();
        assert_eq!(c.len(), c_poly.coeffs().len());
        for i in 0..c.len() {
            let val = if c[i] >= 0 { R::from(c[i] as u32) } else { -R::from(-c[i] as u32) };
            assert!(val == c_poly.coeffs()[i]);
        }

        // TODO: this fails sometimes, but it is not clear why
        // Check that ||c||op is an upper bound on ||cr||_2 / ||r||_2 for many random r
        for _ in 0..NUM_REPETITIONS {
            let r = PolyR::rand(&mut thread_rng());
            let l2 = (c_poly * r).l2_norm();
            assert!((c_poly * r).l2_norm_squared() as f64 <= (op_norm * op_norm * r.l2_norm_squared() as f64) + TOLERANCE, "||cr||_2^2 = {}  is not <= ||c||_op^2 * ||r||_2^2 = {} * {} = {}", (c_poly * r).l2_norm_squared(), op_norm * op_norm, r.l2_norm_squared(), op_norm * op_norm * (r.l2_norm_squared() as f64));
            assert!(l2 / r.l2_norm() <= op_norm + TOLERANCE, "||cr||_2 / ||r||_2 = {}/{} = {} is not <= ||c||_op = {}", l2, r.l2_norm(), l2 / r.l2_norm(), op_norm);
        }
    }
}