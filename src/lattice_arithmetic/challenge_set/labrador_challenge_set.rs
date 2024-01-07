use bitter::BitReader;
use nalgebra::ComplexField;
use num_traits::{One, Zero};

use crate::lattice_arithmetic::challenge_set::challenge::ChallengeSet;
use crate::lattice_arithmetic::matrix::Matrix;
use crate::lattice_arithmetic::poly_ring::PolyRing;
use crate::lattice_arithmetic::pow2_cyclotomic_poly_ring::Pow2CyclotomicPolyRing;
use crate::lattice_arithmetic::pow2_cyclotomic_poly_ring_ntt::Pow2CyclotomicPolyRingNTT;
use crate::lattice_arithmetic::ring::{Ring, Zq};
use crate::lattice_arithmetic::traits::FromRandomBytes;

pub struct LabradorChallengeSet<R: PolyRing> {
    _marker: std::marker::PhantomData<R>,
}


/// Challenge set for Zq[X]/(X^64+1) where entries have 23 zero coefficients, 31 coefficients with value ±1, and 10 coefficients with value ±2
/// There are more than 2^128 such elements, and they all have l2-norm 71.
/// In addition, rejection sampling is used to restrict to challenges with operator norm at most 15.
/// On average, 6 elements need to be sampled to get some c with ||c||_op < 15.
/// Differences of distinct challenges are invertible.
impl<R: PolyRing> LabradorChallengeSet<R> {
    type Ring = R;
    type BaseRing = R::BaseRing;

    pub fn is_valid(r: &R) -> bool {
        // let coeffs = r.coeffs();
        // assert_eq!(coeffs.len(), 64);
        // let num_zero = coeffs.into_iter().filter(|c| c.is_zero()).count();
        // let num_one = coeffs.into_iter().filter(|c| c.is_one() || (-*c).is_one()).count();
        // let num_two = coeffs.into_iter().filter(|c| c == &R::BaseRing::from(2) || (-*c) == R::BaseRing::from(2)).count();
        //
        // num_zero == 23 && num_one == 31 && num_two == 10 && Self::operator_norm(&coeffs) < 15.
        todo!()
    }
}

impl<const Q: u64, const N: usize> FromRandomBytes<Pow2CyclotomicPolyRingNTT<Zq<Q>, N>> for LabradorChallengeSet<Pow2CyclotomicPolyRingNTT<Zq<Q>, N>> {
    fn byte_size() -> usize {
        6 * N * Self::BaseRing::byte_size()
    }

    fn from_random_bytes(bytes: &[u8]) -> Option<Self::Ring> {
        assert!(bytes.len() >= Self::byte_size());
        todo!()
    }
}

impl<const Q: u64, const N: usize> ChallengeSet<Pow2CyclotomicPolyRingNTT<Zq<Q>, N>> for LabradorChallengeSet<Pow2CyclotomicPolyRingNTT<Zq<Q>, N>> {}

impl<const Q: u64, const N: usize> LabradorChallengeSet<Pow2CyclotomicPolyRing<Zq<Q>, N>> {
    const EXPECTED_OPERATOR_NORM_REJECTION_SAMPLES: usize = 6;
    const OPERATOR_NORM_THRESHOLD: f64 = 15.;

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
        let num_zeros = 23;
        let num_pm_ones = 31;
        let num_pm_twos = 10;

        // Sample 31 + 10 sign bits
        let num_sign_bits: usize = num_pm_ones + num_pm_twos;
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
        let sample_bytesize = Self::permutation_byte_size() + Self::sign_bits_bytesize();
        loop {
            let coeffs = Self::unchecked_coeffs_from_random_bytes(&bytes[i * sample_bytesize..(i + 1) * sample_bytesize]);
            if Self::operator_norm(&coeffs) < Self::OPERATOR_NORM_THRESHOLD {
                return coeffs;
            }
            i += 1;
            if i == Self::EXPECTED_OPERATOR_NORM_REJECTION_SAMPLES {
                eprintln!("Could not sample a challenge with operator norm < 15 after {} attempts", Self::EXPECTED_OPERATOR_NORM_REJECTION_SAMPLES);
            }
            if i >= 100 {
                panic!("Could not sample a challenge with operator norm < 15 after 100 attempts");
            }
        }
    }

    fn sign_bits_bytesize() -> usize {
        let num_sign_bits: usize = 31 + 10;
        num_sign_bits.div_ceil(8) as usize
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

impl<const Q: u64, const N: usize> FromRandomBytes<Pow2CyclotomicPolyRing<Zq<Q>, N>> for LabradorChallengeSet<Pow2CyclotomicPolyRing<Zq<Q>, N>> {
    fn byte_size() -> usize {
        Self::EXPECTED_OPERATOR_NORM_REJECTION_SAMPLES * Self::sample_byte_size()
    }

    fn from_random_bytes(bytes: &[u8]) -> Option<Self::Ring> {
        assert_eq!(bytes.len(), Self::byte_size());
        Some(Self::Ring::from(Self::checked_coeffs_from_random_bytes(bytes).iter().cloned().map(|c|
            if c >= 0 { Self::BaseRing::from(c as u32) } else { -Self::BaseRing::from(-c as u32) }
        ).collect::<Vec<Self::BaseRing>>()))
    }
}

impl<const Q: u64, const N: usize> LabradorChallengeSet<Pow2CyclotomicPolyRing<Zq<Q>, N>> {
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
        let mut c_mat = Self::challenge_to_matrix(c);
        let eig = (c_mat.transpose() * c_mat).symmetric_eigenvalues().iter().map(|e| e.abs()).into_iter().max_by(|a, b| a.partial_cmp(b).unwrap()).unwrap();
        eig.sqrt()
    }
}

impl<const Q: u64, const N: usize> ChallengeSet<Pow2CyclotomicPolyRing<Zq<Q>, N>> for LabradorChallengeSet<Pow2CyclotomicPolyRing<Zq<Q>, N>> {}

mod tests {
    use ark_std::UniformRand;
    use rand::{Rng, thread_rng};
    use crate::lattice_arithmetic::matrix::{Matrix, Vector};

    use crate::lattice_arithmetic::poly_ring::PolyRing;
    use crate::lattice_arithmetic::pow2_cyclotomic_poly_ring::Pow2CyclotomicPolyRing;
    use crate::lattice_arithmetic::ring::Zq;
    use crate::lattice_arithmetic::traits::{FromRandomBytes, Normed};
    use crate::nimue::serialization::ToBytes;

    use super::LabradorChallengeSet;

    const Q: u64 = 4294967291;
    const D: usize = 64;

    type R = Zq<Q>;

    type PolyR = Pow2CyclotomicPolyRing<R, D>;
    type LabCS = LabradorChallengeSet<PolyR>;

    const NUM_REPETITIONS: usize = 1000;
    const TOLERANCE: f64 = 1e-60;

    #[test]
    fn test_challenge_to_matrix() {
        let mut bytes = vec![0u8; LabCS::byte_size()];
        thread_rng().fill(bytes.as_mut_slice());

        let c = LabCS::checked_coeffs_from_random_bytes(&bytes.as_slice());
        let op_norm = LabCS::operator_norm(&c);
        assert!(op_norm < LabCS::OPERATOR_NORM_THRESHOLD, "||c||_op = {} is not < {}", op_norm, LabCS::OPERATOR_NORM_THRESHOLD);

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
        let mut norm: f64;
        let z = vec![0i8; D];
        norm = LabCS::operator_norm(&z);
        assert_eq!(norm, 0.);

        let mut bytes = vec![0u8; LabCS::byte_size()];
        thread_rng().fill(bytes.as_mut_slice());

        let c = LabCS::checked_coeffs_from_random_bytes(&bytes.as_slice());

        // Check that the operator norm is less than the threshold
        let op_norm = LabCS::operator_norm(&c);
        assert!(op_norm < LabCS::OPERATOR_NORM_THRESHOLD, "||c||_op = {} is not < {}", op_norm, LabCS::OPERATOR_NORM_THRESHOLD);

        // Check that from_random_bytes() is consistent with check_coeffs_from_random_bytes()
        let c_poly = LabradorChallengeSet::<PolyR>::from_random_bytes(&bytes.as_slice()).unwrap();
        assert_eq!(c.len(), c_poly.coeffs().len());
        for i in 0..c.len() {
            let val = if c[i] >= 0 { R::from(c[i] as u32) } else { -R::from(-c[i] as u32) };
            assert!(val == c_poly.coeffs()[i]);
        }

        // TODO: this fails sometimes, but it is not clear why
        // Check that ||c||op is an upper bound on ||cr||_2 / ||r||_2 for many random r
        for _ in 0..NUM_REPETITIONS {
            let r = PolyR::rand(&mut thread_rng());
            let l2 = (c_poly * r).norm();
            assert!((c_poly * r).norm_squared() as f64 <= (op_norm * op_norm * r.norm_squared() as f64) + TOLERANCE, "||cr||_2^2 = {}  is not <= ||c||_op^2 * ||r||_2^2 = {} * {} = {}", (c_poly * r).norm_squared(), op_norm * op_norm, r.norm_squared(), op_norm * op_norm * (r.norm_squared() as f64));
            assert!(l2 / r.norm() <= op_norm + TOLERANCE, "||cr||_2 / ||r||_2 = {}/{} = {} is not <= ||c||_op = {}", l2, r.norm(), l2 / r.norm(), op_norm);
        }
    }
}