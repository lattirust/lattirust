use std::marker::PhantomData;

use ark_ff::BigInt;
use ark_ff::Field;

use crate::ring::{fq_zero, Fq, Ring};

/// Return x^k mod Q
const fn const_pow_mod<const Q: u64>(x: u64, k: u64) -> u64 {
    let mut res: u128 = 1;
    let mut x: u128 = x as u128;
    let mut k: u128 = k as u128;
    while k > 0 {
        if k % 2 == 1 {
            res = (res * x) % (Q as u128);
        }
        x = (x * x) % (Q as u128);
        k /= 2;
    }
    res as u64
}

/// Return inverse of x modulo Q
const fn const_inv_mod<const Q: u64>(x: u64) -> u64 {
    let mut t: i128 = 0;
    let mut new_t: i128 = 1;
    let mut r: i128 = Q as i128;
    let mut new_r: i128 = x as i128;

    while new_r != 0 {
        let quotient = r / new_r;
        (t, new_t) = (new_t, t - quotient * new_t);
        (r, new_r) = (new_r, r - quotient * new_r);
    }
    if r > 1 {
        panic!("could not find inverse");
    }
    if t < 0 {
        t = t.checked_add(Q as i128).unwrap();
    }
    t as u64
}

/// Return a prime q such that 2^(bit_size-1) <= q < 2^bit_size and q mod 2*N = 1
// noinspection RsAssertEqual
pub const fn ntt_prime<const N: usize>(bit_size: usize) -> u64 {
    assert!(bit_size < 64, "bit_size must be less than 64");
    let mut p = const_primes::next_prime(1 << (bit_size - 1)).unwrap();
    while p % (2 * N as u64) != 1 && p.leading_zeros() >= (64 - bit_size as u32) {
        p = const_primes::next_prime(p + 1).unwrap();
    }
    if p.leading_zeros() >= (64 - bit_size as u32) {
        p
    } else {
        panic!("No prime found for the given bit_size and N");
    }
}

// noinspection RsAssertEqual
const fn primitive_root_of_unity<const Q: u64, const N: usize>() -> u64 {
    assert!(N.is_power_of_two());
    assert!(
        Q % (2 * N as u64) == 1,
        "Q is not NTT-friendly, i.e., not equal to 1 mod 2N"
    );
    let exp = (Q - 1) / N as u64;
    let n_half = N as u64 / 2;
    // Ideally, this would be a randomized algorithm. Randomness in const functions is hard, so we iterate through all integers mod Q instead.
    let mut i = 1;
    while i < Q {
        let g = const_pow_mod::<Q>(i, exp);
        let g_2 = const_pow_mod::<Q>(g, n_half);
        if g_2 != 1 {
            return g;
        }
        i += 1;
    }
    panic!("No primitive root of unity found");
}

//noinspection RsAssertEqual
const fn root_of_unity_pows_bit_reversed<const Q: u64, const N: usize>(psi: u64) -> [Fq<Q>; N] {
    assert!(N.is_power_of_two());
    assert!(
        Q % (2 * N as u64) == 1,
        "Q is not NTT-friendly, i.e., not equal to 1 mod 2N"
    );
    let log_n = N.ilog2();
    let mut pows = [fq_zero(); N];
    let mut i = 0;
    while i < N {
        let iinv = (i as u64).reverse_bits() >> (64 - log_n);
        pows[iinv as usize] = const_fq_from(const_pow_mod::<Q>(psi, i as u64));
        i += 1;
    }
    pows
}

const fn const_fq_from<const Q: u64>(x: u64) -> Fq<Q> {
    Fq::<Q> {
        0: BigInt::<1>([x]),
        1: PhantomData,
    }
}

//noinspection RsAssertEqual
const fn root_of_unity_neg_pows_bit_reversed<const Q: u64, const N: usize>(psi: u64) -> [Fq<Q>; N] {
    assert!(N.is_power_of_two());
    assert!(
        Q % (2 * N as u64) == 1,
        "Q is not NTT-friendly, i.e., not equal to 1 mod 2N"
    );
    let log_n = N.ilog2();
    let mut pows = [fq_zero(); N];
    let mut i = 0;
    let psi_inv = const_inv_mod::<Q>(psi);
    while i < N {
        let iinv = (i as u64).reverse_bits() >> (64 - log_n);
        pows[iinv as usize] = const_fq_from(const_pow_mod::<Q>(psi_inv, i as u64));
        i += 1;
    }
    pows
}

struct RootOfUnity<const Q: u64, const N: usize> {}

impl<const Q: u64, const N: usize> RootOfUnity<Q, N> {
    const ROOT_OF_UNITY: u64 = primitive_root_of_unity::<Q, N>();
    const POWS_ROOT_OF_UNITY: [Fq<Q>; N] = root_of_unity_pows_bit_reversed(Self::ROOT_OF_UNITY);
    const NEG_POWS_ROOT_OF_UNITY: [Fq<Q>; N] =
        root_of_unity_neg_pows_bit_reversed(Self::ROOT_OF_UNITY);
}

impl<const Q: u64, const N: usize> Ntt<N> for Fq<Q> {
    fn ntt(coeffs: &mut [Self; N])
    where
        Self: Sized,
    {
        // TODO: figure out a way to assert these conditions at compile-time
        // assert!(N.is_power_of_two());
        // assert!(BigUint::from(Q).rem(BigUint::from(2 * N)).is_one());
        let mut t = N;
        let mut m = 1;
        let mut j1: usize;
        let mut j2: usize;
        let mut s: Fq<Q>;
        while m < N {
            t /= 2;
            for i in 0..m {
                j1 = 2 * i * t;
                j2 = j1 + t - 1;
                s = RootOfUnity::<Q, N>::POWS_ROOT_OF_UNITY[m + i];
                for j in j1..j2 + 1 {
                    let u = coeffs[j];
                    let v = coeffs[j + t] * s;
                    coeffs[j] = u + v;
                    coeffs[j + t] = u - v;
                }
            }
            m *= 2;
        }
    }

    fn intt(evals: &mut [Self; N])
    where
        Self: Sized,
    {
        // TODO: figure out a way to assert these conditions at compile-time
        // assert!(N.is_power_of_two());
        // assert!(BigUint::from(Q).rem(BigUint::from(2 * N)).is_one());
        let mut t = 1;
        let mut m = N;
        let mut j1: usize;
        let mut j2: usize;
        let mut h: usize;
        let mut s: Fq<Q>;
        while m > 1 {
            j1 = 0;
            h = m / 2;
            for i in 0..h {
                j2 = j1 + t - 1;
                s = RootOfUnity::<Q, N>::NEG_POWS_ROOT_OF_UNITY[h + i];
                for j in j1..j2 + 1 {
                    let u = evals[j];
                    let v = evals[j + t];
                    evals[j] = u + v;
                    evals[j + t] = (u - v) * s;
                }
                j1 += 2 * t;
            }
            t *= 2;
            m /= 2;
        }
        let n_inv = Fq::<Q>::from(N as u128).inverse().unwrap();
        for a_i in evals.iter_mut() {
            *a_i *= n_inv;
        }
    }
}

pub trait Ntt<const N: usize> {
    fn ntt(coeffs: &mut [Self; N])
    where
        Self: Sized;

    fn intt(evals: &mut [Self; N])
    where
        Self: Sized;
}

pub trait NttRing<const N: usize>: Ntt<N> + Ring {}
impl<T, const N: usize> NttRing<N> for T where T: Ntt<N> + Ring {}

#[cfg(test)]
mod tests {
    use num_traits::One;

    use super::*;

    const N: usize = 64;

    // 2^16+1 is the 4th Fermat prime, NTT-friendly for N up to 2**15
    const Q65537: u64 = 2u64.pow(16) + 1;

    // 274177 and 67280421310721 are the prime factors of the LaBRADOR modulus 2^64 + 1
    const Q274177: u64 = 274177;
    const Q67280421310721: u64 = 67280421310721;

    const Q16BITS: u64 = ntt_prime::<N>(16);
    const Q32BITS: u64 = ntt_prime::<N>(32);
    const Q62BITS: u64 = ntt_prime::<N>(62);

    const NUM_SAMPLES: u64 = 128;

    macro_rules! test_const_inv_mod {
        ($($Q:expr),*) => {
            $(
            paste::expr! {
                #[test]
                fn [< test_const_inv_mod_ $Q>] () {
                    for x in (1..$Q).step_by((($Q - 1) / NUM_SAMPLES) as usize) {
                        let x_inv = const_inv_mod::<$Q>(x);
                        assert_eq!(Fq::<$Q>::from(x) * Fq::<$Q>::from(x_inv), Fq::<$Q>::one());
                    }
                }
            }
            )*
        };
    }

    test_const_inv_mod!(Q65537, Q274177, Q67280421310721, Q16BITS, Q32BITS, Q62BITS);

    macro_rules! test_primitive_root_of_unity {
        ($($Q:expr),*) => {
            $(
            paste::expr! {
                #[test]
                fn [< test_primitive_root_of_unity_ $Q>] () {
                    let psi = Fq::<$Q>::from(primitive_root_of_unity::<$Q, N>());
                    assert_eq!(psi.pow([N as u64]), Fq::<$Q>::one());
                    for i in 1..N {
                        assert_ne!(
                            psi.pow([i as u64]),
                            Fq::<$Q>::one(),
                            "psi^i = {}^{} = {} should not be 1",
                            psi,
                            i,
                            psi.pow([i as u64])
                        );
                    }
                }
            }
            )*
        };
    }

    test_primitive_root_of_unity!(Q65537, Q274177, Q67280421310721, Q16BITS, Q32BITS, Q62BITS);

    macro_rules! test_ntt_prime {
        ($B: expr, $N: expr) => {
            paste::expr! {
                #[test]
                fn [< test_ntt_prime_ $B _ $N >] () {
                    let q = ntt_prime::<$N>($B);
                    assert!(q >= 1 << ($B - 1));
                    assert!(q < 1 << ($B as usize));
                    assert!(const_primes::is_prime(q));
                    assert_eq!(q % (2 * $N as u64), 1);
                }
            }
        };
    }

    test_ntt_prime!(14, 64);
    test_ntt_prime!(14, 128);
    test_ntt_prime!(14, 256);
    test_ntt_prime!(14, 1024);
    test_ntt_prime!(14, 2048);
    test_ntt_prime!(30, 64);
    test_ntt_prime!(30, 128);
    test_ntt_prime!(30, 256);
    test_ntt_prime!(30, 1024);
    test_ntt_prime!(30, 2048);
    test_ntt_prime!(30, 4096);
    test_ntt_prime!(30, 8192);
    test_ntt_prime!(62, 64);
    test_ntt_prime!(62, 128);
    test_ntt_prime!(62, 256);
    test_ntt_prime!(62, 1024);
    test_ntt_prime!(62, 2048);
    test_ntt_prime!(62, 4096);
    test_ntt_prime!(62, 8192);

    macro_rules! test_ntt_intt {
        ($($Q:expr),*) => {
            $(
            paste::expr! {
                #[test]
                fn [< test_ntt_intt_ $Q>] () {
                    let mut a: [Fq<$Q>; N] = core::array::from_fn(|i| Fq::<$Q>::from(i as u128));
                    let a_original = a.clone();
                    Fq::<$Q>::ntt(&mut a);
                    Fq::<$Q>::intt(&mut a);
                    assert_eq!(a_original, a);
                }
            }
            )*
        };
    }

    test_ntt_intt!(Q65537, Q274177, Q67280421310721, Q16BITS, Q32BITS, Q62BITS);
}
