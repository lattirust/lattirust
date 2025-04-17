use ark_ff::{BigInt, Fp, Fp64, MontBackend};

use crate::ring::{fq_zero, Fq, FqConfig, Ring};

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

const fn const_pows<const Q: u64, const N: usize>(x: u64) -> [Fq<Q>; N] {
    let mut pows = [0; N];
    let mut pows_fq = [fq_zero(); N];
    pows[0] = 1;
    pows_fq[0] = const_fq_from(1);
    let mut i = 1;
    while i < N {
        let pow: u128 = (pows[i - 1] as u128).checked_mul(x as u128).unwrap();
        pows[i] = (pow % (Q as u128)) as u64;
        pows_fq[i] = const_fq_from(pows[i]);
        i += 1;
    }
    pows_fq
}

/// Return inverse of x modulo Q
const fn const_inv_mod<const Q: u64>(x: u64) -> u64 {
    let mut t: i128 = 0;
    let mut new_t: i128 = 1;
    let mut r: i128 = Q as i128;
    let mut new_r: i128 = x as i128;

    while new_r != 0 {
        let quotient = r / new_r;
        (t, new_t) = (
            new_t,
            t.checked_sub(quotient.checked_mul(new_t).unwrap()).unwrap(),
        );
        (r, new_r) = (
            new_r,
            r.checked_sub(quotient.checked_mul(new_r).unwrap()).unwrap(),
        );
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
// TODO: look up in a precomputed table for N powers of 2 and reasonable bit sizes
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
#[allow(dead_code)]
const fn nth_primitive_root_of_unity<const Q: u64, const N: usize>() -> u64 {
    assert!(N.is_power_of_two());
    assert!(
        Q % (N as u64) == 1,
        "Q is not NTT-friendly, i.e., not equal to 1 mod N"
    );
    let exp = (Q - 1) / N as u64;
    let n_half = (N / 2) as u64;
    // Typically, this would be a randomized algorithm. Randomness in const functions is hard, so we iterate through all integers mod Q instead.
    let mut i = 1;
    while i < Q {
        // g = x^((q-1)/N) is an n-th root of unity, since (x^((q-1)/N))^N = x^(q-1) = 1
        let g = const_pow_mod::<Q>(i, exp);
        // If g_2 = g^(N/2) = x^((q-1)/2) != 1, then g is a primitive n-th root of unity
        let g_2 = const_pow_mod::<Q>(g, n_half);
        if g_2 != 1 {
            return g;
        }
        i += 1;
    }
    panic!("No primitive n-th root of unity found");
}

// noinspection RsAssertEqual
const fn two_nth_primitive_root_of_unity<const Q: u64, const N: usize>() -> u64 {
    assert!(N.is_power_of_two());
    assert!(
        Q % (2 * N as u64) == 1,
        "Q is not NTT-friendly, i.e., not equal to 1 mod 2N"
    );
    let exp = (Q - 1) / (2 * N) as u64;
    let two_n_half = N as u64;
    // Typically, this would be a randomized algorithm. Randomness in const functions is hard, so we iterate through all integers mod Q instead.
    let mut i = 1;
    while i < Q {
        // g = x^((q-1)/(2N)) is an 2n-th root of unity, since (x^((q-1)/(2N)))^2N = x^(q-1) = 1
        let g = const_pow_mod::<Q>(i, exp);
        // If g_2 = g^(2N/2) = x^((q-1)/2) != 1, then g is a primitive 2n-th root of unity
        let g_2 = const_pow_mod::<Q>(g, two_n_half);
        if g_2 != 1 {
            return g;
        }
        i += 1;
    }
    panic!("No primitive 2n-th root of unity found");
}

const fn bit_reversed_index<const N: usize>(i: usize) -> usize {
    let log_n = N.ilog2();
    let iinv = (i as u64).reverse_bits() >> (64 - log_n);
    iinv as usize
}

//noinspection RsAssertEqual
const fn pows_bit_reversed<const Q: u64, const N: usize>(psi: u64) -> [Fq<Q>; N] {
    assert!(N.is_power_of_two());
    assert!(
        Q % (2 * N as u64) == 1,
        "Q is not NTT-friendly, i.e., not equal to 1 mod 2N"
    );
    let pows = const_pows::<Q, N>(psi);
    let mut pows_bit_reversed = [fq_zero(); N];
    let mut i = 0;
    while i < N {
        pows_bit_reversed[i] = pows[bit_reversed_index::<N>(i)];
        i += 1;
    }
    pows_bit_reversed
}

const fn const_fq_from<const Q: u64>(x: u64) -> Fp64<MontBackend<FqConfig<Q>, 1>> {
    if x >= Q {
        panic!("x must be less than Q");
    }
    // Perform Montgomery reduction
    if x == 0 {
        Fp::new_unchecked(BigInt::<1>([x]))
    } else {
        // Put in Montgomery form
        // Hacky, but the methods we need are either non-const or private :(
        let r = Fq::<Q>::R.0[0];
        let x_ = (((x as u128) * (r as u128)) % (Q as u128)) as u64;
        Fp::new_unchecked(BigInt::<1>([x_]))
    }
}

struct RootOfUnity<const Q: u64, const N: usize> {}

impl<const Q: u64, const N: usize> RootOfUnity<Q, N> {
    const ROOT_OF_UNITY: u64 = two_nth_primitive_root_of_unity::<Q, N>();
    const POWS_ROOT_OF_UNITY_BIT_REVERSED: [Fq<Q>; N] = pows_bit_reversed(Self::ROOT_OF_UNITY);
    const NEG_POWS_ROOT_OF_UNITY_BIT_REVERSED: [Fq<Q>; N] =
        pows_bit_reversed(const_inv_mod::<Q>(Self::ROOT_OF_UNITY));
    const N_INV_MOD_Q: Fq<Q> = const_fq_from(const_inv_mod::<Q>(N as u64));
}

impl<const Q: u64, const N: usize> Ntt<N> for Fq<Q> {
    /// Computes the NTT of the given coefficients in place.
    /// Following Algorithm 1 of https://eprint.iacr.org/2016/504.pdf
    fn ntt(coeffs: &mut [Self; N])
    where
        Self: Sized,
    {
        let mut t = N;
        let mut m = 1;
        let mut j1: usize;
        let mut j2: usize;
        let mut s: Fq<Q>;
        let psi_bitrev = RootOfUnity::<Q, N>::POWS_ROOT_OF_UNITY_BIT_REVERSED;
        while m < N {
            t /= 2;
            for i in 0..m {
                j1 = 2 * i * t;
                j2 = j1 + t - 1;
                s = psi_bitrev[m + i];
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
        let mut t = 1;
        let mut m = N;
        let mut j1: usize;
        let mut j2: usize;
        let mut h: usize;
        let mut s: Fq<Q>;
        let psi_inv_bitrev = RootOfUnity::<Q, N>::NEG_POWS_ROOT_OF_UNITY_BIT_REVERSED;
        while m > 1 {
            j1 = 0;
            h = m / 2;
            for i in 0..h {
                j2 = j1 + t - 1;
                s = psi_inv_bitrev[h + i];
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
        for i in 0..N {
            evals[i] *= RootOfUnity::<Q, N>::N_INV_MOD_Q;
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

    // 2^16+1 is the 4th Fermat prime, NTT-friendly for N up to 2**15
    const N65537: usize = 2usize.pow(15);
    const Q65537: u64 = 2u64.pow(16) + 1;

    // 274177 and 67280421310721 are the prime factors of the LaBRADOR modulus 2^64 + 1
    // Both are NTT-friendly for N up to 128
    const N274177: usize = 128;
    const Q274177: u64 = 274177;
    const N67280421310721: usize = 128;
    const Q67280421310721: u64 = 67280421310721;

    const NMAX_16BITS: usize = 2usize.pow(10);
    const Q16BITS: u64 = ntt_prime::<NMAX_16BITS>(16);
    const NMAX_32BITS: usize = 2usize.pow(11);
    #[allow(long_running_const_eval)]
    const Q32BITS: u64 = ntt_prime::<NMAX_32BITS>(32);
    const NMAX_62BITS: usize = 2usize.pow(11);
    #[allow(long_running_const_eval)]
    const Q62BITS: u64 = ntt_prime::<NMAX_62BITS>(62);

    const NUM_SAMPLES: u64 = 128;

    macro_rules! test_const_fq_from {
        ($Q:expr) => {
            paste::expr! {
                #[test]
                fn [< test_const_fq_from_ $Q >] () {
                    use ark_std::UniformRand;
                    use ark_std::Zero;

                    let zero = const_fq_from::<$Q>(0);
                    assert!(zero.is_zero());
                    assert_eq!(zero, Fq::<$Q>::zero());

                    let one = const_fq_from::<$Q>(1);
                    assert!(one.is_one());
                    assert_eq!(one, Fq::<$Q>::one());

                    let rng = &mut ark_std::test_rng();
                    let mut a = u64::rand(rng);
                    if a >= $Q {
                        a = a % $Q;
                    }
                    assert_eq!(const_fq_from::<$Q>(a), Fq::<$Q>::from(a));
                }
            }
        };
    }

    test_const_fq_from!(Q65537);
    test_const_fq_from!(Q274177);
    test_const_fq_from!(Q67280421310721);
    test_const_fq_from!(Q16BITS);
    test_const_fq_from!(Q32BITS);
    test_const_fq_from!(Q62BITS);

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

    macro_rules! test_two_nth_primitive_root_of_unity {
        ($Q:expr, $N:expr) => {
            paste::expr! {
                #[test]
                fn [< test_two_nth_primitive_root_of_unity_ $Q _ $N >] () {
                    use ark_ff::Field;

                    let psi = Fq::<$Q>::from(two_nth_primitive_root_of_unity::<$Q, $N>());
                    assert_eq!(psi.pow([($N) as u64]), -Fq::<$Q>::one());
                    assert_eq!(psi.pow([(2 * $N) as u64]), Fq::<$Q>::one());
                    for i in 1..(2 * $N) {
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
        };
    }

    test_two_nth_primitive_root_of_unity!(Q65537, N65537);
    test_two_nth_primitive_root_of_unity!(Q274177, N274177);
    test_two_nth_primitive_root_of_unity!(Q67280421310721, N67280421310721);
    test_two_nth_primitive_root_of_unity!(Q16BITS, NMAX_16BITS);
    test_two_nth_primitive_root_of_unity!(Q32BITS, NMAX_32BITS);
    test_two_nth_primitive_root_of_unity!(Q62BITS, NMAX_62BITS);

    macro_rules! test_ntt_prime {
        ($B: expr, $($N:expr),*) => {
            paste::expr! {
                $(
                #[test]
                fn [< test_ntt_prime_ $B bits_N $N >] () {
                    let q = ntt_prime::<$N>($B);
                    assert!(q >= 1 << ($B - 1));
                    assert!(q < 1 << ($B as usize));
                    assert!(const_primes::is_prime(q));
                    assert_eq!(q % (2 * $N as u64), 1);
                }
                )*
            }
        };
    }

    test_ntt_prime!(14, 64, 128, 256, 512);
    test_ntt_prime!(30, 64, 128, 256, 512);
    test_ntt_prime!(62, 64, 128, 256, 512);

    macro_rules! test_ntt_intt {
        ($Q:expr,$($N:expr),*) => {
            $(
            paste::expr! {
                #[test]
                fn [< test_ntt_intt_ $Q _N $N >] () {
                    use ark_std::UniformRand;
                    let rng = &mut ark_std::test_rng();
                    let mut a: [Fq<$Q>; $N] = core::array::from_fn(|_| Fq::<$Q>::rand(rng));

                    let a_original = a.clone();
                    Fq::<$Q>::ntt(&mut a);
                    Fq::<$Q>::intt(&mut a);
                    assert_eq!(a_original, a);
                }
            }
            )*
        };
    }

    test_ntt_intt!(Q65537, 64, 128, 256, 512, 1024, 2048, 4096, 8192);
    test_ntt_intt!(Q274177, 64, 128);
    test_ntt_intt!(Q67280421310721, 64, 128);
    test_ntt_intt!(Q16BITS, 64, 128, 256, 512, 1024);
    test_ntt_intt!(Q32BITS, 64, 128, 256, 512, 1024, 2048);
    test_ntt_intt!(Q62BITS, 64, 128, 256, 512, 1024, 2048);

    macro_rules! test_ntt_add{
        ($Q:expr,$($N:expr),*) => {
            $(
            paste::expr! {
                #[test]
                fn [< test_ntt_add_ $Q _N $N >] () {
                    use ark_std::UniformRand;

                    let rng = &mut ark_std::test_rng();
                    let mut a: [Fq<$Q>; $N] = core::array::from_fn(|_| Fq::<$Q>::rand(rng));
                    let mut b: [Fq<$Q>; $N] = core::array::from_fn(|_| Fq::<$Q>::rand(rng));
                    let a_plus_b: [Fq<$Q>; $N] = core::array::from_fn(|i| a[i] + b[i]);

                    Fq::<$Q>::ntt(&mut a);
                    Fq::<$Q>::ntt(&mut b);
                    let mut a_plus_b_: [Fq<$Q>; $N] = core::array::from_fn(|i| a[i] + b[i]);
                    Fq::<$Q>::intt(&mut a_plus_b_);

                    assert_eq!(a_plus_b, a_plus_b_);
                }
            }
            )*
        };
    }

    test_ntt_add!(Q65537, 64, 128, 256, 512, 1024, 2048, 4096, 8192);
    test_ntt_add!(Q274177, 64, 128);
    test_ntt_add!(Q67280421310721, 64, 128);
    test_ntt_add!(Q16BITS, 64, 128, 256, 512, 1024);
    test_ntt_add!(Q32BITS, 64, 128, 256, 512, 1024, 2048);
    test_ntt_add!(Q62BITS, 64, 128, 256, 512, 1024, 2048);

    macro_rules! test_ntt_mul{
        ($Q:expr,$($N:expr),*) => {
            $(
            paste::expr! {
                #[test]
                fn [< test_ntt_mul_ $Q _N $N >] () {
                    use ark_std::UniformRand;

                    let rng = &mut ark_std::test_rng();
                    let mut a: [Fq<$Q>; $N] = core::array::from_fn(|_| Fq::<$Q>::rand(rng));
                    let mut b: [Fq<$Q>; $N] = core::array::from_fn(|_| Fq::<$Q>::rand(rng));

                    let mut a_mul_b: [Fq<$Q>; $N] = core::array::from_fn(|_| fq_zero());
                    for i in 0..$N {
                        for j in 0..$N {
                            if i+j < $N {
                                a_mul_b[i+j] += a[i] * b[j];
                            } else {
                                a_mul_b[i+j-$N] -= a[i] * b[j];
                            }
                        }
                    }

                    Fq::<$Q>::ntt(&mut a);
                    Fq::<$Q>::ntt(&mut b);
                    let mut a_mul_b_: [Fq<$Q>; $N] = core::array::from_fn(|i| a[i] * b[i]);
                    Fq::<$Q>::intt(&mut a_mul_b_);

                    assert_eq!(a_mul_b, a_mul_b_);
                }
            }
            )*
        };
    }

    test_ntt_mul!(Q65537, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192);
    test_ntt_mul!(Q274177, 32, 64, 128);
    test_ntt_mul!(Q67280421310721, 32, 64, 128);
    test_ntt_mul!(Q16BITS, 32, 64, 128, 256, 512, 1024);
    test_ntt_mul!(Q32BITS, 32, 64, 128, 256, 512, 1024, 2048);
    test_ntt_mul!(Q62BITS, 32, 64, 128, 256, 512, 1024, 2048);
}
