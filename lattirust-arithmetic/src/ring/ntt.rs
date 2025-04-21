use ark_ff::{BigInt, Fp, Fp64, MontBackend};

use crate::ring::{Fq, fq_zero, FqConfig, Ring};

/// Return x^k mod Q
pub const fn const_pow_mod<const Q: u64>(x: u64, k: u64) -> u64 {
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

pub const fn const_pows<const Q: u64, const N: usize>(x: u64) -> [Fq<Q>; N] {
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
pub const fn const_inv_mod<const Q: u64>(x: u64) -> u64 {
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
pub const fn nth_primitive_root_of_unity<const Q: u64, const N: usize>() -> u64 {
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
pub const fn two_nth_primitive_root_of_unity<const Q: u64, const N: usize>() -> u64 {
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

pub const fn bit_reversed_index<const N: usize>(i: usize) -> usize {
    let log_n = N.ilog2();
    let iinv = (i as u64).reverse_bits() >> (64 - log_n);
    iinv as usize
}

//noinspection RsAssertEqual
pub const fn pows_bit_reversed<const Q: u64, const N: usize>(psi: u64) -> [Fq<Q>; N] {
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

pub const fn const_fq_from<const Q: u64>(x: u64) -> Fp64<MontBackend<FqConfig<Q>, 1>> {
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

pub const fn is_primitive_root_of_unity<const Q: u64>(x: u64, n: u64) -> bool {
    assert!(n.is_power_of_two());
    if (const_pow_mod::<Q>(x, n/2) !=  Q-1) || (const_pow_mod::<Q>(x, n) != 1) {
        return false;
    }
    let mut i = 1;
    while i < n {
        if const_pow_mod::<Q>(x, i) == 1{
            return false;
        }
        i += 1;
    }
    true
}

pub const fn largest_power_of_two_dividing(mut n: u64) -> u32 {
    let mut k = 0;
    while n % 2 == 0 {
        n /= 2;
        k += 1;
    }
    k
}

/// Returns the smallest 2-adic root of unity modulo Q.
/// That is, ω such that ω^(2^k) ≡ 1 mod Q and ω^(2^(k - 1)) ≠ 1 mod Q,
/// where 2^k is the largest power of two dividing Q - 1.
pub const fn two_adic_root_of_unity<const Q: u64>() -> u64 {
    let k = largest_power_of_two_dividing(Q - 1);
    let generator = generator::<Q>();
    const_pow_mod::<Q>(generator, (Q - 1) >> k)
}

pub const fn is_generator<const Q: u64>(g: u64, factors: &[u64; 64]) -> bool {
    let order = Q - 1;
    let mut i = 0;
    while i < 64 {
        let p = factors[i];
        if p == 0 {
            break;
        }
        if const_pow_mod::<Q>(g, order / p) == 1 {
            return false;
        }
        i += 1;
    }
    true
}

/// Returns an array containing the prime factors of `n`.
/// The length of the array is fixed to 64, with remaining slots filled with 0s.
/// This is necessary since we cannot return a dynamically-sized array in a `const fn`.
pub const fn prime_factors(mut n: u64) -> [u64; 64] {
    let mut factors = [0u64; 64];
    let mut index = 0;
    let mut divisor = 2;

    while n > 1 {
        if n % divisor == 0 {
            factors[index] = divisor;
            index += 1;
            n /= divisor;
        } else {
            divisor = const_primes::next_prime(divisor + 1).unwrap();
        }
    }

    factors
}

pub const fn generator<const Q: u64>() -> u64 {
    let factors = prime_factors(Q - 1);
    let mut g = 2;
    while g < Q {
        if is_generator::<Q>(g, &factors) {
            return g;
        }
        g += 1;
    }
    panic!("no generator found");
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
    fn ntt_inplace(coeffs: &mut [Self; N])
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

    fn intt_inplace(evals: &mut [Self; N])
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
        for evals_i in evals.iter_mut() {
            *evals_i *= RootOfUnity::<Q, N>::N_INV_MOD_Q;
        }
    }
}

pub trait Ntt<const N: usize> {
    fn ntt_inplace(coeffs: &mut [Self; N])
    where
        Self: Sized;

    fn intt_inplace(evals: &mut [Self; N])
    where
        Self: Sized;

    #[must_use]
    fn ntt(coeffs: [Self; N]) -> [Self; N]
    where
        Self: Sized + Clone,
    {
        let mut evals = coeffs;
        Self::ntt_inplace(&mut evals);
        evals
    }

    #[must_use]
    fn intt(evals: [Self; N]) -> [Self; N]
    where
        Self: Sized + Clone,
    {
        let mut coeffs = evals;
        Self::intt_inplace(&mut coeffs);
        coeffs
    }
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

    #[macro_export]
    macro_rules! test_ntt_intt {
        ($Tname:ident, $T:ty, $($N:expr),*) => {
            $(
                paste::expr! {
                    #[test]
                    fn [< test_ntt_intt_ $Tname _N $N >] () {
                        use ark_std::UniformRand;
                        let rng = &mut ark_std::test_rng();
                        let mut a: [$T; $N] = core::array::from_fn(|_| $T::rand(rng));

                        let a_original = a.clone();
                        $T::ntt_inplace(&mut a);
                        $T::intt_inplace(&mut a);
                        assert_eq!(a_original, a);
                    }
                }
            )*
        };
    }

    #[macro_export]
    macro_rules! test_ntt_add{
        ($Tname:ident, $T:ty, $($N:expr),*) => {
            $(
                paste::expr! {
                    #[test]
                    fn [< test_ntt_add_ $Tname _N $N >] () {
                        use ark_std::UniformRand;

                        let rng = &mut ark_std::test_rng();
                        let a: [$T; $N] = core::array::from_fn(|_| $T::rand(rng));
                        let b: [$T; $N] = core::array::from_fn(|_| $T::rand(rng));
                        let a_plus_b_naive: [$T; $N] = core::array::from_fn(|i| a[i] + b[i]);

                        let a_ntt = $T::ntt(a);
                        let b_ntt = $T::ntt(b);
                        let a_plus_b_ntt: [$T; $N] = core::array::from_fn(|i| a_ntt[i] + b_ntt[i]);
                        let a_plus_b = $T::intt(a_plus_b_ntt);

                        assert_eq!(a_plus_b, a_plus_b_naive);
                    }
                }
            )*
        };
    }

    #[macro_export]
    macro_rules! test_ntt_mul{
        ($Tname:ident, $T:ty, $($N:expr),*) => {
            $(
                paste::expr! {
                    #[test]
                    fn [< test_ntt_mul_ $Tname _N $N >] () {
                        use ark_std::UniformRand;
                        use num_traits::Zero;

                        let rng = &mut ark_std::test_rng();
                        let a: [$T; $N] = core::array::from_fn(|_| $T::rand(rng));
                        let b: [$T; $N] = core::array::from_fn(|_| $T::rand(rng));

                        let mut a_mul_b_naive: [$T; $N] = core::array::from_fn(|_| $T::zero());
                        for i in 0..$N {
                            for j in 0..$N {
                                if i+j < $N {
                                    a_mul_b_naive[i+j] += a[i] * b[j];
                                } else {
                                    a_mul_b_naive[i+j-$N] -= a[i] * b[j];
                                }
                            }
                        }

                        let a_ntt = $T::ntt(a);
                        let b_ntt = $T::ntt(b);
                        let a_mul_b_ntt: [$T; $N] = core::array::from_fn(|i| a_ntt[i] * b_ntt[i]);
                        let a_mul_b = $T::intt(a_mul_b_ntt);

                        assert_eq!(a_mul_b, a_mul_b_naive);
                    }
                }
            )*
        };
    }

    test_ntt_intt!(Fq65537, Fq::<Q65537>, 64, 128, 256, 512, 1024, 2048, 4096);
    test_ntt_add!(Fq65537, Fq::<Q65537>, 64, 128, 256, 512, 1024, 2048, 4096);
    test_ntt_mul!(Fq65537, Fq::<Q65537>, 64, 128, 256, 512, 1024, 2048, 4096);

    test_ntt_intt!(Fq274177, Fq::<Q274177>, 64, 128);
    test_ntt_add!(Fq274177, Fq::<Q274177>, 64, 128);
    test_ntt_mul!(Fq274177, Fq::<Q274177>, 64, 128);

    test_ntt_intt!(Fq67280421310721, Fq::<Q67280421310721>, 64, 128);
    test_ntt_add!(Fq67280421310721, Fq::<Q67280421310721>, 64, 128);
    test_ntt_mul!(Fq67280421310721, Fq::<Q67280421310721>, 64, 128);

    test_ntt_intt!(Fq16bits, Fq::<Q16BITS>, 64, 128, 256, 512, 1024);
    test_ntt_add!(Fq16bits, Fq::<Q16BITS>, 64, 128, 256, 512, 1024);
    test_ntt_mul!(Fq16bits, Fq::<Q16BITS>, 64, 128, 256, 512, 1024);

    test_ntt_intt!(Fq32bits, Fq::<Q32BITS>, 64, 128, 256, 512, 1024, 2048);
    test_ntt_add!(Fq32bits, Fq::<Q32BITS>, 64, 128, 256, 512, 1024, 2048);
    test_ntt_mul!(Fq32bits, Fq::<Q32BITS>, 64, 128, 256, 512, 1024, 2048);

    test_ntt_intt!(Fq62bits, Fq::<Q62BITS>, 64, 128, 256, 512, 1024, 2048);
    test_ntt_add!(Fq62bits, Fq::<Q62BITS>, 64, 128, 256, 512, 1024, 2048);
    test_ntt_mul!(Fq62bits, Fq::<Q62BITS>, 64, 128, 256, 512, 1024, 2048);
}
