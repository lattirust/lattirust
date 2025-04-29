use ark_ff::Field;
use icicle_core::traits::FieldImpl;
use num_bigint::BigUint;
use num_traits::One;
use std::ops::Rem;

use crate::gpu_context::{convert_to_Zq, convert_to_babybear, copy_from_host, copy_to_host, get_default_ntt_config, init_ntt_context_on_device, try_load_and_set_GPU_backend_device};
use crate::ring::{const_fq_from, Zq};
use crate::traits::Modulus;

use log::{info, warn, debug, error};
use pretty_env_logger::env_logger;

use icicle_core::ntt::{self};
use icicle_babybear::field::ScalarField as BabybearField;

use clap::Parser;

use crate::gpu_context::Error::*;


// -----------------------------------------------------------------------------------------------


#[derive(Parser, Debug)]
struct Args {
    size: u8,
    device_type: String,
}

// -----------------------------------------------------------------------------------------------

//noinspection RsAssertEqual
/// Return q such that 2^(bit_size-1) <= q < 2^bit_size and q mod 2*N = 1
pub const fn ntt_modulus<const N: usize>(bit_size: usize) -> u64 {
    assert!(bit_size < 64, "bit_size must be less than 64");
    let n_bit_size = 64 - (N as u64).leading_zeros();
    let res = 1 + (N as u64) * (1 << (bit_size - n_bit_size as usize));
    assert!(64 - res.leading_zeros() == bit_size as u32);
    res
}

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
    let mut t: i64 = 0;
    let mut new_t: i64 = 1;
    let mut r: i64 = Q as i64;
    let mut new_r: i64 = x as i64;

    while new_r != 0 {
        let quotient = r / new_r;
        (t, new_t) = (new_t, t - quotient * new_t);
        (r, new_r) = (new_r, r - quotient * new_r);
    }
    if r > 1 {
        panic!("could not find inverse");
    }
    if t < 0 {
        t = t + (Q as i64);
    }
    t as u64
}

//noinspection RsAssertEqual
const fn primitive_root_of_unity<const Q: u64, const N: usize>() -> u64 {
    assert!(N.is_power_of_two());
    assert!(
        Q % (2 * N as u64) == 1,
        "Q is not NTT-friendly, i.e., not equal to 1 mod 2N"
    );
    let exp = (Q - 1) / N as u64;
    let n_half = N as u64 / 2;
    // Ideally, this would be a randomized algorithm, randomness in const functions is hard, so we iterate through all integers mod Q instead.
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
const fn root_of_unity_pows_bit_reversed<const Q: u64, const N: usize>(psi: u64) -> [Zq<Q>; N] {
    assert!(N.is_power_of_two());
    assert!(
        Q % (2 * N as u64) == 1,
        "Q is not NTT-friendly, i.e., not equal to 1 mod 2N"
    );
    let log_n = N.ilog2();
    let mut pows = [Zq::<Q>::ZERO; N];
    let mut i = 0;
    while i < N {
        let iinv = (i as u64).reverse_bits() >> (64 - log_n);
        pows[iinv as usize] = const_fq_from(const_pow_mod::<Q>(psi, i as u64));
        i += 1;
    }
    pows
}

//noinspection RsAssertEqual
const fn root_of_unity_neg_pows_bit_reversed<const Q: u64, const N: usize>(psi: u64) -> [Zq<Q>; N] {
    assert!(N.is_power_of_two());
    assert!(
        Q % (2 * N as u64) == 1,
        "Q is not NTT-friendly, i.e., not equal to 1 mod 2N"
    );
    let log_n = N.ilog2();
    let mut pows = [Zq::<Q>::ZERO; N];
    let mut i = 0;
    let psi_inv = const_inv_mod::<Q>(psi);
    while i < N {
        let iinv = (i as u64).reverse_bits() >> (64 - log_n);
        pows[iinv as usize] = const_fq_from(const_pow_mod::<Q>(psi_inv, i as u64));
        i += 1;
    }
    pows
}


fn gpu_ntt_acceleration<const Q: u64, const N: usize>(
    size: usize,
    input: &mut [Zq<Q>; N],
    direction: ntt::NTTDir,
) -> Result<(), crate::gpu_context::Error> {

    #[cfg(feature = "GPU")]
    {

        if try_load_and_set_GPU_backend_device().val() != ErrNone.val() {
            return Err(GpuNotAvailable);
        }

        let (mut d_input, mut d_output) = init_ntt_context_on_device(size).map_err(|_| IcicleBackendNotFound)?;
        let data = convert_to_babybear(*input);

        copy_from_host(&mut d_input, data).map_err(|_| CopyFailed)?;
        let config = get_default_ntt_config();

        ntt::ntt(&d_input, direction, &config, &mut d_output[..]).map_err(|_| InvalidInput)?;
        let mut host_babybear_results = vec![BabybearField::zero(); size];

        copy_to_host(&mut host_babybear_results, d_output).map_err(|_| CopyFailed)?;
        let res: Result<[Zq<Q>; N], _> = convert_to_Zq(host_babybear_results);

        res.map_err(|_| InvalidInput)?.iter().enumerate().for_each(|(i, val)| input[i] = *val);

        return Ok(());
    }
    #[cfg(not(feature = "GPU"))]
    {
        warn!("GPU feature not enabled at compile time, using CPU implementation.");
    }

    // Default to CPU
    return Err(GpuNotAvailable)

}
// -----------------------------------------------------------------------------------------------

pub trait NTT<const Q: u64, const N: usize> {
    const ROOT_OF_UNITY: u64 = primitive_root_of_unity::<Q, N>();
    const POWS_ROOT_OF_UNITY: [Zq<Q>; N] = root_of_unity_pows_bit_reversed(Self::ROOT_OF_UNITY);
    const NEG_POWS_ROOT_OF_UNITY: [Zq<Q>; N] =
        root_of_unity_neg_pows_bit_reversed(Self::ROOT_OF_UNITY);

    fn ntt(a: &mut [Zq<Q>; N]) {
        assert!(N.is_power_of_two());
        assert!(Zq::<Q>::modulus().rem(BigUint::from(2 * N)).is_one());
        let mut t = N;
        let mut m = 1;
        let mut j1: usize;
        let mut j2: usize;
        let mut s: Zq<Q>;

        let direction = ntt::NTTDir::kForward;
        let res = gpu_ntt_acceleration::<Q, N>(t, a, direction);
        if res.is_ok() {
            return;
        }

        while m < N {
            t = t / 2;
            for i in 0..m {
                j1 = 2 * i * t;
                j2 = j1 + t - 1;
                s = Self::POWS_ROOT_OF_UNITY[m + i];
                for j in j1..j2 + 1 {
                    let u = a[j];
                    let v = a[j + t] * s;
                    a[j] = u + v;
                    a[j + t] = u - v;
                }
            }
            m *= 2;
        }
    }

    fn intt(a: &mut [Zq<Q>; N]) {
        assert!(N.is_power_of_two());
        assert!(Zq::<Q>::modulus().rem(BigUint::from(2 * N)).is_one());
        let mut t = 1;
        let mut m = N;
        let mut j1: usize;
        let mut j2: usize;
        let mut h: usize;
        let mut s: Zq<Q>;

        let direction = ntt::NTTDir::kInverse;
        let res = gpu_ntt_acceleration::<Q, N>(m, a, direction);
        if res.is_ok() {
            return;
        }

        while m > 1 {
            j1 = 0;
            h = m / 2;
            for i in 0..h {
                j2 = j1 + t - 1;
                s = Self::NEG_POWS_ROOT_OF_UNITY[h + i];
                for j in j1..j2 + 1 {
                    let u = a[j];
                    let v = a[j + t];
                    a[j] = u + v;
                    a[j + t] = (u - v) * s;
                }
                j1 += 2 * t;
            }
            t *= 2;
            m /= 2;
        }
        let n_inv = Zq::<Q>::from(N as u128).inverse().unwrap();
        for i in 0..N {
            a[i] *= n_inv;
        }
    }

    fn ntt_coeffs(&self) -> Vec<Zq<Q>>;
}

#[cfg(test)]
mod tests {
    use num_traits::One;

    use crate::ring::Zq;

    use super::*;

    // const Q: u64 = 2u64.pow(16) + 1;
    const Q: u64 = 0x78000001;

    const N: usize = 64; // as used by labrador

    struct NttStruct {}

    impl NTT<Q, N> for NttStruct {
        fn ntt_coeffs(&self) -> Vec<Zq<Q>> {
            unimplemented!()
        }
    }

    #[test]
    fn test_const_inv_mod() {
        for x in 1..Q {
            let x_inv = const_inv_mod::<Q>(x);
            assert_eq!(x * x_inv % Q, 1);
        }
    }

    #[test]
    fn test_primitive_root_of_unity() { 
        let psi = Zq::<Q>::from(primitive_root_of_unity::<Q, N>());
        assert_eq!(psi.pow([N as u64]), Zq::<Q>::one());
        for i in 1..N {
            assert_ne!(
                psi.pow([i as u64]),
                Zq::<Q>::one(),
                "psi^i = {}^{} = {} should not be 1",
                psi,
                i,
                psi.pow([i as u64])
            );
        }
    }
    
    #[test]
    fn test_ntt_intt() {
        let mut a: [Zq<Q>; N] = core::array::from_fn(|i| Zq::<Q>::from(i as u128));

        let a_original = a.clone();
        NttStruct::ntt(&mut a);
        NttStruct::intt(&mut a);
        assert_eq!(a_original, a);
    }
}
