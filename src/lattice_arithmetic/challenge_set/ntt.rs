use rand::thread_rng;

use crate::lattice_arithmetic::ring::Ring;
use crate::lattice_arithmetic::traits::Modulus;

pub trait NTT<const N: usize, R: Ring + Modulus> {
    fn primitive_root_of_unity() -> R {
        assert!(N.is_power_of_two());
        assert_eq!(R::modulus() % (2 * N as u64), 1);
        let rng = &mut thread_rng();
        let q = R::modulus();
        let exp = (q - 1) / N as u64;
        let n_half = N as u64 / 2;
        loop {
            let x = R::rand(rng);
            let g = x.pow(&[exp]);
            let g_2 = g.pow(&[n_half]);
            if g_2 != R::one() {
                return g;
            }
        }
    }

    fn root_of_unity_pows_bit_reversed(psi: &R) -> [R; N] {
        assert!(N.is_power_of_two());
        assert_eq!(R::modulus() % (2 * N as u64), 1);
        let logN = N.ilog2();
        let mut pows = [R::zero(); N];
        for i in 0..N {
            let iinv = (i as u64).reverse_bits() >> (64 - logN);
            pows[iinv as usize] = psi.pow(&[i as u64]);
        }
        pows
    }
    fn root_of_unity_neg_pows_bit_reversed(psi: &R) -> [R; N] {
        assert!(N.is_power_of_two());
        assert_eq!(R::modulus() % (2 * N as u64), 1);
        let logN = N.ilog2();
        let mut pows = [R::zero(); N];
        for i in 0..N {
            let iinv = (i as u64).reverse_bits() >> (64 - logN);
            pows[iinv as usize] = psi.pow(&[i as u64]).inverse().unwrap();
        }
        pows
    }

    fn ntt(a: &mut [R], root_of_unity_pows_bit_reversed: &[R; N]) {
        assert!(N.is_power_of_two());
        assert_eq!(R::modulus() % (2 * N as u64), 1);
        let mut t = N;
        let mut m = 1;
        let mut j1: usize;
        let mut j2: usize;
        let mut s: R;
        while m < N {
            t = t / 2;
            for i in 0..m {
                j1 = 2 * i * t;
                j2 = j1 + t - 1;
                s = root_of_unity_pows_bit_reversed[m + i];
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

    fn intt(a: &mut [R], root_of_unity_neg_pows_bit_reversed: &[R; N]) {
        assert!(N.is_power_of_two());
        assert_eq!(R::modulus() % (2 * N as u64), 1);
        let mut t = 1;
        let mut m = N;
        let mut j1: usize;
        let mut j2: usize;
        let mut h: usize;
        let mut s: R;
        while m > 1 {
            j1 = 0;
            h = m / 2;
            for i in 0..h {
                j2 = j1 + t - 1;
                s = root_of_unity_neg_pows_bit_reversed[h + i];
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
        let n_inv = R::from(N as u128).inverse().unwrap();
        for i in 0..N {
            a[i] *= n_inv;
        }
    }
}

#[cfg(test)]
mod tests {
    use num_traits::One;
    use crate::lattice_arithmetic::ring::Zq;
    use super::*;

    const Q: u64 = 2u64.pow(16) + 1;

    // 4th Fermat prime, NTT_friendly for N up to 2**15
    type R = Zq<Q>;

    const N: usize = 64; // as used by labrador

    struct NTT_struct {}

    impl NTT<N, R> for NTT_struct {}


    #[test]
    fn test_primitive_root_of_unity() {
        let psi = NTT_struct::primitive_root_of_unity();
        assert_eq!(psi.pow([N as u64]), R::one());
        for i in 1..N {
            assert_ne!(psi.pow([i as u64]), R::one(), "psi^i = {}^{} = {} should not be 1", psi, i, psi.pow([i as u64]));
        }
    }

    #[test]
    fn test_ntt_intt() {
        let psi = NTT_struct::primitive_root_of_unity();
        let psi_pows = NTT_struct::root_of_unity_pows_bit_reversed(&psi);
        let psi_neg_pows = NTT_struct::root_of_unity_neg_pows_bit_reversed(&psi);
        let mut a = (0..N).into_iter().map(|i| R::from(i as u128)).collect::<Vec<R>>();
        let a_original = a.clone();
        NTT_struct::ntt(a.as_mut_slice(), &psi_pows);
        NTT_struct::intt(a.as_mut_slice(), &psi_neg_pows);
        assert_eq!(a_original, a);
    }
}