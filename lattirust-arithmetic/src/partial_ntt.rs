// Partial NTT

use crate::ring::{Ring, Zq};

pub trait PartialNTT<
    const Q: u64,
    const N: usize, // Euler totient function of m for an m-th cyclotomic ring
    const D: usize, // Max degree of the ring splits
    const Z: usize, // Such that ord_m(p) = m / z and p = 1 mod z
    const PHI_Z: usize,
>
{
    fn ntt(a: &mut [Zq<Q>; N], rou: Zq<Q>) {
        // Precompute the root of unity
        // 1. Split array in chunks of size D
        let components = coprimes_set::<Z, PHI_Z>().map(|power| {
            let mut temp = a.clone();
            let rj = rou.pow([power as u64]);
            for i in (D..N).rev() {
                let high_coeff = a[i];

                if high_coeff != Zq::<Q>::ZERO {
                    let target_idx = i - D;
                    temp[target_idx] = temp[target_idx] + high_coeff * rj;
                }
            }
            let mut res = [Zq::<Q>::ZERO; D];
            res.copy_from_slice(&temp[..D]);
            let _ = temp;
            res
        });
        for (j, component) in components.iter().enumerate() {
            for i in 0..D {
                a[j * D + i] = component[i];
            }
        }
    }
}

const fn coprimes_set<const Z: usize, const PHI_Z: usize>() -> [usize; PHI_Z] {
    let mut i = 0;
    let mut j = 0;
    let mut set = [0; PHI_Z];
    while i < Z {
        if gcd::<Z>(i) {
            set[j] = i;
            j += 1;
        }
        i += 1;
    }
    set
}

const fn gcd<const Z: usize>(i: usize) -> bool {
    let mut a = i;
    let mut b = Z;
    while b != 0 {
        let temp = b;
        b = a % b;
        a = temp
    }
    a == 1
}
