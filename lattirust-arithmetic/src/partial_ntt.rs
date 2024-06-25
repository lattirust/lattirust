// Partial NTT

use crate::ring::{Ring, Zq};

const fn primitive_root_of_unity<const Q: u64, const Z: usize>() -> u64 {
    // Find the root of unity such that w^z = 1
    assert!((Q - 1) % Z as u64 == 0);
    let resized_power = (Q - 1) / Z as u64;
    todo!()
}

const fn root_of_unity_relative_powers<const Q: u64, const Z: usize, const PHI_Z: usize>(
    rou: u64,
) -> [Zq<Q>; PHI_Z] {
    // 1. Generate set of relative primes to Z
    // 2. raise rou to powers in the set
    todo!()
}

pub trait PartialNTT<
    const Q: u64,
    const N: usize,
    const D: usize,
    const Z: usize,
    const PHI_Z: usize,
>
{
    const ROOT_OF_UNITY: u64 = primitive_root_of_unity::<Q, Z>();
    const REL_POWS_ROOT_OF_UNITY: [Zq<Q>; PHI_Z] =
        root_of_unity_relative_powers::<Q, Z, PHI_Z>(Self::ROOT_OF_UNITY);

    fn ntt(a: &mut [Zq<Q>; N]) {
        let rou = Zq::<Q>::from(Self::ROOT_OF_UNITY);
        // 1. Split array in chunks of size D
        let components = coprimes_set::<Z, PHI_Z>().map(|power| {
            let rj = rou.pow([power as u64]);
        });
        // 2. reduce polynomial for each chunk with the corresponding relative power of the root of
        //    unity
        todo!()
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
