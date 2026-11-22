//! Reference implementation retained for before/after transform benchmarks.
use ark_std::Zero;

use lattirust_arithmetic::ring::{Ring, Zq1};

pub trait PartialNTT<
    const Q: u64,
    const N: usize,     // Ring degree
    const D: usize,     // Degree of each split factor x^D - r_j
    const Z: usize,     // Split parameter, with q = 1 mod Z
    const PHI_Z: usize, // Number of split factors
>
{
    fn ntt(a: &mut [Zq1<Q>; N], primitive_zth_root: Zq1<Q>) {
        debug_assert_eq!(N, D * PHI_Z);
        debug_assert_eq!(Q % Z as u64, 1);

        let coeffs = *a;
        *a = [Zq1::<Q>::ZERO; N];
        let roots = split_roots::<Q, Z, PHI_Z>(primitive_zth_root);

        for (component_idx, root) in roots.into_iter().enumerate() {
            let mut root_power = Zq1::<Q>::ONE;
            for block in 0..PHI_Z {
                for local_idx in 0..D {
                    a[component_idx * D + local_idx] += coeffs[block * D + local_idx] * root_power;
                }
                root_power *= root;
            }
        }
    }

    fn intt(a: &mut [Zq1<Q>; N], primitive_zth_root: Zq1<Q>) {
        debug_assert_eq!(N, D * PHI_Z);
        debug_assert_eq!(Q % Z as u64, 1);

        let evals = *a;
        let inverse_vandermonde =
            inverse_vandermonde::<Q, PHI_Z>(&vandermonde::<Q, Z, PHI_Z>(primitive_zth_root));

        let mut coeffs = [Zq1::<Q>::ZERO; N];
        for local_idx in 0..D {
            let mut values = [Zq1::<Q>::ZERO; PHI_Z];
            for component_idx in 0..PHI_Z {
                values[component_idx] = evals[component_idx * D + local_idx];
            }

            for block in 0..PHI_Z {
                let mut coeff = Zq1::<Q>::ZERO;
                for component_idx in 0..PHI_Z {
                    coeff += inverse_vandermonde[block][component_idx] * values[component_idx];
                }
                coeffs[block * D + local_idx] = coeff;
            }
        }

        *a = coeffs;
    }
}

pub(crate) const fn coprimes_set<const Z: usize, const PHI_Z: usize>() -> [usize; PHI_Z] {
    let mut i = 0;
    let mut j = 0;
    let mut set = [0; PHI_Z];
    while i < Z {
        if is_coprime_to::<Z>(i) {
            set[j] = i;
            j += 1;
        }
        i += 1;
    }
    set
}

pub(crate) fn split_roots<const Q: u64, const Z: usize, const PHI_Z: usize>(
    primitive_zth_root: Zq1<Q>,
) -> [Zq1<Q>; PHI_Z] {
    coprimes_set::<Z, PHI_Z>().map(|power| primitive_zth_root.pow(power as u64))
}

pub(crate) fn vandermonde<const Q: u64, const Z: usize, const PHI_Z: usize>(
    primitive_zth_root: Zq1<Q>,
) -> [[Zq1<Q>; PHI_Z]; PHI_Z] {
    let roots = split_roots::<Q, Z, PHI_Z>(primitive_zth_root);
    core::array::from_fn(|row| {
        let root = roots[row];
        let mut powers = [Zq1::<Q>::ZERO; PHI_Z];
        let mut current = Zq1::<Q>::ONE;
        let mut col = 0;
        while col < PHI_Z {
            powers[col] = current;
            current *= root;
            col += 1;
        }
        powers
    })
}

pub(crate) fn inverse_vandermonde<const Q: u64, const PHI_Z: usize>(
    matrix: &[[Zq1<Q>; PHI_Z]; PHI_Z],
) -> [[Zq1<Q>; PHI_Z]; PHI_Z] {
    let mut left = *matrix;
    let mut right =
        core::array::from_fn(|row| core::array::from_fn(|col| Zq1::<Q>::from(row == col)));

    for pivot_idx in 0..PHI_Z {
        if left[pivot_idx][pivot_idx].is_zero() {
            let swap_idx = (pivot_idx + 1..PHI_Z)
                .find(|&row| !left[row][pivot_idx].is_zero())
                .expect("split Vandermonde matrix is singular");
            left.swap(pivot_idx, swap_idx);
            right.swap(pivot_idx, swap_idx);
        }

        let pivot_inv = left[pivot_idx][pivot_idx]
            .inverse()
            .expect("non-zero pivot must be invertible");
        for col in 0..PHI_Z {
            left[pivot_idx][col] *= pivot_inv;
            right[pivot_idx][col] *= pivot_inv;
        }

        for row in 0..PHI_Z {
            if row == pivot_idx {
                continue;
            }
            let factor = left[row][pivot_idx];
            if factor.is_zero() {
                continue;
            }

            for col in 0..PHI_Z {
                left[row][col] -= factor * left[pivot_idx][col];
                right[row][col] -= factor * right[pivot_idx][col];
            }
        }
    }

    right
}

const fn is_coprime_to<const Z: usize>(i: usize) -> bool {
    let mut a = i;
    let mut b = Z;
    while b != 0 {
        let temp = b;
        b = a % b;
        a = temp;
    }
    a == 1
}
