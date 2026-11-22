//! Exact batched negacyclic products for the degree-256, 38-bit proof ring.
//! A batch shares transformed rows. Auxiliary CRT primes are arithmetic
//! scratch space; the result and its representation stay in the proof ring.
use std::any::{Any, TypeId};

use ark_ff::PrimeField;
use rayon::prelude::*;

use lattirust_arithmetic::ring::{
    f_p::{Fq, FqConfig},
    ntt::{const_inv_mod, Ntt},
    Pow2CyclotomicPolyRing, Ring, Zq1,
};

const DEGREE: usize = 256;
const MODULUS: u64 = 274_877_906_837;
type Coefficient = Zq1<MODULUS>;
type Polynomial = Pow2CyclotomicPolyRing<Coefficient, DEGREE>;

const P: u64 = 2_013_265_921;
const Q: u64 = 2_113_929_217;
const R: u64 = 1_811_939_329;
const S: u64 = 469_762_049;

fn centered(value: Coefficient) -> i64 {
    use lattirust_arithmetic::ring::{Zq1ConfigImpl, ZqConfig};
    let residue = Zq1ConfigImpl::<MODULUS>::into_bigint(value).0[0];
    if residue > MODULUS / 2 {
        residue as i64 - MODULUS as i64
    } else {
        residue as i64
    }
}

fn norm_max(rows: &[Vec<&Polynomial>]) -> u128 {
    rows.par_iter()
        .flat_map(|row| row.par_iter())
        .filter(|&&p| *p != Polynomial::ZERO)
        .map(|p| {
            reference_coefficients(p)
                .into_iter()
                .map(|c| centered(c).unsigned_abs())
                .max()
                .unwrap_or_default()
        })
        .max()
        .unwrap_or_default()
        .into()
}

type Transform = [u32; DEGREE];

fn evals<const M: u64>(rows: &[Vec<&Polynomial>], start: usize, end: usize) -> Vec<Vec<Transform>> {
    rows.par_iter()
        .map(|row| {
            row[start.min(row.len())..end.min(row.len())]
                .par_iter()
                .map(|&p| {
                    if *p == Polynomial::ZERO {
                        return [0; DEGREE];
                    }
                    let coefficients = reference_coefficients(p)
                        .map(|c| FqConfig::<M>::from_u64(centered(c).rem_euclid(M as i64) as u64));
                    Fq::<M>::ntt(coefficients).map(|c| c.value as u32)
                })
                .collect()
        })
        .collect()
}

fn accumulate(left: &[Transform], right: &[Transform]) -> [u128; DEGREE] {
    let mut sum = [0; DEGREE];
    for (a, b) in left.iter().zip(right) {
        for ((result, &x), &y) in sum.iter_mut().zip(a).zip(b) {
            *result += u128::from(u64::from(x) * u64::from(y));
        }
    }
    sum
}

fn prime_products<const M: u64>(
    left: &[Vec<&Polynomial>],
    right: Option<&[Vec<&Polynomial>]>,
    pairs: &[(usize, usize)],
) -> Vec<[u64; DEGREE]> {
    // Transform a column tile, use it for every requested dot, and release
    // it before the next tile or the next CRT prime.
    let width = pairs
        .iter()
        .map(|&(i, _)| left[i].len())
        .max()
        .unwrap_or_default();
    let mut sums = vec![[0u128; DEGREE]; pairs.len()];
    const TILE: usize = 512;
    for start in (0..width).step_by(TILE) {
        let end = (start + TILE).min(width);
        let a = evals::<M>(left, start, end);
        let b = right.map(|b| evals::<M>(b, start, end));
        let b = b.as_ref().unwrap_or(&a);
        sums.par_iter_mut().zip(pairs).for_each(|(out, &(i, j))| {
            let l = &a[i];
            let r = &b[j];
            let subtotal = if pairs.len() < rayon::current_num_threads() && l.len() >= 256 {
                l.par_chunks(128)
                    .zip(r.par_chunks(128))
                    .map(|(l, r)| accumulate(l, r))
                    .reduce(
                        || [0; DEGREE],
                        |mut l, r| {
                            for (x, y) in l.iter_mut().zip(r) {
                                *x += y;
                            }
                            l
                        },
                    )
            } else {
                accumulate(l, r)
            };
            for (x, y) in out.iter_mut().zip(subtotal) {
                *x += y;
            }
        });
    }
    sums.into_par_iter()
        .map(|sum| {
            // Raw inputs carry one Montgomery factor each. Remove one
            // factor before interpreting the fused dot as an NTT value.
            let eval =
                sum.map(|v| Fq::<M>::from_raw((v % u128::from(M)) as u64) * Fq::<M>::from_raw(1));
            Fq::<M>::intt(eval).map(|v| v.into_bigint().0[0])
        })
        .collect()
}

fn merge<const M: u64>(
    values: &mut [[u128; DEGREE]],
    residues: &[[u64; DEGREE]],
    product: u128,
    inverse: u64,
) {
    for (row, fresh) in values.iter_mut().zip(residues) {
        for (x, &y) in row.iter_mut().zip(fresh) {
            let digit = (y + M - (*x % u128::from(M)) as u64) * inverse % M;
            *x += product * u128::from(digit);
        }
    }
}

pub(super) fn supports<T: Ring>() -> bool {
    TypeId::of::<T>() == TypeId::of::<Polynomial>()
}

fn borrow_rows<'a, T: Ring>(rows: &[&'a [T]]) -> Vec<Vec<&'a Polynomial>> {
    rows.iter()
        .map(|row| {
            row.iter()
                .map(|p| (p as &dyn Any).downcast_ref::<Polynomial>().unwrap())
                .collect()
        })
        .collect()
}

/// Try the accelerated backend and return None for other rings or small batches.
pub(super) fn products<T: Ring>(
    left: &[&[T]],
    right: Option<&[&[T]]>,
    pairs: &[(usize, usize)],
) -> Option<Vec<T>> {
    if !supports::<T>() || pairs.is_empty() {
        return None;
    }
    let width = pairs.iter().map(|&(i, _)| left[i].len()).max()?;
    if width < 16 || width.saturating_mul(pairs.len()) < 128 {
        return None;
    }
    let a = borrow_rows(left);
    let b = right.map(borrow_rows);
    let other = b.as_ref().unwrap_or(&a);
    assert!(pairs.iter().all(|&(i, j)| a[i].len() == other[j].len()));
    let norm_a = norm_max(&a);
    let norm_b = b.as_ref().map_or(norm_a, |b| norm_max(b));
    let bound = (width as u128)
        .checked_mul(DEGREE as u128)?
        .checked_mul(norm_a)?
        .checked_mul(norm_b)?
        .checked_mul(2)?;
    if bound == 0 {
        return Some(vec![T::ZERO; pairs.len()]);
    }
    const PQ: u128 = P as u128 * Q as u128;
    const PQR: u128 = PQ * R as u128;
    const PQRS: u128 = PQR * S as u128;
    if bound >= PQRS {
        return None;
    }
    let mut values: Vec<[u128; DEGREE]> = prime_products::<P>(&a, b.as_deref(), pairs)
        .into_iter()
        .map(|row| row.map(u128::from))
        .collect();
    let mut product = P as u128;
    if bound >= product {
        let fresh = prime_products::<Q>(&a, b.as_deref(), pairs);
        merge::<Q>(&mut values, &fresh, product, const_inv_mod::<Q>(P % Q));
        product = PQ;
    }
    if bound >= product {
        let fresh = prime_products::<R>(&a, b.as_deref(), pairs);
        merge::<R>(
            &mut values,
            &fresh,
            product,
            const_inv_mod::<R>((PQ % R as u128) as u64),
        );
        product = PQR;
    }
    if bound >= product {
        let fresh = prime_products::<S>(&a, b.as_deref(), pairs);
        merge::<S>(
            &mut values,
            &fresh,
            product,
            const_inv_mod::<S>((PQR % S as u128) as u64),
        );
        product = PQRS;
    }
    Some(
        values
            .into_iter()
            .map(|row| {
                let polynomial = Polynomial::from_fn(|i| {
                    let value = if row[i] > product / 2 {
                        row[i] as i128 - product as i128
                    } else {
                        row[i] as i128
                    };
                    Coefficient::try_from(value.rem_euclid(i128::from(MODULUS)) as u64).unwrap()
                });
                *(&polynomial as &dyn Any).downcast_ref::<T>().unwrap()
            })
            .collect(),
    )
}

// The original method is crate-private. Reproduce its array copy through
// the existing public, allocation-free conversion into an SVector.
#[inline(always)]
fn reference_coefficients(p: &Polynomial) -> [Coefficient; DEGREE] {
    let coefficients: lattirust_arithmetic::linear_algebra::SVector<Coefficient, DEGREE> =
        (*p).into();
    core::array::from_fn(|i| coefficients[i])
}
