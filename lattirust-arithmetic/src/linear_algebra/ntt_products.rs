//! Exact batched products with an optional, statically selected ring backend.
//! Negacyclic products use auxiliary NTT primes and checked integer CRT
//! reconstruction. All temporary transforms are local to one call.
//!
//! Existing linear-algebra entry points select the backend without extra bounds:
//! ```
//! use lattirust_arithmetic::{
//!     linear_algebra::inner_products::dot_product,
//!     ring::{Pow2CyclotomicPolyRing, Ring, Zq1},
//! };
//! type P = Pow2CyclotomicPolyRing<Zq1<67280421310721>, 128>;
//! let row = vec![P::ONE; 128];
//! assert_eq!(dot_product(&row, &row), P::try_from(128u64).unwrap());
//! ```
use std::marker::PhantomData;

use ark_ff::PrimeField;
use rayon::prelude::*;

use crate::ring::{
    f_p::{fq_zero, Fq, FqConfig},
    ntt::{
        const_fq_from, const_inv_mod, intt_inplace_with_roots, ntt_inplace_with_roots,
        pows_bit_reversed, two_nth_primitive_root_of_unity,
    },
    Pow2CyclotomicPolyRing, Ring,
};

/// Optional acceleration for a ring. Returning `None` selects ordinary ring
/// arithmetic. Existing ring implementations need not supply a backend.
/// Implementations must return exact products in the input ring, in pair order.
pub trait BatchedProductBackend<R: Ring> {
    /// Whether the backend can accelerate batches in `R`.
    const SUPPORTED: bool = false;
    /// Whether it can accelerate coefficient-domain polynomials over `R`.
    const NEGACYCLIC_SUPPORTED: bool = false;

    fn products(
        _left: &[&[R]],
        _right: Option<&[&[R]]>,
        _pairs: &[(usize, usize)],
    ) -> Option<Vec<R>> {
        None
    }

    /// Extension point for coefficient rings that support polynomial products.
    fn negacyclic_products<const N: usize>(
        _left: &[&[Pow2CyclotomicPolyRing<R, N>]],
        _right: Option<&[&[Pow2CyclotomicPolyRing<R, N>]]>,
        _pairs: &[(usize, usize)],
    ) -> Option<Vec<Pow2CyclotomicPolyRing<R, N>>> {
        None
    }
}

/// Default backend: use the existing ring operations.
pub struct NoBatchedProducts;
impl<R: Ring> BatchedProductBackend<R> for NoBatchedProducts {}

/// A scalar residue ring with a canonical, single-word integer representation.
/// `MODULUS = None` disables acceleration. If it is `Some(q)`, canonical values
/// must be in `0..q`, and reconstruction must invert that representation.
pub trait IntegerCoefficient: Ring {
    const MODULUS: Option<u64>;
    fn canonical_u64(self) -> u64;
    fn from_canonical_u64(value: u64) -> Self {
        Self::try_from(value).expect("canonical coefficient must belong to the ring")
    }
}

/// Coefficient-ring backend for exact negacyclic products.
pub struct IntegerCoefficientBackend;
impl<C: IntegerCoefficient> BatchedProductBackend<C> for IntegerCoefficientBackend {
    // Prefer ordinary multiplication for word moduli that are cheap to reduce.
    // Smaller moduli can still use the explicit CRT entry point.
    const NEGACYCLIC_SUPPORTED: bool = matches!(C::MODULUS, Some(q) if q > u32::MAX as u64);

    fn negacyclic_products<const N: usize>(
        left: &[&[Pow2CyclotomicPolyRing<C, N>]],
        right: Option<&[&[Pow2CyclotomicPolyRing<C, N>]]>,
        pairs: &[(usize, usize)],
    ) -> Option<Vec<Pow2CyclotomicPolyRing<C, N>>> {
        compute_negacyclic_products(left, right, pairs, true)
    }
}

/// Forward a polynomial batch to the coefficient ring's backend.
pub struct NegacyclicProductBackend<C, const N: usize>(PhantomData<C>);
impl<C: Ring, const N: usize> BatchedProductBackend<Pow2CyclotomicPolyRing<C, N>>
    for NegacyclicProductBackend<C, N>
{
    // At small degrees, conversion and auxiliary transforms cost more than
    // the existing multiplication path. The explicit CRT API remains available.
    const SUPPORTED: bool =
        C::ProductBackend::NEGACYCLIC_SUPPORTED && supported_degree(N) && N >= 128;

    fn products(
        left: &[&[Pow2CyclotomicPolyRing<C, N>]],
        right: Option<&[&[Pow2CyclotomicPolyRing<C, N>]]>,
        pairs: &[(usize, usize)],
    ) -> Option<Vec<Pow2CyclotomicPolyRing<C, N>>> {
        if !Self::SUPPORTED {
            return None;
        }
        C::ProductBackend::negacyclic_products(left, right, pairs)
    }
}

const fn supported_degree(n: usize) -> bool {
    n >= 2 && n.is_power_of_two() && n <= (1 << 24)
}

pub(super) fn supports<R: Ring>() -> bool {
    R::ProductBackend::SUPPORTED
}

/// Try a batch of `sum_k left[i][k] * right[j][k]` for each requested `(i, j)`.
/// `right = None` reuses `left`. A backend may decline with `None`; callers
/// should then use ordinary ring products. Invalid indices or dimensions panic
/// when the backend accepts the batch.
pub fn products<R: Ring>(
    left: &[&[R]],
    right: Option<&[&[R]]>,
    pairs: &[(usize, usize)],
) -> Option<Vec<R>> {
    R::ProductBackend::products(left, right, pairs)
}

const P: u64 = 2_013_265_921;
const Q: u64 = 2_113_929_217;
const R: u64 = 1_811_939_329;
const S: u64 = 469_762_049;

// Invalid degrees must compile and decline acceleration. Guard constant
// generation as well as execution, so no invalid NTT is instantiated.
struct AuxiliaryNtt<const M: u64, const N: usize>;
impl<const M: u64, const N: usize> AuxiliaryNtt<M, N> {
    const ROOT: u64 = if supported_degree(N) {
        two_nth_primitive_root_of_unity::<M, N>()
    } else {
        1
    };
    const FORWARD: [Fq<M>; N] = if supported_degree(N) {
        pows_bit_reversed(Self::ROOT)
    } else {
        [fq_zero(); N]
    };
    const INVERSE: [Fq<M>; N] = if supported_degree(N) {
        pows_bit_reversed(const_inv_mod::<M>(Self::ROOT))
    } else {
        [fq_zero(); N]
    };
    const INVERSE_LENGTH: Fq<M> = if supported_degree(N) {
        const_fq_from(const_inv_mod::<M>(N as u64))
    } else {
        fq_zero()
    };

    fn forward(mut coefficients: [Fq<M>; N]) -> [Fq<M>; N] {
        ntt_inplace_with_roots(&mut coefficients, &Self::FORWARD);
        coefficients
    }
    fn inverse(mut values: [Fq<M>; N]) -> [Fq<M>; N] {
        intt_inplace_with_roots(&mut values, &Self::INVERSE, Self::INVERSE_LENGTH);
        values
    }
}

#[inline]
fn centered<C: IntegerCoefficient>(value: C) -> i64 {
    let modulus = C::MODULUS.expect("integer coefficient modulus");
    let residue = value.canonical_u64();
    if residue > modulus / 2 {
        -((modulus - residue) as i64)
    } else {
        residue as i64
    }
}

fn norm_max<C: IntegerCoefficient, const N: usize>(
    rows: &[&[Pow2CyclotomicPolyRing<C, N>]],
) -> u128 {
    rows.par_iter()
        .flat_map(|row| row.par_iter())
        .filter(|p| **p != Pow2CyclotomicPolyRing::ZERO)
        .map(|p| {
            p.coefficient_slice()
                .iter()
                .copied()
                .map(|c| centered(c).unsigned_abs())
                .max()
                .unwrap_or_default()
        })
        .max()
        .unwrap_or_default()
        .into()
}

type Transform<const N: usize> = [u32; N];

fn evals<C: IntegerCoefficient, const N: usize, const M: u64>(
    rows: &[&[Pow2CyclotomicPolyRing<C, N>]],
    start: usize,
    end: usize,
) -> Vec<Vec<Transform<N>>> {
    rows.par_iter()
        .map(|row| {
            row[start.min(row.len())..end.min(row.len())]
                .par_iter()
                .map(|p| {
                    if *p == Pow2CyclotomicPolyRing::ZERO {
                        return [0; N];
                    }
                    let input = p.coefficient_slice();
                    let coefficients = core::array::from_fn(|i| {
                        FqConfig::<M>::from_canonical_u64(
                            centered(input[i]).rem_euclid(M as i64) as u64
                        )
                        .expect("reduced auxiliary coefficient")
                    });
                    AuxiliaryNtt::<M, N>::forward(coefficients).map(|c| c.value as u32)
                })
                .collect()
        })
        .collect()
}

// Each residue is below 2^31. Even usize::MAX terms on a 64-bit target
// have sum below 2^126, so the fused accumulator fits in u128.
fn accumulate<const N: usize>(left: &[Transform<N>], right: &[Transform<N>]) -> [u128; N] {
    let mut sum = [0; N];
    for (a, b) in left.iter().zip(right) {
        for ((result, &x), &y) in sum.iter_mut().zip(a).zip(b) {
            *result += u128::from(u64::from(x) * u64::from(y));
        }
    }
    sum
}

fn prime_products<C: IntegerCoefficient, const N: usize, const M: u64>(
    left: &[&[Pow2CyclotomicPolyRing<C, N>]],
    right: Option<&[&[Pow2CyclotomicPolyRing<C, N>]]>,
    pairs: &[(usize, usize)],
) -> Vec<[u64; N]> {
    // Transform a column tile, use it for every requested dot, and release
    // it before the next tile or the next CRT prime.
    let width = pairs
        .iter()
        .map(|&(i, _)| left[i].len())
        .max()
        .unwrap_or_default();
    let mut sums = vec![[0u128; N]; pairs.len()];
    // Preserve the existing batch schedule and bound scratch by one tile.
    const TILE: usize = 512;
    for start in (0..width).step_by(TILE) {
        let end = (start + TILE).min(width);
        let a = evals::<C, N, M>(left, start, end);
        let b = right.map(|b| evals::<C, N, M>(b, start, end));
        let b = b.as_ref().unwrap_or(&a);
        sums.par_iter_mut().zip(pairs).for_each(|(out, &(i, j))| {
            let l = &a[i];
            let r = &b[j];
            let subtotal = if pairs.len() < rayon::current_num_threads() && l.len() >= 256 {
                l.par_chunks(128)
                    .zip(r.par_chunks(128))
                    .map(|(l, r)| accumulate(l, r))
                    .reduce(
                        || [0; N],
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
            AuxiliaryNtt::<M, N>::inverse(eval).map(|v| v.into_bigint().0[0])
        })
        .collect()
}

fn merge<const N: usize, const M: u64>(
    values: &mut [[u128; N]],
    residues: &[[u64; N]],
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

/// Exact auxiliary-prime products, with conservative batch-size thresholds.
/// Supports power-of-two degrees from 2 through 2^24 and single-word moduli.
/// Returns `None` if the batch is small or integer reconstruction is unsafe.
/// If `A` and `B` bound the centered coefficients, `w` is the largest dot
/// width, and `M` is the product of the selected auxiliary primes, require
/// `2 * w * N * A * B < M`. The check uses overflow-checked arithmetic.
pub fn negacyclic_products<C: IntegerCoefficient, const N: usize>(
    left: &[&[Pow2CyclotomicPolyRing<C, N>]],
    right: Option<&[&[Pow2CyclotomicPolyRing<C, N>]]>,
    pairs: &[(usize, usize)],
) -> Option<Vec<Pow2CyclotomicPolyRing<C, N>>> {
    compute_negacyclic_products(left, right, pairs, false)
}

fn compute_negacyclic_products<C: IntegerCoefficient, const N: usize>(
    left: &[&[Pow2CyclotomicPolyRing<C, N>]],
    right: Option<&[&[Pow2CyclotomicPolyRing<C, N>]]>,
    pairs: &[(usize, usize)],
    automatic: bool,
) -> Option<Vec<Pow2CyclotomicPolyRing<C, N>>> {
    if !supported_degree(N) || !matches!(C::MODULUS, Some(q) if q > 1) || pairs.is_empty() {
        return None;
    }
    let width = pairs.iter().map(|&(i, _)| left[i].len()).max()?;
    if width < 16 || width.saturating_mul(pairs.len()) < 128 {
        return None;
    }
    let a = left;
    let b = right;
    let other = b.unwrap_or(a);
    assert!(
        pairs.iter().all(|&(i, j)| a[i].len() == other[j].len()),
        "inner product dimension mismatch"
    );
    // Decline low-reuse batches over wide coefficient moduli before scanning
    // input norms. The explicit entry point bypasses this conservative policy.
    let families = left
        .len()
        .saturating_add(right.map_or(0, |rows| rows.len()));
    if automatic && C::MODULUS.unwrap() > (1 << 40) && pairs.len() < families.saturating_mul(8) {
        return None;
    }
    let norm_a = norm_max(a);
    let norm_b = b.map_or(norm_a, norm_max);
    let bound = (width as u128)
        .checked_mul(N as u128)?
        .checked_mul(norm_a)?
        .checked_mul(norm_b)?
        .checked_mul(2)?;
    if bound == 0 {
        return Some(vec![Pow2CyclotomicPolyRing::ZERO; pairs.len()]);
    }
    const PQ: u128 = P as u128 * Q as u128;
    const PQR: u128 = PQ * R as u128;
    const PQRS: u128 = PQR * S as u128;
    if bound >= PQRS {
        return None;
    }
    let mut values: Vec<[u128; N]> = prime_products::<C, N, P>(a, b, pairs)
        .into_iter()
        .map(|row| row.map(u128::from))
        .collect();
    let mut product = P as u128;
    if bound >= product {
        let fresh = prime_products::<C, N, Q>(a, b, pairs);
        merge::<N, Q>(&mut values, &fresh, product, const_inv_mod::<Q>(P % Q));
        product = PQ;
    }
    if bound >= product {
        let fresh = prime_products::<C, N, R>(a, b, pairs);
        merge::<N, R>(
            &mut values,
            &fresh,
            product,
            const_inv_mod::<R>((PQ % R as u128) as u64),
        );
        product = PQR;
    }
    if bound >= product {
        let fresh = prime_products::<C, N, S>(a, b, pairs);
        merge::<N, S>(
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
                Pow2CyclotomicPolyRing::from_fn(|i| {
                    let value = if row[i] > product / 2 {
                        row[i] as i128 - product as i128
                    } else {
                        row[i] as i128
                    };
                    C::from_canonical_u64(value.rem_euclid(i128::from(C::MODULUS.unwrap())) as u64)
                })
            })
            .collect(),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ring::Zq1;
    const MODULUS: u64 = 274_877_906_837;
    type Coefficient = Zq1<MODULUS>;
    type Polynomial = Pow2CyclotomicPolyRing<Coefficient, 256>;
    use ark_std::{
        rand::{rngs::StdRng, SeedableRng},
        UniformRand,
    };

    #[test]
    fn reconstructs_signed_dots_for_small_and_full_width_coefficients() {
        let mut rng = StdRng::seed_from_u64(734);
        let small: Vec<_> = (0..529)
            .map(|i| {
                Polynomial::from_fn(|k| {
                    if (i + k) % 7 == 0 {
                        -Coefficient::ONE
                    } else {
                        Coefficient::from((i + k) % 11 == 0)
                    }
                })
            })
            .collect();
        let large: Vec<_> = (0..529).map(|_| Polynomial::rand(&mut rng)).collect();
        let rows = [small.as_slice(), large.as_slice()];
        let pairs = [(0, 0), (1, 0), (1, 1)];
        let expected: Vec<Polynomial> = pairs
            .iter()
            .map(|&(i, j)| rows[i].iter().zip(rows[j]).map(|(&a, &b)| a * b).sum())
            .collect();
        assert_eq!(products(&rows, None, &pairs).unwrap(), expected);
        let left = [small.as_slice()];
        assert_eq!(products(&left, None, &[(0, 0)]).unwrap(), vec![expected[0]]);
        let right = [large.as_slice()];
        assert_eq!(
            products(&left, Some(&right), &[(0, 0)]).unwrap(),
            vec![expected[1]]
        );
    }
    #[test]
    fn reconstructs_a_wide_dot_with_four_primes() {
        let high = Coefficient::try_from(MODULUS / 2).unwrap();
        let x = Polynomial::from_fn(|_| high);
        let y = -x;
        let left = vec![x; 1031];
        let right = vec![y; 1031];
        let actual = products(&[&left], Some(&[&right]), &[(0, 0)]).unwrap();
        let scale = Coefficient::try_from(1031u64).unwrap();
        // For (1 + ... + X^255)^2 modulo X^256+1, coefficient k
        // is 2(k+1)-256. The largest signed integer coefficient exceeds
        // half the three-prime product, so the fourth prime is necessary.
        let expected = Polynomial::from_fn(|k| {
            let factor = Coefficient::try_from((2 * (k + 1)) as u64).unwrap()
                - Coefficient::try_from(256u64).unwrap();
            -(factor * high * high * scale)
        });
        assert_eq!(actual, vec![expected]);
    }
    fn schoolbook<C: IntegerCoefficient, const N: usize>(
        left: &[Pow2CyclotomicPolyRing<C, N>],
        right: &[Pow2CyclotomicPolyRing<C, N>],
    ) -> Pow2CyclotomicPolyRing<C, N> {
        assert_eq!(left.len(), right.len());
        let mut result = [C::ZERO; N];
        for (a, b) in left.iter().zip(right) {
            for (i, &x) in a.coefficient_slice().iter().enumerate() {
                for (j, &y) in b.coefficient_slice().iter().enumerate() {
                    if i + j < N {
                        result[i + j] += x * y;
                    } else {
                        result[i + j - N] -= x * y;
                    }
                }
            }
        }
        Pow2CyclotomicPolyRing::from(result)
    }

    fn check_generic<const Q: u64, const N: usize>() {
        type C<const Q: u64> = Zq1<Q>;
        let mut rng = StdRng::seed_from_u64(972);
        let rows: Vec<Vec<Pow2CyclotomicPolyRing<C<Q>, N>>> = (0..3)
            .map(|_| {
                (0..129)
                    .map(|_| Pow2CyclotomicPolyRing::rand(&mut rng))
                    .collect()
            })
            .collect();
        let left = [&rows[0][..128], &rows[1][..129]];
        let right = [&rows[2][..128], &rows[2][..129]];
        let pairs = [(0, 0), (1, 1), (0, 0)];
        let expected: Vec<_> = pairs
            .iter()
            .map(|&(i, j)| schoolbook(left[i], right[j]))
            .collect();
        assert_eq!(
            negacyclic_products(&left, Some(&right), &pairs).unwrap(),
            expected
        );
        let paired = [&rows[0][..], &rows[1][..]];
        let expected: Vec<_> = [(0, 0), (1, 0), (1, 1)]
            .iter()
            .map(|&(i, j)| schoolbook(paired[i], paired[j]))
            .collect();
        assert_eq!(
            negacyclic_products(&paired, None, &[(0, 0), (1, 0), (1, 1)]).unwrap(),
            expected
        );
    }

    #[test]
    fn generic_degrees_moduli_and_prime_counts_match_schoolbook() {
        check_generic::<17, 8>(); // One auxiliary prime.
        check_generic::<65537, 16>(); // Two auxiliary primes.
        check_generic::<274877906837, 32>(); // Three auxiliary primes.
        check_generic::<3329, 64>();
    }

    #[test]
    fn declines_unsupported_degrees_and_excessive_reconstruction_bounds() {
        type Wide = Pow2CyclotomicPolyRing<Zq1<18446744069414584321>, 8>;
        let high = Wide::from_fn(|_| Zq1::try_from(9223372034707292160u64).unwrap());
        let rows = vec![high; 128];
        assert!(negacyclic_products(&[&rows], None, &[(0, 0)]).is_none());
        // This bound fits in u128 but exceeds the product of all four primes.
        type Bounded = Pow2CyclotomicPolyRing<Zq1<2305843009213693951>, 2>;
        let high = Bounded::from_fn(|_| Zq1::try_from(1152921504606846975u64).unwrap());
        let bounded = vec![high; 16];
        assert!(negacyclic_products(&[&bounded], None, &[(0, 0); 8]).is_none());
        // A large target modulus remains supported when the actual lift is small.
        let negative = Wide::from_fn(|_| -Zq1::ONE);
        let small = vec![negative; 128];
        assert_eq!(
            negacyclic_products(&[&small], None, &[(0, 0)]).unwrap(),
            vec![schoolbook(&small, &small)]
        );
        // Degree three is representable by the existing coefficient type;
        // it must fall back without instantiating an invalid NTT.
        type Odd = Pow2CyclotomicPolyRing<Zq1<17>, 3>;
        let odd = vec![Odd::ONE; 128];
        assert!(products(&[&odd], None, &[(0, 0)]).is_none());
        type Constant = Pow2CyclotomicPolyRing<Zq1<17>, 1>;
        let constant = vec![Constant::ONE; 128];
        assert!(products(&[&constant], None, &[(0, 0)]).is_none());
    }

    #[test]
    fn zero_empty_and_small_batches_preserve_fallbacks() {
        let zeros = vec![Polynomial::ZERO; 128];
        assert_eq!(
            products(&[&zeros], None, &[(0, 0)]),
            Some(vec![Polynomial::ZERO])
        );
        assert!(products::<Polynomial>(&[], None, &[]).is_none());
        assert!(products(&[&zeros[..15]], None, &[(0, 0)]).is_none());
        assert!(products(&[&zeros[..16]], None, &[(0, 0)]).is_none());
        assert_eq!(
            crate::linear_algebra::inner_products::dot_product(&zeros[..7], &zeros[..7]),
            Polynomial::ZERO
        );
        assert!(!supports::<crate::ring::Z2>());
        type Small = Pow2CyclotomicPolyRing<Zq1<65537>, 64>;
        let small = vec![Small::ONE; 128];
        assert!(products(&[&small], None, &[(0, 0)]).is_none());
        type SmallModulus = Pow2CyclotomicPolyRing<Zq1<65537>, 256>;
        let small_modulus = vec![SmallModulus::ONE; 128];
        assert!(products(&[&small_modulus], None, &[(0, 0)]).is_none());
        assert_eq!(
            negacyclic_products(&[&small], None, &[(0, 0)]).unwrap(),
            vec![Small::try_from(128u64).unwrap()]
        );
    }

    #[test]
    #[should_panic(expected = "inner product dimension mismatch")]
    fn rejects_incompatible_row_lengths() {
        let left = vec![Polynomial::ONE; 128];
        let right = vec![Polynomial::ONE; 129];
        let _ = products(&[&left], Some(&[&right]), &[(0, 0)]);
    }

    #[test]
    fn existing_linear_algebra_interfaces_preserve_results_across_dispatch_and_fallback() {
        use crate::linear_algebra::{inner_products::*, Matrix, Vector};
        type C = Zq1<67280421310721>;
        type P = Pow2CyclotomicPolyRing<C, 128>;
        assert!(supports::<P>());
        let mut rng = StdRng::seed_from_u64(143);
        let s: Vec<Vector<P>> = (0..3).map(|_| Vector::rand(129, &mut rng)).collect();
        let t: Vec<Vector<P>> = (0..3).map(|_| Vector::rand(129, &mut rng)).collect();
        let gram = inner_products(&s);
        let cross = inner_products2(&s, &t);
        let symmetric = symmetric_cross_products(&s, &t);
        for i in 0..3 {
            assert_eq!(
                dot_product(s[i].as_slice(), t[i].as_slice()),
                schoolbook(s[i].as_slice(), t[i].as_slice())
            );
            for j in 0..=i {
                let st = schoolbook(s[i].as_slice(), t[j].as_slice());
                assert_eq!(*gram.at(i, j), schoolbook(s[i].as_slice(), s[j].as_slice()));
                assert_eq!(*cross.at(i, j), st);
                assert_eq!(
                    *symmetric.at(i, j),
                    if i == j {
                        st
                    } else {
                        st + schoolbook(s[j].as_slice(), t[i].as_slice())
                    }
                );
            }
        }
        let matrix = Matrix::from_fn(3, 129, |i, j| s[i][j]);
        let result = matrix_vector_product(&matrix, &t[0]);
        for i in 0..3 {
            assert_eq!(result[i], schoolbook(s[i].as_slice(), t[0].as_slice()));
        }
    }
}
