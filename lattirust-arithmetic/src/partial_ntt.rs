use ark_std::Zero;
use std::fmt::{Display, Formatter};

use crate::ring::{Ring, Zq1};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum PartialNttError {
    InvalidParameters,
    InvalidRoot,
}

impl Display for PartialNttError {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Self::InvalidParameters => "invalid power-of-two split NTT parameters",
            Self::InvalidRoot => "root does not have the required order",
        })
    }
}

impl std::error::Error for PartialNttError {}

/// Public transform parameters for power-of-two split orders.
/// The factorization is `X^N + 1 = product_j (X^D - r_j)`.
///
/// The factors retain the existing order `r_j = root^(2*j + 1)`.
/// Forward transforms twist the coefficients by `root^k`, then apply a
/// radix-two DFT across the factors. Inverse transforms undo both steps.
/// The plan stores O(PHI_Z) field elements and contains no polynomial data.
/// Construct it once and pass a shared reference to each transform.
/// Existing root-taking APIs also support non-power-of-two orders through
/// their original dense transform; this plan requires a power-of-two order.
///
/// ```
/// use lattirust_arithmetic::{partial_ntt::PartialNttPlan,
///     ring::{CyclotomicPolyRingSplittedNTT, Ring, Zq1}};
/// type Split = CyclotomicPolyRingSplittedNTT<3329, 256, 2, 256, 128>;
/// let plan = PartialNttPlan::new(Split::primitive_zth_root()).unwrap();
/// let coefficients = [Zq1::<3329>::ONE; 256];
/// let transformed = Split::from_coefficient_array_with_plan(coefficients, &plan);
/// assert_eq!(transformed.coefficients_with_plan(&plan), coefficients);
/// ```
#[derive(Clone, Debug)]
pub struct PartialNttPlan<
    const Q: u64,
    const N: usize,
    const D: usize,
    const Z: usize,
    const PHI_Z: usize,
> {
    primitive_zth_root: Zq1<Q>,
    roots: [Zq1<Q>; PHI_Z],
    twists: [Zq1<Q>; PHI_Z],
    inverse_scales: [Zq1<Q>; PHI_Z],
    twiddles: [Zq1<Q>; PHI_Z],
    inverse_twiddles: [Zq1<Q>; PHI_Z],
    bit_reversed: [usize; PHI_Z],
}

impl<const Q: u64, const N: usize, const D: usize, const Z: usize, const PHI_Z: usize>
    PartialNttPlan<Q, N, D, Z, PHI_Z>
{
    pub fn new(primitive_zth_root: Zq1<Q>) -> Result<Self, PartialNttError> {
        if Z < 2
            || !Z.is_power_of_two()
            || PHI_Z != Z / 2
            || D == 0
            || D.checked_mul(PHI_Z) != Some(N)
            || Q % Z as u64 != 1
        {
            return Err(PartialNttError::InvalidParameters);
        }
        if primitive_zth_root.pow(Z as u64) != Zq1::<Q>::ONE
            || primitive_zth_root.pow((Z / 2) as u64) != -Zq1::<Q>::ONE
        {
            return Err(PartialNttError::InvalidRoot);
        }
        let inverse_root = primitive_zth_root
            .inverse()
            .ok_or(PartialNttError::InvalidRoot)?;
        let inverse_length = Zq1::<Q>::try_from(PHI_Z as u64)
            .ok()
            .and_then(|length| length.inverse())
            .ok_or(PartialNttError::InvalidParameters)?;
        let powers = |base| {
            let mut value = Zq1::<Q>::ONE;
            core::array::from_fn(|_| {
                let current = value;
                value *= base;
                current
            })
        };
        let twists = powers(primitive_zth_root);
        let twiddles = powers(primitive_zth_root * primitive_zth_root);
        let inverse_twiddles = powers(inverse_root * inverse_root);
        let inverse_scales = powers(inverse_root).map(|value| value * inverse_length);
        let roots = twiddles.map(|value| primitive_zth_root * value);
        let bits = PHI_Z.trailing_zeros();
        let bit_reversed = core::array::from_fn(|i| {
            if bits == 0 {
                0
            } else {
                i.reverse_bits() >> (usize::BITS - bits)
            }
        });
        Ok(Self {
            primitive_zth_root,
            roots,
            twists,
            inverse_scales,
            twiddles,
            inverse_twiddles,
            bit_reversed,
        })
    }

    pub fn primitive_zth_root(&self) -> Zq1<Q> {
        self.primitive_zth_root
    }

    pub fn roots(&self) -> &[Zq1<Q>; PHI_Z] {
        &self.roots
    }

    pub fn ntt(&self, coefficients: &mut [Zq1<Q>; N]) {
        self.scale(coefficients, &self.twists);
        self.butterflies(coefficients, &self.twiddles);
    }

    pub fn intt(&self, values: &mut [Zq1<Q>; N]) {
        self.butterflies(values, &self.inverse_twiddles);
        self.scale(values, &self.inverse_scales);
    }

    fn scale(&self, values: &mut [Zq1<Q>; N], factors: &[Zq1<Q>; PHI_Z]) {
        for (block, &factor) in values.chunks_exact_mut(D).zip(factors) {
            for value in block {
                *value *= factor;
            }
        }
    }

    fn butterflies(&self, values: &mut [Zq1<Q>; N], twiddles: &[Zq1<Q>; PHI_Z]) {
        for (i, &j) in self.bit_reversed.iter().enumerate() {
            if i < j {
                for local in 0..D {
                    values.swap(i * D + local, j * D + local);
                }
            }
        }
        let mut length = 2;
        while length <= PHI_Z {
            let half = length / 2;
            let stride = PHI_Z / length;
            for start in (0..PHI_Z).step_by(length) {
                for j in 0..half {
                    let weight = twiddles[j * stride];
                    let left = (start + j) * D;
                    let right = (start + j + half) * D;
                    for local in 0..D {
                        let a = values[left + local];
                        let b = values[right + local] * weight;
                        values[left + local] = a + b;
                        values[right + local] = a - b;
                    }
                }
            }
            length *= 2;
        }
    }
}

pub trait PartialNTT<
    const Q: u64,
    const N: usize,
    const D: usize,
    const Z: usize,
    const PHI_Z: usize,
>
{
    fn ntt(a: &mut [Zq1<Q>; N], primitive_zth_root: Zq1<Q>) {
        if Z.is_power_of_two() {
            PartialNttPlan::<Q, N, D, Z, PHI_Z>::new(primitive_zth_root)
                .expect("invalid split NTT parameters or root")
                .ntt(a);
        } else {
            debug_assert_eq!(N, D * PHI_Z);
            debug_assert_eq!(Q % Z as u64, 1);

            let coeffs = *a;
            *a = [Zq1::<Q>::ZERO; N];
            let roots = split_roots::<Q, Z, PHI_Z>(primitive_zth_root);

            for (component_idx, root) in roots.into_iter().enumerate() {
                let mut root_power = Zq1::<Q>::ONE;
                for block in 0..PHI_Z {
                    for local_idx in 0..D {
                        a[component_idx * D + local_idx] +=
                            coeffs[block * D + local_idx] * root_power;
                    }
                    root_power *= root;
                }
            }
        }
    }

    fn intt(a: &mut [Zq1<Q>; N], primitive_zth_root: Zq1<Q>) {
        if Z.is_power_of_two() {
            PartialNttPlan::<Q, N, D, Z, PHI_Z>::new(primitive_zth_root)
                .expect("invalid split NTT parameters or root")
                .intt(a);
        } else {
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ring::{CyclotomicPolyRingSplittedNTT, Pow2CyclotomicPolyRing};
    use ark_std::{test_rng, UniformRand};

    fn check<const Q: u64, const N: usize, const D: usize, const Z: usize, const F: usize>() {
        let root = CyclotomicPolyRingSplittedNTT::<Q, N, D, Z, F>::primitive_zth_root();
        let rng = &mut test_rng();
        for root in [root, root.pow(3u64)] {
            let plan = PartialNttPlan::<Q, N, D, Z, F>::new(root).unwrap();
            assert_eq!(plan.primitive_zth_root(), root);
            for (j, &factor) in plan.roots().iter().enumerate() {
                assert_eq!(factor, root.pow((2 * j + 1) as u64));
            }
            for _ in 0..4 {
                let a = Pow2CyclotomicPolyRing::<Zq1<Q>, N>::rand(rng);
                let b = Pow2CyclotomicPolyRing::<Zq1<Q>, N>::rand(rng);
                let original = a.coefficient_array();
                let mut transformed = original;
                plan.ntt(&mut transformed);
                // Polynomial remainders give an independent definition of
                // the transform and its existing factor ordering.
                for (j, &factor) in plan.roots().iter().enumerate() {
                    for local in 0..D {
                        let mut expected = Zq1::<Q>::ZERO;
                        for block in (0..F).rev() {
                            expected = expected * factor + original[block * D + local];
                        }
                        assert_eq!(transformed[j * D + local], expected);
                    }
                }
                plan.intt(&mut transformed);
                assert_eq!(transformed, original);
                let left = CyclotomicPolyRingSplittedNTT::from_coefficient_array_with_plan(
                    original, &plan,
                );
                let right = CyclotomicPolyRingSplittedNTT::from_coefficient_array_with_plan(
                    b.coefficient_array(),
                    &plan,
                );
                let legacy = CyclotomicPolyRingSplittedNTT::<Q, N, D, Z, F>::from_coefficient_array(
                    original, root,
                );
                assert_eq!(left, legacy);
                assert_eq!(legacy.coefficients(root), original);
                let product = left.mul_with_plan(&right, &plan);
                assert_eq!(product, left.mul_with_root(&right, root));
                assert_eq!(
                    product.coefficients_with_plan(&plan),
                    (a * b).coefficient_array()
                );
            }
        }
    }

    #[test]
    fn plan_matches_factor_remainders_and_ring_products() {
        check::<17, 16, 8, 4, 2>();
        check::<3329, 256, 2, 256, 128>();
        check::<12289, 256, 1, 512, 256>();
        check::<274877906837, 256, 128, 4, 2>();
        check::<17, 4, 4, 2, 1>();
    }

    #[test]
    fn root_api_preserves_general_split_orders() {
        type Split = CyclotomicPolyRingSplittedNTT<19, 6, 3, 3, 2>;
        type F = Zq1<19>;
        let root = F::try_from(7u64).unwrap();
        assert!(matches!(
            PartialNttPlan::<19, 6, 3, 3, 2>::new(root),
            Err(PartialNttError::InvalidParameters)
        ));
        let rng = &mut test_rng();
        for _ in 0..16 {
            let a: [F; 6] = core::array::from_fn(|_| F::rand(rng));
            let b: [F; 6] = core::array::from_fn(|_| F::rand(rng));
            let left = Split::from_coefficient_array(a, root);
            let right = Split::from_coefficient_array(b, root);
            assert_eq!(left.coefficients(root), a);
            let values = left.ntt_array();
            for (j, factor) in [root, root * root].into_iter().enumerate() {
                for local in 0..3 {
                    assert_eq!(values[j * 3 + local], a[local] + factor * a[3 + local]);
                }
            }
            let mut expected = [F::ZERO; 11];
            for i in 0..6 {
                for j in 0..6 {
                    expected[i + j] += a[i] * b[j];
                }
            }
            // Reduce independently modulo Phi_3(X^3) = X^6 + X^3 + 1.
            for i in (6..11).rev() {
                let leading = expected[i];
                expected[i - 3] -= leading;
                expected[i - 6] -= leading;
                expected[i] = F::ZERO;
            }
            assert_eq!(
                left.mul_with_root(&right, root).coefficients(root),
                expected[..6]
            );
        }
    }

    #[test]
    fn plan_rejects_invalid_parameters_and_roots() {
        type Plan = PartialNttPlan<17, 16, 8, 4, 2>;
        assert!(matches!(
            Plan::new(Zq1::<17>::ZERO),
            Err(PartialNttError::InvalidRoot)
        ));
        assert!(matches!(
            Plan::new(Zq1::<17>::ONE),
            Err(PartialNttError::InvalidRoot)
        ));
        assert!(matches!(
            Plan::new(-Zq1::<17>::ONE),
            Err(PartialNttError::InvalidRoot)
        ));
        assert!(matches!(
            PartialNttPlan::<17, 16, 4, 4, 4>::new(Zq1::<17>::ONE),
            Err(PartialNttError::InvalidParameters)
        ));
        assert!(matches!(
            PartialNttPlan::<17, 16, 1, 32, 16>::new(Zq1::<17>::ONE),
            Err(PartialNttError::InvalidParameters)
        ));
        assert!(matches!(
            PartialNttPlan::<17, 16, 8, 0, 2>::new(Zq1::<17>::ONE),
            Err(PartialNttError::InvalidParameters)
        ));
        assert!(matches!(
            PartialNttPlan::<65, 4, 2, 4, 2>::new(Zq1::<65>::try_from(12u64).unwrap()),
            Err(PartialNttError::InvalidRoot)
        ));
    }
}
