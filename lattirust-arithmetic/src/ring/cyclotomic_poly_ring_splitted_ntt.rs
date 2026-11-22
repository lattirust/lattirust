use std::fmt::{Display, Formatter};
use std::ops::{Add, AddAssign, Neg, Sub, SubAssign};

use crate::linear_algebra::SVector;
use crate::partial_ntt::{coprimes_set, PartialNTT, PartialNttPlan};
use crate::ring::{Ring, Zq1};

/// Transform-domain representation for rings that split only into factors
/// `x^D - r_j` instead of linear factors.
///
/// For power-of-two cyclotomics this corresponds to the factorization
/// `x^N + 1 = \prod_j (x^D - r_j)` where the `r_j` are the primitive
/// `Z`-th roots of unity and `N = D * phi(Z)`.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
pub struct CyclotomicPolyRingSplittedNTT<
    const Q: u64,
    const N: usize,
    const D: usize,
    const Z: usize,
    const PHI_Z: usize,
>(SVector<Zq1<Q>, N>);

impl<const Q: u64, const N: usize, const D: usize, const Z: usize, const PHI_Z: usize> Display
    for CyclotomicPolyRingSplittedNTT<Q, N, D, Z, PHI_Z>
{
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        Display::fmt(&self.0, f)
    }
}

impl<const Q: u64, const N: usize, const D: usize, const Z: usize, const PHI_Z: usize>
    CyclotomicPolyRingSplittedNTT<Q, N, D, Z, PHI_Z>
{
    pub fn primitive_zth_root() -> Zq1<Q> {
        assert!(Z.is_power_of_two(), "Z must be a power of two");
        assert_eq!(Q % Z as u64, 1, "Q must be congruent to 1 modulo Z");

        let two_adicity = (Q - 1).trailing_zeros();
        let z_adicity = Z.trailing_zeros();
        let odd_part = (Q - 1) >> two_adicity;

        for candidate in 2..Q {
            let candidate = Zq1::<Q>::try_from(candidate).unwrap();
            let two_adic_root = candidate.pow(odd_part);

            if two_adic_root.pow(1u64 << (two_adicity - 1)) != Zq1::<Q>::ONE {
                return two_adic_root.pow(1u64 << (two_adicity - z_adicity));
            }
        }

        panic!("could not find a primitive Z-th root of unity");
    }

    pub fn from_coefficient_array(mut coeffs: [Zq1<Q>; N], primitive_zth_root: Zq1<Q>) -> Self {
        Self::ntt(&mut coeffs, primitive_zth_root);
        Self::from_ntt_array(coeffs)
    }

    /// Construct using reusable public transform parameters.
    pub fn from_coefficient_array_with_plan(
        mut coeffs: [Zq1<Q>; N],
        plan: &PartialNttPlan<Q, N, D, Z, PHI_Z>,
    ) -> Self {
        plan.ntt(&mut coeffs);
        Self::from_ntt_array(coeffs)
    }

    pub fn from_ntt_array(coeffs_ntt: [Zq1<Q>; N]) -> Self {
        Self(SVector::const_from_array(coeffs_ntt))
    }

    pub fn from_fn<F>(f: F, primitive_zth_root: Zq1<Q>) -> Self
    where
        F: FnMut(usize) -> Zq1<Q>,
    {
        Self::from_coefficient_array(core::array::from_fn(f), primitive_zth_root)
    }

    pub fn coefficients(&self, primitive_zth_root: Zq1<Q>) -> [Zq1<Q>; N] {
        let mut coeffs = self.ntt_array();
        Self::intt(&mut coeffs, primitive_zth_root);
        coeffs
    }

    /// Invert using the same root used to transform this element.
    pub fn coefficients_with_plan(&self, plan: &PartialNttPlan<Q, N, D, Z, PHI_Z>) -> [Zq1<Q>; N] {
        let mut coeffs = self.ntt_array();
        plan.intt(&mut coeffs);
        coeffs
    }

    pub fn ntt_array(&self) -> [Zq1<Q>; N] {
        self.0 .0.data.0[0]
    }

    pub fn ntt_values(&self) -> Vec<Zq1<Q>> {
        self.0.iter().cloned().collect()
    }

    pub fn mul_with_root(&self, rhs: &Self, primitive_zth_root: Zq1<Q>) -> Self {
        let roots = coprimes_set::<Z, PHI_Z>().map(|power| primitive_zth_root.pow(power as u64));
        self.mul_with_factors(rhs, &roots)
    }

    /// Multiply using precomputed factor roots. The plan must use the root
    /// used to transform both operands.
    pub fn mul_with_plan(&self, rhs: &Self, plan: &PartialNttPlan<Q, N, D, Z, PHI_Z>) -> Self {
        self.mul_with_factors(rhs, plan.roots())
    }

    fn mul_with_factors(&self, rhs: &Self, roots: &[Zq1<Q>; PHI_Z]) -> Self {
        debug_assert_eq!(N, D * PHI_Z);
        let lhs = self.ntt_array();
        let rhs = rhs.ntt_array();
        let mut products = [Zq1::<Q>::ZERO; N];
        for (component_idx, &root) in roots.iter().enumerate() {
            let base = component_idx * D;

            for lhs_idx in 0..D {
                for rhs_idx in 0..D {
                    let mut coeff = lhs[base + lhs_idx] * rhs[base + rhs_idx];
                    let target = lhs_idx + rhs_idx;
                    let out_idx = base + (target % D);
                    if target >= D {
                        coeff *= root;
                    }
                    products[out_idx] += coeff;
                }
            }
        }

        Self::from_ntt_array(products)
    }
}

impl<const Q: u64, const N: usize, const D: usize, const Z: usize, const PHI_Z: usize>
    PartialNTT<Q, N, D, Z, PHI_Z> for CyclotomicPolyRingSplittedNTT<Q, N, D, Z, PHI_Z>
{
}

impl<const Q: u64, const N: usize, const D: usize, const Z: usize, const PHI_Z: usize> Add
    for CyclotomicPolyRingSplittedNTT<Q, N, D, Z, PHI_Z>
{
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self::from_ntt_array(core::array::from_fn(|idx| self.0[idx] + rhs.0[idx]))
    }
}

impl<const Q: u64, const N: usize, const D: usize, const Z: usize, const PHI_Z: usize> AddAssign
    for CyclotomicPolyRingSplittedNTT<Q, N, D, Z, PHI_Z>
{
    fn add_assign(&mut self, rhs: Self) {
        for idx in 0..N {
            self.0[idx] += rhs.0[idx];
        }
    }
}

impl<const Q: u64, const N: usize, const D: usize, const Z: usize, const PHI_Z: usize> Sub
    for CyclotomicPolyRingSplittedNTT<Q, N, D, Z, PHI_Z>
{
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self::from_ntt_array(core::array::from_fn(|idx| self.0[idx] - rhs.0[idx]))
    }
}

impl<const Q: u64, const N: usize, const D: usize, const Z: usize, const PHI_Z: usize> SubAssign
    for CyclotomicPolyRingSplittedNTT<Q, N, D, Z, PHI_Z>
{
    fn sub_assign(&mut self, rhs: Self) {
        for idx in 0..N {
            self.0[idx] -= rhs.0[idx];
        }
    }
}

impl<const Q: u64, const N: usize, const D: usize, const Z: usize, const PHI_Z: usize> Neg
    for CyclotomicPolyRingSplittedNTT<Q, N, D, Z, PHI_Z>
{
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self::from_ntt_array(core::array::from_fn(|idx| -self.0[idx]))
    }
}

#[cfg(test)]
mod tests {
    use ark_std::{test_rng, UniformRand};

    use super::CyclotomicPolyRingSplittedNTT;
    use crate::ring::{Pow2CyclotomicPolyRing, Zq1};

    const Q: u64 = 17;
    const N: usize = 16;
    const Z: usize = 4;
    const D: usize = 8;
    const PHI_Z: usize = 2;

    type BR = Zq1<Q>;
    type SplitNTT = CyclotomicPolyRingSplittedNTT<Q, N, D, Z, PHI_Z>;
    type CoeffRing = Pow2CyclotomicPolyRing<BR, N>;

    #[test]
    fn roundtrip_matches_coefficients() {
        let root = SplitNTT::primitive_zth_root();
        let rng = &mut test_rng();

        for _ in 0..32 {
            let coeff_ring = CoeffRing::rand(rng);
            let coeffs: [BR; N] = coeff_ring.coefficient_array();
            let split = SplitNTT::from_coefficient_array(coeffs, root);
            assert_eq!(split.coefficients(root), coeffs);
        }
    }

    #[test]
    fn multiplication_matches_pow2_cyclotomic_ring() {
        let root = SplitNTT::primitive_zth_root();
        let rng = &mut test_rng();

        for _ in 0..32 {
            let lhs = CoeffRing::rand(rng);
            let rhs = CoeffRing::rand(rng);

            let lhs_split = SplitNTT::from_coefficient_array(lhs.coefficient_array(), root);
            let rhs_split = SplitNTT::from_coefficient_array(rhs.coefficient_array(), root);

            let product_split = lhs_split.mul_with_root(&rhs_split, root);
            let product_coeffs = product_split.coefficients(root);
            let expected = lhs * rhs;

            assert_eq!(product_coeffs, expected.coefficient_array());
        }
    }
}
