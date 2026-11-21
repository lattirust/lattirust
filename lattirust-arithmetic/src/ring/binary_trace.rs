//! Binarity of packed coefficients in (Z/qZ)[X]/(X^N+1).
//!
//! Set d = nN for n polynomials. If the squared centered integer norm is at
//! most d, then the integer sum S = sum_i a_i(a_i - 1) is in [0, 2d].
//! For 2d < q, S = 0 mod q exactly when every coefficient is a bit.
//! The constant coefficient of sigma(a) * (a - (1 + ... + X^(N-1)))
//! gives this sum for one polynomial, where sigma maps X to X^(-1).
use super::{PolyRing, Pow2CyclotomicPolyRing, Ring, Zq1};
use crate::traits::WithL2Norm;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct BinaryShape<const Q: u64, const N: usize> {
    poly_count: usize,
}

impl<const Q: u64, const N: usize> BinaryShape<Q, N> {
    pub fn new(poly_count: usize) -> Option<Self> {
        let slots = poly_count.checked_mul(N)?;
        (poly_count != 0 && N.is_power_of_two() && 2 * (slots as u128) < u128::from(Q))
            .then_some(Self { poly_count })
    }

    pub fn polynomial_count(self) -> usize {
        self.poly_count
    }

    pub fn norm_bound_squared(self) -> u64 {
        (self.poly_count * N) as u64
    }

    pub fn slot_count(self) -> usize {
        self.poly_count * N
    }
}

/// Public bounded-quadratic relation. The norm and the trace must refer to
/// the same witness. Evaluating the relation privately emits no proof.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct BinaryTrace<const Q: u64, const N: usize> {
    pub shape: BinaryShape<Q, N>,
}

impl<const Q: u64, const N: usize> BinaryTrace<Q, N> {
    /// Compute the required constant coefficient in linear time. The centered
    /// norm is accumulated as an integer, independently of the ring sum.
    pub fn is_satisfied(self, polys: &[Pow2CyclotomicPolyRing<Zq1<Q>, N>]) -> bool {
        if polys.len() != self.shape.poly_count {
            return false;
        }
        let mut remaining = self.shape.norm_bound_squared();
        let mut residual = Zq1::<Q>::ZERO;
        for poly in polys {
            for v in poly.coefficients() {
                if v == Zq1::<Q>::ZERO {
                    continue;
                }
                let squared = if v == Zq1::<Q>::ONE {
                    1
                } else {
                    residual += v * (v - Zq1::<Q>::ONE);
                    let Ok(squared) = u64::try_from(v.l2_norm_squared()) else {
                        return false;
                    };
                    squared
                };
                let Some(next) = remaining.checked_sub(squared) else {
                    return false;
                };
                remaining = next;
            }
        }
        residual == Zq1::<Q>::ZERO
    }
}

/// Little-endian bits, in coefficient order. The last word has zero padding.
pub fn from_packed_bits<const Q: u64, const N: usize>(
    words: &[u64],
) -> Option<Pow2CyclotomicPolyRing<Zq1<Q>, N>> {
    if !N.is_power_of_two() || words.len() != N.div_ceil(64) {
        return None;
    }
    if N % 64 != 0 && (words[words.len() - 1] >> (N % 64)) != 0 {
        return None;
    }
    Some(Pow2CyclotomicPolyRing::from_fn(|i| {
        if words[i / 64] >> (i % 64) & 1 == 0 {
            Zq1::<Q>::ZERO
        } else {
            Zq1::<Q>::ONE
        }
    }))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::traits::WithConjugationAutomorphism;

    type F = Zq1<97>;
    type R = Pow2CyclotomicPolyRing<F, 4>;

    #[test]
    fn trace_matches_negacyclic_product_on_signed_coefficients() {
        let ones = R::from([F::ONE; 4]);
        for bits in 0_u64..256 {
            let poly = R::from_fn(|i| {
                let magnitude = F::try_from((bits >> (2 * i)) & 3).unwrap();
                if i % 2 == 0 {
                    magnitude
                } else {
                    -magnitude
                }
            });
            let reference = (poly.apply_automorphism() * (poly - ones)).coefficients()[0];
            let direct = poly
                .coefficients()
                .into_iter()
                .fold(F::ZERO, |r, v| r + v * (v - F::ONE));
            assert_eq!(reference, direct);
        }
    }

    #[test]
    fn norm_and_shape_are_enforced_without_modular_wrap() {
        let trace = BinaryTrace {
            shape: BinaryShape::<97, 4>::new(2).unwrap(),
        };
        assert!(trace.is_satisfied(&[
            from_packed_bits(&[10]).unwrap(),
            from_packed_bits(&[15]).unwrap()
        ]));
        assert!(!trace.is_satisfied(&[R::ONE]));
        assert!(!trace.is_satisfied(&[-R::ONE, R::ONE]));
        assert!(!trace.is_satisfied(&[R::ONE + R::ONE, R::ONE]));
        // 12*4 slots fit the strict modulus bound, 13*4 do not.
        assert!(BinaryShape::<97, 4>::new(12).is_some());
        assert!(BinaryShape::<97, 4>::new(13).is_none());
        assert!(BinaryShape::<97, 0>::new(2).is_none());
        assert!(from_packed_bits::<97, 4>(&[16]).is_none());
        assert!(from_packed_bits::<97, 4>(&[]).is_none());
    }

    #[test]
    fn packs_a_full_degree_256_witness() {
        let words = [u64::MAX, u64::MAX, u64::MAX, u64::MAX];
        let poly = from_packed_bits::<274_877_906_837, 256>(&words).unwrap();
        assert_eq!(
            poly.coefficients()
                .into_iter()
                .filter(|&v| v == Zq1::ONE)
                .count(),
            256
        );
        assert_eq!(poly.coefficients().len(), 256);

        let shape = BinaryShape::<274_877_906_837, 256>::new(1).unwrap();
        assert_eq!(shape.slot_count(), 256);
        assert!(BinaryTrace { shape }.is_satisfied(&[poly]));

        assert!(from_packed_bits::<274_877_906_837, 256>(&[0; 3]).is_none());
    }
}
