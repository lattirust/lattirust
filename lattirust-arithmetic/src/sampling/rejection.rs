//! Rejection sampling for centered Gaussian shifts.
//!
//! Implements the standard and sign-leak LNP-Lite rejection probabilities.
//! The threshold uses exact integer comparisons and an adaptive Bernoulli
//! exponential sampler. It does not fix a projection dimension or a ring.
use ark_std::rand::{CryptoRng, RngCore};
use num_bigint::BigInt;
use num_traits::{One, Zero};
use std::cmp::Ordering;

/// x=log(M)+(2<z,v>-<v,v>)/(2 sigma²) has the form
/// (a*sqrt(t²)+b)/d, where t² is the public projection bound.
struct Exponent {
    a: BigInt,
    b: BigInt,
    d: BigInt,
    t2: BigInt,
}

impl Exponent {
    /// Compare x with n/d0 exactly. Callers supply d0>0 and a,t²≥0.
    fn cmp(&self, n: &BigInt, d0: &BigInt) -> Ordering {
        let right = n * &self.d - d0 * &self.b;
        if right < BigInt::zero() {
            return Ordering::Greater;
        }
        let left = d0 * &self.a;
        (&left * &left * &self.t2).cmp(&(&right * &right))
    }
    /// Bernoulli((x-base)/divisor), with the argument in [0,1]. Refine
    /// the uniform prefix whenever the algebraic threshold falls inside it.
    fn coin<G: RngCore>(&self, base: &BigInt, divisor: usize, rng: &mut G) -> bool {
        let mut lo = BigInt::zero();
        let mut den = BigInt::one();
        loop {
            lo = (lo << 64) + rng.next_u64();
            den <<= 64;
            let lower = base * &den + &lo * divisor;
            if self.cmp(&lower, &den) != Ordering::Greater {
                return false;
            }
            if self.cmp(&(lower + divisor), &den) != Ordering::Less {
                return true;
            }
        }
    }
    /// If a∈[0,1], successive continuation coins have probabilities a/j.
    /// Returning the parity of the first failed continuation gives exp(-a)
    /// by the alternating exponential series.
    fn exp_fraction<G: RngCore>(&self, base: &BigInt, rng: &mut G) -> bool {
        let mut j = 1usize;
        while self.coin(base, j, rng) {
            j += 1;
        }
        j % 2 == 1
    }
    fn accept<G: RngCore>(&self, rng: &mut G) -> bool {
        let zero = BigInt::zero();
        let one = BigInt::one();
        if self.cmp(&zero, &one) != Ordering::Greater {
            return true;
        }
        let unit = Self {
            a: zero.clone(),
            b: one.clone(),
            d: one.clone(),
            t2: zero.clone(),
        };
        let mut base = zero.clone();
        while self.cmp(&(&base + 1), &one) == Ordering::Greater {
            if !unit.exp_fraction(&zero, rng) {
                return false;
            }
            base += 1;
        }
        self.exp_fraction(&base, rng)
    }
}

#[derive(Clone, Copy)]
pub enum RejectionMode {
    Standard,
    SignLeak,
}

/// Accept the centered Gaussian response with the LNP-Lite rejection
/// probability. The public parameter bound applies to the full shift vector.
/// With sigma = (31/20)*2^log_scale and T^2 = shift_bound_squared,
/// Standard uses M = exp(14*T/sigma + T^2/(2*sigma^2)). SignLeak uses
/// M = exp(T^2/(2*sigma^2)) and rejects negative <response, shift>.
/// Acceptance probability is min(1, exp((-2<response,shift>+||shift||^2)
/// / (2*sigma^2)) / M). `None` means invalid shape, scale, bound, or overflow.
pub fn accept_masked_shift<G: CryptoRng + RngCore>(
    response: &[i64],
    shift: &[i64],
    log_scale: u32,
    shift_bound_squared: u128,
    mode: RejectionMode,
    rng: &mut G,
) -> Option<bool> {
    if log_scale > 32 || response.len() != shift.len() || shift.is_empty() {
        return None;
    }
    let vv = shift.iter().try_fold(0u128, |sum, &v| {
        let a = u128::from(v.unsigned_abs());
        sum.checked_add(a * a)
    })?;
    if vv > shift_bound_squared {
        return None;
    }
    let zv = response.iter().zip(shift).try_fold(0i128, |s, (&z, &v)| {
        s.checked_add(i128::from(z) * i128::from(v))
    })?;
    if matches!(mode, RejectionMode::SignLeak) && zv < 0 {
        return Some(false);
    }
    let sigma = 31u128 << log_scale;
    let t2 = BigInt::from(shift_bound_squared);
    let b = (t2.clone() + 2 * BigInt::from(zv) - BigInt::from(vv)) * 400;
    let a = match mode {
        RejectionMode::Standard => BigInt::from(560 * sigma),
        RejectionMode::SignLeak => BigInt::zero(),
    };
    Some(
        Exponent {
            a,
            b,
            d: BigInt::from(2 * sigma * sigma),
            t2,
        }
        .accept(rng),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_std::rand::{rngs::StdRng, SeedableRng};
    #[test]
    fn irrational_threshold_and_alternating_sampler() {
        let sqrt_two = Exponent {
            a: 1.into(),
            b: 0.into(),
            d: 1.into(),
            t2: 2.into(),
        };
        assert_eq!(
            sqrt_two.cmp(&1414213.into(), &1000000.into()),
            Ordering::Greater
        );
        assert_eq!(
            sqrt_two.cmp(&1414214.into(), &1000000.into()),
            Ordering::Less
        );
        let mut rng = StdRng::seed_from_u64(4242);
        let accepted = (0..4096).filter(|_| sqrt_two.accept(&mut rng)).count();
        assert!((880..1105).contains(&accepted));
        let negative = Exponent {
            a: 0.into(),
            b: (-1).into(),
            d: 1.into(),
            t2: 1.into(),
        };
        assert!(negative.accept(&mut rng));
    }

    #[test]
    fn rejects_invalid_shapes_bounds_and_overflow() {
        let mut rng = StdRng::seed_from_u64(1);
        assert_eq!(
            accept_masked_shift(&[], &[], 12, 1, RejectionMode::Standard, &mut rng),
            None
        );
        assert_eq!(
            accept_masked_shift(&[0], &[0], 33, 1, RejectionMode::Standard, &mut rng),
            None
        );
        assert_eq!(
            accept_masked_shift(&[0], &[2], 12, 3, RejectionMode::Standard, &mut rng),
            None
        );
        assert_eq!(
            accept_masked_shift(
                &[0],
                &[i64::MIN; 2],
                12,
                u128::MAX,
                RejectionMode::Standard,
                &mut rng
            ),
            None
        );
        assert_eq!(
            accept_masked_shift(
                &[i64::MAX; 3],
                &[i64::MAX; 3],
                12,
                u128::MAX,
                RejectionMode::Standard,
                &mut rng
            ),
            None
        );
        assert_eq!(
            accept_masked_shift(
                &[0; 4],
                &[i64::MIN; 4],
                12,
                u128::MAX,
                RejectionMode::Standard,
                &mut rng
            ),
            None
        );
    }

    #[test]
    fn sign_leak_mode_refuses_negative_inner_product() {
        let mut rng = StdRng::seed_from_u64(15);
        assert_eq!(
            accept_masked_shift(&[-1], &[1], 12, 1, RejectionMode::SignLeak, &mut rng),
            Some(false)
        );
        assert_eq!(
            accept_masked_shift(&[0, 0], &[0, 0], 12, 0, RejectionMode::Standard, &mut rng),
            Some(true)
        );
    }
}
