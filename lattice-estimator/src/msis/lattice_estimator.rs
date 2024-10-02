use crate::errors::LatticeEstimatorError;
use crate::msis::MSIS;

/// Return the smallest `h` such that `MSIS\[h, w, d, q, length_bound\]` is $2^\lambda$-hard (for a given norm).
pub fn find_optimal_h(msis: &MSIS, lambda: usize) -> Result<usize, LatticeEstimatorError> {
    let mut lo_sis: usize = 1;
    let mut hi_sis: usize = msis.to_sis().upper_bound_h();

    // Make sure we have a valid initial interval
    let msis = msis.with_h(hi_sis.div_floor(msis.d));

    let lambda_hi = msis.security_level();
    debug_assert!(
        lambda_hi >= lambda as f64,
        "{msis} has sec. param. {lambda_hi}  < target lambda = {lambda}"
    );
    // Loop invariant: SIS_{hi_sis * d, q, m*d, length_bound} is 2^lambda_hi-hard with lambda_hi >= lambda
    while hi_sis > lo_sis {
        let mid_sis = lo_sis + (hi_sis - lo_sis) / 2;
        // Use closest multiple of d for the MSIS instance
        let msis = msis.with_h((mid_sis as f64 / msis.d as f64).round() as usize);
        let mid_lambda = msis.security_level();
        if mid_lambda >= lambda as f64 {
            // Search for smaller n in [lo, mid]
            hi_sis = mid_sis;
        } else {
            // Search for smaller n in [mid+1, hi]
            lo_sis = mid_sis + 1;
        }
    }
    assert_eq!(hi_sis, lo_sis);
    Ok(hi_sis)
}

/// Return the smallest h such that `MSIS\[h, w, d, q, length_bound(h)\]` is $2^\lambda$-hard (for a given norm), where `length_bound` is a function of `h`.
pub fn find_optimal_h_dynamic<F>(
    msis: &MSIS,
    length_bound: F,
    lambda: usize,
) -> Result<usize, LatticeEstimatorError>
where
    F: Fn(usize) -> f64,
{
    // We can't assume that length_bound is monotonic, so we can't use binary search.
    // Instead, exhaustively search powers of 2 until we find a suitable n.
    // TODO: use a better search algorithm / return a more fine-grained result
    let hi = msis.d * msis.to_sis().upper_bound_h();
    let candidates = (1..hi).map(|i| 2usize.pow(i as u32));
    candidates
        .map(|h| msis.with_h(h).with_length_bound(length_bound(h)))
        .find(|sis| sis.security_level() >= lambda as f64)
        .map(|sis| sis.h)
        .ok_or(LatticeEstimatorError::from(
            "no suitable h found".to_string(),
        ))
}

#[cfg(test)]
mod test {
    use crate::{norms::Norm, reduction::Estimates};

    use super::*;

    const Q: u128 = 2147483649;
    const SQRT_Q: f64 = 46340.95;

    #[test]
    fn test_msis_security_level_l2() {
        let test_l2: MSIS = MSIS {
            h: 5,
            d: 64,
            q: Q.into(),
            length_bound: SQRT_Q,
            w: 512,
            norm: Norm::L2,
        };
        let lambda: f64 = test_l2.security_level();
        println!("External {test_l2} -> lambda: {lambda}");

        let lambda: f64 = test_l2.security_level_internal(Estimates::LLL);
        println!("Internal  {test_l2} -> lambda: {lambda}");
    }

    #[test]
    fn test_msis_security_level_linf() {
        let test_linf: MSIS = MSIS {
            h: 2,
            d: 64,
            q: Q.into(),
            length_bound: 1.,
            w: 512,
            norm: Norm::Linf,
        };
        let lambda = test_linf.security_level();
        println!("{test_linf} -> lambda: {lambda}");
    }

    #[test]
    fn test_find_optimal_h_l2() {
        let test_l2: MSIS = MSIS {
            h: 5,
            d: 64,
            q: Q.into(),
            length_bound: SQRT_Q,
            w: 512,
            norm: Norm::L2,
        };
        let h_opt =
            crate::msis::security_estimates::find_optimal_h(&test_l2.with_h(0), 128).unwrap();
        let msis = test_l2.with_h(h_opt);
        println!("{test_l2} -> lambda: {}", test_l2.security_level());
        println!("{msis} -> lambda: {}", msis.security_level());
    }

    #[test]
    fn test_find_optimal_h_linf() {
        let test_linf: MSIS = MSIS {
            h: 2,
            d: 64,
            q: Q.into(),
            length_bound: 1.,
            w: 512,
            norm: Norm::Linf,
        };
        let h_opt =
            crate::msis::security_estimates::find_optimal_h(&test_linf.with_h(0), 128).unwrap();
        let msis = test_linf.with_h(h_opt);
        println!("{test_linf} -> lambda: {}", test_linf.security_level());
        println!("{msis} -> lambda: {}", msis.security_level());
    }
}
