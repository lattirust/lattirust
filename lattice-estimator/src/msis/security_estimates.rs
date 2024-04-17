use log::debug;

use crate::errors::LatticeEstimatorError;
use crate::msis::MSIS;

/// Return the smallest `h` such that `MSIS\[h, w, d, q, length_bound\]` is $2^\lambda$-hard (for a given norm).
pub fn find_optimal_h(msis: &MSIS, lambda: usize) -> Result<usize, LatticeEstimatorError> {
    let lo: usize = 1;
    let hi: usize = msis.upper_bound_h();

    debug!("Searching optimal h in [{lo}, {hi}] for {msis} with target lambda={lambda}...");
    // Start from h=1 and exhaustively test upwards, since the upper bound hi is often not very large, and computing the security level for larger h is expensive.
    for h in lo..=hi {
        let curr = msis.with_h(h);
        let lambda_curr = curr.security_level();
        debug!("\t{curr} -> {lambda_curr} bits of security");
        if lambda_curr >= lambda as f64 {
            return Ok(h);
        }
    }
    Err(LatticeEstimatorError::from(format!(
        "no suitable n found in range [{lo}, {hi}] for {msis}"
    )))
}

/// Return the smallest h such that `MSIS\[h, w, d, q, length_bound(h)\]` is $2^\lambda$-hard (for a given norm), where `length_bound` is a function of `h`.
pub fn find_optimal_n_dynamic<F>(
    msis: &MSIS,
    length_bound: F,
    lambda: usize,
) -> Result<usize, LatticeEstimatorError>
where
    F: Fn(usize) -> f64,
{
    let lo: usize = 1;
    let hi: usize = msis.upper_bound_h();

    // Start from h=1 and exhaustively test upwards, since the upper bound hi is often not very large, and computing the security level for larger h is expensive.
    debug!("Searching optimal h in [{lo}, {hi}] for {msis} with target lambda={lambda}...");
    for h in lo..=hi {
        let curr = msis.with_h(h).with_length_bound(length_bound(h));
        let lambda_curr = curr.security_level();
        debug!("\t{curr} -> {lambda_curr} bits of security");
        if lambda_curr >= lambda as f64 {
            return Ok(h);
        }
    }
    Err(LatticeEstimatorError::from(format!(
        "no suitable h found in range [{lo}, {hi}] for {msis}"
    )))
}

#[cfg(test)]
mod test {
    use crate::norms::Norm;

    use super::*;

    const Q: u128 = 2147483649;
    const SQRT_Q: f64 = 46340.95;
    const TEST_L2: MSIS = MSIS {
        h: 5,
        d: 64,
        q: Q,
        length_bound: SQRT_Q,
        w: 512,
        norm: Norm::L2,
    };
    const TEST_LINF: MSIS = MSIS {
        h: 2,
        d: 64,
        q: Q,
        length_bound: 1.,
        w: 512,
        norm: Norm::Linf,
    };

    #[test]
    fn test_msis_security_level_l2() {
        let lambda = TEST_L2.security_level();
        println!("{TEST_L2} -> lambda: {lambda}");
    }

    #[test]
    fn test_msis_security_level_linf() {
        let lambda = TEST_LINF.security_level();
        println!("{TEST_LINF} -> lambda: {lambda}");
    }

    #[test]
    fn test_find_optimal_n_l2() {
        let h_opt = find_optimal_h(&TEST_L2.with_h(0), 128).unwrap();
        let msis = TEST_L2.with_h(h_opt);
        println!("{TEST_L2} -> lambda: {}", TEST_L2.security_level());
        println!("{msis} -> lambda: {}", msis.security_level());
    }

    #[test]
    fn test_find_optimal_n_linf() {
        let h_opt = find_optimal_h(&TEST_LINF.with_h(0), 128).unwrap();
        let msis = TEST_LINF.with_h(h_opt);
        println!("{TEST_L2} -> lambda: {}", TEST_LINF.security_level());
        println!("{msis} -> lambda: {}", msis.security_level());
    }
}
