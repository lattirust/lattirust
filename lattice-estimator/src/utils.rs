use core::f64;
use std::cmp::Ordering;

const LOG_PRECISION: usize = 100;

#[derive(Debug, Default, Clone, Copy)]
pub struct CostParameters {
    pub rop: f64,
    pub red: f64,
    pub sieve: f64,
    pub block_size: usize,
    pub eta: Option<usize>,
    pub zeta: usize,
    pub dim: Option<usize>,
}

impl CostParameters {
    // Optional: Initialize with custom values using a builder pattern
    pub fn with_values(
        rop: f64,
        red: f64,
        sieve: f64,
        block_size: usize,
        eta: usize,
        zeta: usize,
        dim: usize,
    ) -> Self {
        CostParameters { rop, red, sieve, block_size, eta: Some(eta), zeta, dim: Some(dim) }
    }
}

pub fn binary_param_search<F>(function: F, min: usize, max: usize) -> (usize, CostParameters)
where
    F: Fn(usize) -> CostParameters,
{
    let mut left = min;
    let mut right = max;

    while left < right {
        let mid = left + (right - left) / 2;
        
        let f_mid: CostParameters = function(mid);
        let f_mid_next: CostParameters = function(mid + 1);

        // Compare the function values at mid and mid + 1
        match f_mid.rop.partial_cmp(&f_mid_next.rop) {
            Some(Ordering::Greater) => left = mid + 1, // Move left boundary right
            Some(Ordering::Less) => right = mid,       // Move right boundary left
            _ => return (mid, function(mid))                          // Stop if they are equal
        }
    }
    (left, function(left))
}

pub fn grid_param_search<F>(function: F, x_range: (usize, usize), y_range: (usize, usize)) -> (usize, usize, CostParameters)
where
    F: Fn(usize, usize) -> CostParameters
{
    let (x_min, x_max) = x_range;
    let (y_min, y_max) = y_range;

    let mut best_x = x_min;
    let mut best_y = y_min;

    let mut best_value = CostParameters::with_values(f64::INFINITY, 0.0, 0.0, 0, 0, 0, 0);

    for x in x_min..=x_max {
        for y in y_min..=y_max {
            let r = function(x, y);

            // Choose the parameter set (x, y) that minimizes `rop`
            if r.rop < best_value.rop {
                best_x = x;
                best_y = y;
                best_value = r;
            }
        }
    }

    (best_x, best_y, best_value)
}


pub fn adaptive_grid_search<F>(function: F, x_range: (usize, usize), y_range: (usize, usize), step_size: usize) -> (usize, usize, CostParameters)
where
    F: Fn(usize, usize) -> CostParameters
{
    let (x_min, x_max) = x_range;
    let (y_min, y_max) = y_range;

    let mut best_x = x_min;
    let mut best_y = y_min;

    let mut best_value: CostParameters = CostParameters::with_values(f64::INFINITY, 0.0, 0.0, 0, 0, 0, 0);

    //broad search
    for x in (x_min..=x_max).step_by(step_size) {
        for y in (y_min..=y_max).step_by(step_size) {
            let r: CostParameters = function(x, y);

            // Choose the parameter set (x, y) that minimizes `rop`
            if r.rop < best_value.rop {
                best_x = x;
                best_y = y;
                best_value = r;
            }
        }
    }

    //narrow search
    for x in best_x..=best_y + step_size {
        for y in best_y..=best_y + step_size {
            let r: CostParameters = function(x, y);

            // Choose the parameter set (x, y) that minimizes `rop`
            if r.rop < best_value.rop {
                best_x = x;
                best_y = y;
                best_value = r;
            }
        }
    }

    (best_x, best_y, best_value)
}



 /// The number of trials required to reach the target success probability.
pub fn amplify_via_trials(target_p: f64, success_p: f64) -> f64 {
    

    if target_p < success_p {
        return 1.0;
    } else if success_p == 0.0 {
        return f64::INFINITY;
    } else {
        // target_success_probability = 1 - (1 - success_probability)^trials
        return ((1.0 - target_p).ln() / (1.0 - success_p).ln()).ceil();
    }
}

//Taylor series of ln(1+x)
pub const fn taylor(y:f64, terms: usize) -> f64 {
    let mut result = 0.0;
        let mut term = y; // Start with the first term, which is just y
        for n in 1..=terms {
            if n % 2 == 1 {
                result += term / n as f64;
            } else {
                result -= term / n as f64;
            }
            term *= y; // Move to the next power of y
        }
    result
}

pub const fn ln_estimate(x: f64) -> f64 {
    if x <= 0.0 {
        panic!("ln is not defined for non-positive values.");
    }
        
    // Factor out powers of 2
    let mut n = 0;
    let mut reduced_x = x;
    while reduced_x > 2.0 {
        reduced_x /= 2.0;
        n += 1;
    }
    while reduced_x < 1.0 {
        reduced_x *= 2.0;
        n -= 1;
    }
    
    // Use the Taylor series on the "reduced" value
    let y = reduced_x - 1.0;
    let ln_approx = taylor(y, LOG_PRECISION);
    
    // Approximate ln(x) = n * ln(2) + ln(reduced_x)
    ln_approx + (n as f64) * f64::consts::LN_2
}

pub const fn log2_estimate(x: f64) -> f64 {
    ln_estimate(x) / f64::consts::LN_2
}


mod tests {
    use super::*;
    const EPSILON: f64 = 0.01;
    
    #[test]
    fn test_approximate_ln() {

        let test_cases = [
            (1.0, 0.0),              // ln(1) = 0
            (2.0, 0.693147),         // ln(2) ≈ 0.693147
            (10.0, 2.302585),        // ln(10) ≈ 2.302585
            (50.0, 3.912023),        // ln(50) ≈ 3.912023
            (0.5, -0.693147),        // ln(0.5) ≈ -0.693147
        ];

        for &(x, expected) in test_cases.iter() {
            let result = ln_estimate(x);
            assert!(
                (result - expected).abs() < EPSILON,
                "approximate_ln({}) ≈ {}, expected {}",
                x, result, expected
            );
        }
    }

    // Test cases for approximate_log2 function
    #[test]
    fn test_approximate_log2() {

        let test_cases = [
            (1.0, 0.0),             // log2(1) = 0
            (2.0, 1.0),             // log2(2) = 1
            (8.0, 3.0),             // log2(8) = 3
            (10.0, 3.321928),       // log2(10) ≈ 3.321928
            (50.0, 5.643856),       // log2(50) ≈ 5.643856
            (0.5, -1.0),            // log2(0.5) = -1
        ];

        for &(x, expected) in test_cases.iter() {
            let result = log2_estimate(x);
            assert!(
                (result - expected).abs() < EPSILON,
                "approximate_log2({}) ≈ {}, expected {}",
                x, result, expected
            );
        }
    }

    #[test]
    fn test_increasing_function() {
        // Test with a strictly increasing function, where the optimal value should be the minimum
        let (result, _) = binary_param_search(
            |x| CostParameters::with_values(x as f64, 0.0, 0.0, 0, 0, 0, 0),
            0,
            100,
        );
        assert_eq!(result, 0);
    }
    
    #[test]
    fn test_minimum_in_middle() {
        // Test with a function that has a minimum in the middle (like a parabola)
        let (result, _) = binary_param_search(
            |x| CostParameters::with_values(((x as isize - 50).pow(2)) as f64, 0.0, 0.0, 0, 0, 0, 0),
            0,
            100,
        );
        assert_eq!(result, 50);
    }

    #[test]
    fn test_plateau_function() {
        // Test with a flat (plateau) function where all values are the same
        let (result, _) = binary_param_search(
            |_| CostParameters::with_values(1.0, 0.0, 0.0, 0, 0, 0, 0),
            0,
            100,
        );
        // Should return the middle of the range for a plateau
        assert_eq!(result, 50);
    }

    #[test]
    fn test_amplify_basic_case() {
        // Test a basic case
        let trials = amplify_via_trials(0.95, 0.6);
        assert_eq!(trials, 4.0);
    }

    #[test]
    fn test_amplify_target_less_than_success() {
        // When target success probability is less than current success probability, expect 1 trial
        let trials = amplify_via_trials(0.5, 0.6);
        assert_eq!(trials, 1.0);
    }

    #[test]
    fn test_amplify_zero_success_probability() {
        // When success probability is 0, expect infinite trials (usize::MAX)
        let trials = amplify_via_trials(0.95, 0.0);
        assert_eq!(trials, f64::INFINITY);
    }
}