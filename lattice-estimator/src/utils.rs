use std::cmp::Ordering;

use statrs::function;

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


mod tests {
    use super::*;
    use std::cmp::Ordering;

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
        assert_eq!(trials, f64::MAX);
    }
}