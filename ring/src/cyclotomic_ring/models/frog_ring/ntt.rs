//!
//! A CRT implementation for the ring Fq[X]/(X^16+1).
//!
use ark_ff::BigInt;
use ark_std::vec::*;

use super::{Fq, Fq4};

/// The degree of the cyclotomic polynomial.
pub(super) const D: usize = 16;

/// The number of splits of the cyclotomic ring in the CRT-form.
pub(super) const N: usize = 4;

/// All 8 roots of unity of degree 8.
const ROOTS_OF_UNITY_8: &[Fq] = &[
    Fq::new(BigInt([1u64])),
    Fq::new(BigInt([2755067726615789629u64])),
    Fq::new(BigInt([13238044465818905414u64])),
    Fq::new(BigInt([8043592722274778300u64])),
    Fq::new(BigInt([15912092521325583640u64])), // -1
    Fq::new(BigInt([13157024794709794012u64])),
    Fq::new(BigInt([2674048055506678227u64])),
    Fq::new(BigInt([7868499799050805341u64])),
];

/// 1 / 4
const FOUR_INV: Fq = Fq::new(BigInt([11934069390994187731u64]));

/// Given `coefficients` of a polynoimial `f mod X^16 + 1`
/// returns its CRT:
/// *  `f mod X^3-NONRESIDUE`,
/// *  `f mod X^3-NONRESIDUE^5`,
/// *  `f mod X^3-NONRESIDUE^3`,
/// *  `f mod X^3-NONRESIDUE^7`.
///
/// Each of the components is transformed into an element of `Fq4`
/// by the corresponding unique isomorphism `Fq[X]/(X^4-NONRESIDUE^i) -> Fq[X]/(X^4-NONRESIDUE)`.
///
/// # Panics
///
/// Panics if `coefficients.len() != 16`.
#[inline(always)]
pub fn frog_crt(coefficients: Vec<Fq>) -> Vec<Fq4> {
    // TODO: once we have parallelization set up
    //       this should be either parallel or serial.
    serial_frog_crt(coefficients)
}

/// Same as `frog_crt` but performs the CRT in place.
/// Takes coefficients of a polynomial and outputs
/// components of the corresponding CRT factors in each
/// quadriple `coefficients[4*i]`, `coefficients[4*i + 1]`, `coefficients[4*i+2]`, `coefficients[4*i+3]`.
///
/// # Panics
///
/// Panics if `coefficients.len() != 16`.
#[inline(always)]
pub fn frog_crt_in_place(coefficients: &mut [Fq]) {
    // TODO: once we have parallelization set up
    //       this should be either parallel or serial.
    serial_frog_crt_in_place(coefficients)
}

/// The inverse CRT.
/// Takes the CRT representation in the order:
/// *  `f mod X^3-NONRESIDUE`,
/// *  `f mod X^3-NONRESIDUE^5`,
/// *  `f mod X^3-NONRESIDUE^3`,
/// *  `f mod X^3-NONRESIDUE^7`.
///
/// Each of the components is in its isomorphic form in the `Fq[X]/(X^4-NONRESIDUE)`.
/// Returns the coefficients of the polynomial encoded by this CRT form.
///
/// # Panics
///
/// Panics if `evaluations.len() != 16`.
#[inline(always)]
pub fn frog_icrt(evaluations: Vec<Fq4>) -> Vec<Fq> {
    // TODO: once we have parallelization set up
    //       this should be either parallel or serial.
    serial_frog_icrt(evaluations)
}

/// Same as `frog_icrt` but performs the inverse CRT in place.
/// Each triple `evaluations[4*i]`, `evaluations[4*i + 1]`, `evaluations[4*i+2]`, `evaluations[4*i+3]`
/// has to be an `Fq4` element. In the order described in `frog_icrt`'s docstring.
///
/// # Panics
///
/// Panics if `evaluations.len() != 16`.
#[inline(always)]
pub fn frog_icrt_in_place(evaluations: &mut [Fq]) {
    // TODO: once we have parallelization set up
    //       this should be either parallel or serial.
    serial_frog_icrt_in_place(evaluations)
}

fn serial_frog_crt(mut coefficients: Vec<Fq>) -> Vec<Fq4> {
    serial_frog_crt_in_place(&mut coefficients);

    super::utils::fq_vec_to_fq4_vec(coefficients)
}

fn serial_frog_crt_in_place(coefficients: &mut [Fq]) {
    assert!(coefficients.len() == D);

    // Compute f mod X^8-NONRESIDUE^2, f mod X^8+NONRESIDUE^2
    for i in 0..(D / 2) {
        // f_1
        {
            let (coeff_i, coeff_d_div_4_plus_i) = (coefficients[i], coefficients[D / 2 + i]);
            let sigma_coeff_d_div_2_plus_i = ROOTS_OF_UNITY_8[2] * coeff_d_div_4_plus_i;

            coefficients[i] = coeff_i + sigma_coeff_d_div_2_plus_i;
            coefficients[i + D / 2] = coeff_i - sigma_coeff_d_div_2_plus_i;
        }
    }

    for i in 0..(D / 4) {
        // f_1
        {
            let (coeff_i, coeff_d_div_4_plus_i) = (coefficients[i], coefficients[D / 4 + i]);
            let sigma_coeff_d_div_2_plus_i = ROOTS_OF_UNITY_8[1] * coeff_d_div_4_plus_i;

            coefficients[i] = coeff_i + sigma_coeff_d_div_2_plus_i;
            coefficients[i + D / 4] = coeff_i - sigma_coeff_d_div_2_plus_i;
        }

        // f_2
        {
            let (coeff_i, coeff_d_div_4_plus_i) =
                (coefficients[D / 2 + i], coefficients[3 * D / 4 + i]);
            let sigma_coeff_d_div_4_plus_i = ROOTS_OF_UNITY_8[3] * coeff_d_div_4_plus_i;

            coefficients[D / 2 + i] = coeff_i + sigma_coeff_d_div_4_plus_i;
            coefficients[3 * D / 4 + i] = coeff_i - sigma_coeff_d_div_4_plus_i;
        }
    }

    homogenize_fq4(coefficients);
}

fn serial_frog_icrt(evaluations: Vec<Fq4>) -> Vec<Fq> {
    assert_eq!(evaluations.len(), N);

    let mut evaluations = super::utils::fq4_vec_to_fq_vec(evaluations);

    serial_frog_icrt_in_place(&mut evaluations);

    evaluations
}

fn serial_frog_icrt_in_place(evaluations: &mut [Fq]) {
    assert_eq!(evaluations.len(), D);

    dehomogenize_fq4(evaluations);

    for i in 0..(D / 4) {
        // f_1
        {
            let (coeff_i, coeff_d_div_4_plus_i) = (evaluations[i], evaluations[D / 4 + i]);

            evaluations[i] = coeff_i + coeff_d_div_4_plus_i;
            evaluations[i + D / 4] = ROOTS_OF_UNITY_8[7] * (coeff_i - coeff_d_div_4_plus_i);
        }

        // f_2
        {
            let (coeff_i, coeff_d_div_4_plus_i) =
                (evaluations[D / 2 + i], evaluations[3 * D / 4 + i]);

            evaluations[D / 2 + i] = coeff_i + coeff_d_div_4_plus_i;
            evaluations[3 * D / 4 + i] = ROOTS_OF_UNITY_8[5] * (coeff_i - coeff_d_div_4_plus_i);
        }
    }

    for i in 0..(D / 2) {
        // f_1
        {
            let (coeff_i, coeff_d_div_2_plus_i) = (evaluations[i], evaluations[D / 2 + i]);

            evaluations[i] = FOUR_INV * (coeff_i + coeff_d_div_2_plus_i);
            evaluations[i + D / 2] =
                FOUR_INV * (ROOTS_OF_UNITY_8[6] * (coeff_i - coeff_d_div_2_plus_i));
        }
    }
}

/// At the end of CRT we get elements in different (although isomorphic)
/// degree-4 extensions of Fq,
/// this function converts each triple `c[4 * i]`, `c[4 * i + 1]`, `c[4 * i + 2]`, `c[4 * i + 3]`
/// into coefficients of an element from Fq[X]/(X^4-NONRESIDUE)=Fq4.
#[inline(always)]
fn homogenize_fq4(c: &mut [Fq]) {
    nonresidue_to_nonresidue(&mut c[0..4]);
    nonresidue_to_5_to_nonresidue(&mut c[4..8]);
    nonresidue_to_3_to_nonresidue(&mut c[8..12]);
    nonresidue_to_7_to_nonresidue(&mut c[12..16]);
}

/// The inverse of the above.
#[inline(always)]
fn dehomogenize_fq4(c: &mut [Fq]) {
    nonresidue_to_nonresidue(&mut c[0..4]);
    nonresidue_to_nonresidue_to_5(&mut c[4..8]);
    nonresidue_to_nonresidue_to_3(&mut c[8..12]);
    nonresidue_to_nonresidue_to_7(&mut c[12..16]);
}

#[inline(always)]
fn nonresidue_to_nonresidue(c: &mut [Fq]) {
    c.swap(2, 1);
}

#[inline(always)]
fn nonresidue_to_nonresidue_to_5(c: &mut [Fq]) {
    let c2 = c[2];
    c[2] = ROOTS_OF_UNITY_8[6] * c[1];
    c[1] = ROOTS_OF_UNITY_8[7] * c2;
    c[3] *= ROOTS_OF_UNITY_8[5];
}

// Different automorphisms with the target Fp(NONRESIDUE) and their inverses.
#[inline(always)]
fn nonresidue_to_5_to_nonresidue(c: &mut [Fq]) {
    let c2 = c[2];
    c[2] = ROOTS_OF_UNITY_8[1] * c[1];
    c[1] = ROOTS_OF_UNITY_8[2] * c2;
    c[3] *= ROOTS_OF_UNITY_8[3];
}

#[inline(always)]
fn nonresidue_to_3_to_nonresidue(c: &mut [Fq]) {
    let c3 = c[3];
    c[3] = -c[1];
    c[1] = ROOTS_OF_UNITY_8[1] * c[2];
    c[2] = ROOTS_OF_UNITY_8[6] * c3;
}

#[inline(always)]
fn nonresidue_to_nonresidue_to_3(c: &mut [Fq]) {
    let c3 = c[3];
    c[3] = ROOTS_OF_UNITY_8[2] * c[2];
    c[2] = ROOTS_OF_UNITY_8[7] * c[1];
    c[1] = -c3;
}

#[inline(always)]
fn nonresidue_to_7_to_nonresidue(c: &mut [Fq]) {
    let c3 = c[3];
    c[3] = ROOTS_OF_UNITY_8[1] * c[1];
    c[1] = ROOTS_OF_UNITY_8[3] * c[2];
    c[2] = ROOTS_OF_UNITY_8[5] * c3;
}

#[inline(always)]
fn nonresidue_to_nonresidue_to_7(c: &mut [Fq]) {
    let c3 = c[3];
    c[3] = ROOTS_OF_UNITY_8[3] * c[2];
    c[2] = ROOTS_OF_UNITY_8[5] * c[1];
    c[1] = ROOTS_OF_UNITY_8[7] * c3;
}

#[cfg(test)]
mod tests {
    use ark_ff::{Field, MontFp, UniformRand};
    use ark_std::Zero;
    use rand::SeedableRng;
    use rand_chacha::ChaCha8Rng;

    use super::*;
    use crate::cyclotomic_ring::models::frog_ring::{Fq2, Fq4};

    #[test]
    fn test_roots_of_unity() {
        // Check if they all are roots of unity of degree 24.
        for root in ROOTS_OF_UNITY_8 {
            assert_eq!(root.pow([8]), Fq::ONE);
        }

        // Check if they are pairwise distinct.
        for (i, elem_i) in ROOTS_OF_UNITY_8.iter().enumerate().take(8) {
            for elem_j in ROOTS_OF_UNITY_8.iter().take(8).skip(i + 1) {
                assert_ne!(elem_i, elem_j);
            }
        }

        // Check if they come in order.
        for i in 0..8u64 {
            assert_eq!(ROOTS_OF_UNITY_8[i as usize], ROOTS_OF_UNITY_8[1].pow([i]));
        }
    }

    macro_rules! test_inverses {
        ($($isomorphism:expr, $inverse:expr),+) => {
            let mut rng = ChaCha8Rng::seed_from_u64(0);
            $({
                let x = Fq4::rand(&mut rng);
                let mut x_prime: Vec<Fq> = vec![x.c0.c0, x.c0.c1, x.c1.c0, x.c1.c1];
                $isomorphism(&mut x_prime);
                $inverse(&mut x_prime);
                assert_eq!(Fq4::new(Fq2::new(x_prime[0], x_prime[1]), Fq2::new(x_prime[2], x_prime[3])), x);
            })+
        };
    }

    #[test]
    fn test_isomorphisms_inverses() {
        test_inverses! {
            nonresidue_to_nonresidue, nonresidue_to_nonresidue,
            nonresidue_to_nonresidue_to_7, nonresidue_to_7_to_nonresidue,
            nonresidue_to_nonresidue_to_5, nonresidue_to_5_to_nonresidue,
            nonresidue_to_nonresidue_to_3, nonresidue_to_3_to_nonresidue
        };
    }

    macro_rules! test_powers_of_x {
        ($($inverse:expr),+) => {
            $({
                let mut x_prime: Vec<Fq> = vec![Fq::ZERO, Fq::ONE, Fq::ZERO, Fq::ZERO];
                let mut x_prime_2: Vec<Fq> = vec![Fq::ZERO, Fq::ZERO, Fq::ONE, Fq::ZERO];
                let mut x_prime_3: Vec<Fq> = vec![Fq::ZERO, Fq::ZERO, Fq::ZERO, Fq::ONE];
                $inverse(&mut x_prime);
                $inverse(&mut x_prime_2);
                $inverse(&mut x_prime_3);

                let x_prime = Fq4::new(Fq2::new(x_prime[0], x_prime[1]), Fq2::new(x_prime[2], x_prime[3]));
                let x_prime_2 = Fq4::new(Fq2::new(x_prime_2[0], x_prime_2[1]), Fq2::new(x_prime_2[2], x_prime_2[3]));
                let x_prime_3 = Fq4::new(Fq2::new(x_prime_3[0], x_prime_3[1]), Fq2::new(x_prime_3[2], x_prime_3[3]));

                assert_eq!(x_prime * x_prime, x_prime_2);
                assert_eq!(x_prime * x_prime * x_prime, x_prime_3);
            })+
        };
    }

    #[test]
    fn test_powers() {
        test_powers_of_x! {
            nonresidue_to_nonresidue,
            nonresidue_to_7_to_nonresidue,
            nonresidue_to_5_to_nonresidue,
            nonresidue_to_3_to_nonresidue
        };
    }

    macro_rules! test_extension_equations {
        ($($root_of_unity:expr, $inverse:expr),+) => {
            $({
                let mut x_prime: Vec<Fq> = vec![Fq::ZERO, Fq::ONE, Fq::ZERO, Fq::ZERO];
                $inverse(&mut x_prime);

                let x_prime = Fq4::new(Fq2::new(x_prime[0], x_prime[1]), Fq2::new(x_prime[2], x_prime[3]));

                assert_eq!(x_prime * x_prime * x_prime * x_prime - Fq4::from_base_prime_field($root_of_unity), Fq4::ZERO);
            })+
        };
    }

    #[test]
    fn test_equations() {
        test_extension_equations! {
            ROOTS_OF_UNITY_8[7], nonresidue_to_7_to_nonresidue,
            ROOTS_OF_UNITY_8[5], nonresidue_to_5_to_nonresidue,
            ROOTS_OF_UNITY_8[3], nonresidue_to_3_to_nonresidue
        };
    }

    #[test]
    fn test_normalize_denormalize() {
        let mut rng = ChaCha8Rng::seed_from_u64(0);
        let mut x: Vec<Fq> = (0..16).map(|_| Fq::rand(&mut rng)).collect();

        let expected = x.clone();

        homogenize_fq4(&mut x);
        dehomogenize_fq4(&mut x);

        assert_eq!(x, expected);
    }

    #[test]
    fn test_crt() {
        // 1 + 2 * X + 3 * X^2 + 15 * X^15 + X^23
        let mut test_poly: Vec<Fq> = vec![
            Fq::from(1),
            Fq::from(2),
            Fq::from(3),
            Fq::zero(),
            Fq::zero(),
            Fq::zero(),
            Fq::zero(),
            Fq::zero(),
            Fq::zero(),
            Fq::zero(),
            Fq::zero(),
            Fq::zero(),
            Fq::zero(),
            Fq::zero(),
            Fq::zero(),
            Fq::from(15),
        ];

        let expected: Vec<Fq> = vec![
            MontFp!("1"),
            MontFp!("2"),
            MontFp!("3"),
            MontFp!("9269243184842589013"),
            MontFp!("1"),
            MontFp!("2"),
            MontFp!("3"),
            MontFp!("6642849336482994628"),
            MontFp!("1"),
            MontFp!("2"),
            MontFp!("3"),
            MontFp!("9501830856585677153"),
            MontFp!("1"),
            MontFp!("2"),
            MontFp!("3"),
            MontFp!("6410261664739906488"),
        ];

        serial_frog_crt_in_place(&mut test_poly);

        dehomogenize_fq4(&mut test_poly);

        assert_eq!(test_poly, expected);
    }

    #[test]
    fn test_crt2() {
        // 2342 + 543543 * X + 3 * X^2 + 325*X^3 + 235325325 * X^5 + 765568568 * X^6
        let mut test_poly: Vec<Fq> = vec![
            MontFp!("2342"),
            MontFp!("543543"),
            MontFp!("3"),
            MontFp!("325"),
            MontFp!("0"),
            MontFp!("235325325"),
            MontFp!("765568568"),
        ];

        test_poly.resize_with(D, Fq::zero);

        let expected: Vec<Fq> = vec![
            MontFp!("2342"),
            MontFp!("843289782635822351"),
            MontFp!("9218688620283687143"),
            MontFp!("325"),
            MontFp!("2342"),
            MontFp!("15068802738690848376"),
            MontFp!("6693403901041896504"),
            MontFp!("325"),
            MontFp!("2342"),
            MontFp!("12113166087288599489"),
            MontFp!("3527640652310596771"),
            MontFp!("325"),
            MontFp!("2342"),
            MontFp!("3798926434038071238"),
            MontFp!("12384451869014986876"),
            MontFp!("325"),
        ];

        serial_frog_crt_in_place(&mut test_poly);

        dehomogenize_fq4(&mut test_poly);

        assert_eq!(test_poly, expected);
    }

    #[test]
    fn test_icrt() {
        // 1 + 2 * X + 3 * X^2 + 15 * X^15 + X^23
        let expected: Vec<Fq> = vec![
            Fq::from(1),
            Fq::from(2),
            Fq::from(3),
            Fq::zero(),
            Fq::zero(),
            Fq::zero(),
            Fq::zero(),
            Fq::zero(),
            Fq::zero(),
            Fq::zero(),
            Fq::zero(),
            Fq::zero(),
            Fq::zero(),
            Fq::zero(),
            Fq::zero(),
            Fq::from(15),
        ];

        let mut evaluations: Vec<Fq> = vec![
            MontFp!("1"),
            MontFp!("2"),
            MontFp!("3"),
            MontFp!("9269243184842589013"),
            MontFp!("1"),
            MontFp!("2"),
            MontFp!("3"),
            MontFp!("6642849336482994628"),
            MontFp!("1"),
            MontFp!("2"),
            MontFp!("3"),
            MontFp!("9501830856585677153"),
            MontFp!("1"),
            MontFp!("2"),
            MontFp!("3"),
            MontFp!("6410261664739906488"),
        ];

        homogenize_fq4(&mut evaluations);

        serial_frog_icrt_in_place(&mut evaluations);

        assert_eq!(evaluations, expected);
    }

    #[test]
    fn test_icrt_2() {
        // 2342 + 543543 * X + 3 * X^2 + 325*X^3 + 235325325 * X^5 + 765568568 * X^6
        let mut expected: Vec<Fq> = vec![
            MontFp!("2342"),
            MontFp!("543543"),
            MontFp!("3"),
            MontFp!("325"),
            MontFp!("0"),
            MontFp!("235325325"),
            MontFp!("765568568"),
        ];

        expected.resize_with(D, Fq::zero);

        let mut evaluations: Vec<Fq> = vec![
            MontFp!("2342"),
            MontFp!("843289782635822351"),
            MontFp!("9218688620283687143"),
            MontFp!("325"),
            MontFp!("2342"),
            MontFp!("15068802738690848376"),
            MontFp!("6693403901041896504"),
            MontFp!("325"),
            MontFp!("2342"),
            MontFp!("12113166087288599489"),
            MontFp!("3527640652310596771"),
            MontFp!("325"),
            MontFp!("2342"),
            MontFp!("3798926434038071238"),
            MontFp!("12384451869014986876"),
            MontFp!("325"),
        ];

        homogenize_fq4(&mut evaluations);

        serial_frog_icrt_in_place(&mut evaluations);

        assert_eq!(evaluations, expected);
    }

    fn test_crt_icrt() {
        let mut rng = ChaCha8Rng::seed_from_u64(0);
        let coefficients: Vec<Fq> = (0..16).map(|_| Fq::rand(&mut rng)).collect();

        let mut input = coefficients.clone();

        serial_frog_crt_in_place(&mut input);
        serial_frog_icrt_in_place(&mut input);

        assert_eq!(coefficients, input);
    }

    #[test]
    fn test_crt_icrt_1000000_times() {
        for _i in 0..1000000 {
            test_crt_icrt();
        }
    }

    #[test]
    fn test_crt_non_inplace() {
        let mut rng = ChaCha8Rng::seed_from_u64(0);
        let mut coefficients: Vec<Fq> = (0..16).map(|_| Fq::rand(&mut rng)).collect();

        let non_in_place = serial_frog_crt(coefficients.clone());

        serial_frog_crt_in_place(&mut coefficients);

        for i in 0..N {
            assert_eq!(
                non_in_place[i],
                Fq4::new(
                    Fq2::new(coefficients[4 * i], coefficients[4 * i + 1]),
                    Fq2::new(coefficients[4 * i + 2], coefficients[4 * i + 3])
                )
            );
        }
    }

    #[test]
    fn test_icrt_non_inplace() {
        let mut rng = ChaCha8Rng::seed_from_u64(0);
        let mut evaluations: Vec<Fq> = (0..16).map(|_| Fq::rand(&mut rng)).collect();

        let mut non_in_place: Vec<Fq4> = vec![Fq4::ZERO; N];

        for i in 0..N {
            non_in_place[i] = Fq4::new(
                Fq2::new(evaluations[4 * i], evaluations[4 * i + 1]),
                Fq2::new(evaluations[4 * i + 2], evaluations[4 * i + 3]),
            );
        }

        let non_in_place = serial_frog_icrt(non_in_place);

        serial_frog_icrt_in_place(&mut evaluations);

        assert_eq!(non_in_place, evaluations);
    }
}
