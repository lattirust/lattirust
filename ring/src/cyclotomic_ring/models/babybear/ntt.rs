//!
//! A CRT implementation for the ring Fq[X]/(X^72 - X^36 + 1).
//!
use ark_ff::BigInt;

use crate::cyclotomic_ring::models::babybear::fq9::fq9_vec_to_fq_vec;

use super::{fq9::fq_vec_to_fq9_vec, Fq, Fq9};

// The degree of the cyclotomic polynomial
pub(super) const D: usize = 72;

/// The number of splits of the cyclotomic ring in the CRT-form
pub(super) const N: usize = 8;

/// All the 24 roots of unity of degree 24.
const ROOTS_OF_UNITY_24: &[Fq] = &[
    Fq::new(BigInt([1u64])),          // power = 0
    Fq::new(BigInt([503591070u64])),  // power = 1
    Fq::new(BigInt([782862608u64])),  // power = 2
    Fq::new(BigInt([1592366214u64])), // power = 3
    Fq::new(BigInt([1314723124u64])), // power = 4
    Fq::new(BigInt([715314264u64])),  // power = 5
    Fq::new(BigInt([1728404513u64])), // power = 6
    Fq::new(BigInt([1398021245u64])), // power = 7
    Fq::new(BigInt([1314723123u64])), // power = 8
    Fq::new(BigInt([211723194u64])),  // power = 9
    Fq::new(BigInt([945541905u64])),  // power = 10
    Fq::new(BigInt([1818920952u64])), // power = 11
    Fq::new(BigInt([2013265920u64])), // power = 12
    Fq::new(BigInt([1509674851u64])), // power = 13
    Fq::new(BigInt([1230403313u64])), // power = 14
    Fq::new(BigInt([420899707u64])),  // power = 15
    Fq::new(BigInt([698542797u64])),  // power = 16
    Fq::new(BigInt([1297951657u64])), // power = 18
    Fq::new(BigInt([284861408u64])),  // power = 17
    Fq::new(BigInt([615244676u64])),  // power = 19
    Fq::new(BigInt([698542798u64])),  // power = 20
    Fq::new(BigInt([1801542727u64])), // power = 21
    Fq::new(BigInt([1067724016u64])), // power = 22
    Fq::new(BigInt([194344969u64])),  // power = 23
];

/// Given `coefficients` of a polynoimial `f mod X^72 - X^36 + 1`
/// returns its CRT:
/// *  `f mod X^9-NONRESIDUE`,
/// *  `f mod X^9-NONRESIDUE^13`,
/// *  `f mod X^9-NONRESIDUE^7`,
/// *  `f mod X^9-NONRESIDUE^19`,
/// *  `f mod X^9-NONRESIDUE^5`,
/// *  `f mod X^9-NONRESIDUE^17`,
/// *  `f mod X^9-NONRESIDUE^11`,
/// *  `f mod X^9-NONRESIDUE^23`.
///
/// Each of the components is transformed into an element of `Fq9`
/// by the corresponding unique isomorphism `Fq[X]/(X^9-NONRESIDUE^i) -> Fq[X]/(X^9-NONRESIDUE)`.
///
/// # Panics
///
/// Panics if `coefficients.len() != 24`.
#[inline(always)]
pub fn babybear_crt(coefficients: Vec<Fq>) -> Vec<Fq9> {
    // TODO: once we have parallelization set up
    //       this should be either parallel or serial.
    serial_babybear_crt(coefficients)
}

/// Same as `babybear_crt` but performs the CRT in place.
/// Takes coefficients of a polynomial and outputs
/// components of the corresponding CRT factors in each
/// nonuple from `coefficients[9*i]` to `coefficients[9*i + 8]`.
///
/// # Panics
///
/// Panics if `coefficients.len() != 24`.
#[inline(always)]
pub fn babybear_crt_in_place(coefficients: &mut [Fq]) {
    // TODO: once we have parallelization set up
    //       this should be either parallel or serial.
    serial_babybear_crt_in_place(coefficients)
}

/// The inverse CRT.
/// Takes the CRT representation in the order:
/// *  `f mod X^9-NONRESIDUE`,
/// *  `f mod X^9-NONRESIDUE^13`,
/// *  `f mod X^9-NONRESIDUE^7`,
/// *  `f mod X^9-NONRESIDUE^19`,
/// *  `f mod X^9-NONRESIDUE^5`,
/// *  `f mod X^9-NONRESIDUE^17`,
/// *  `f mod X^9-NONRESIDUE^11`,
/// *  `f mod X^9-NONRESIDUE^23`.
///
/// Each of the components is in its isomorphic form in the `Fq[X]/(X^9-NONRESIDUE)`.
/// Returns the coefficients of the polynomial encoded by this CRT form.
///
/// # Panics
///
/// Panics if `evaluations.len() != 8`.
#[inline(always)]
pub fn babybear_icrt(evaluations: Vec<Fq9>) -> Vec<Fq> {
    // TODO: once we have parallelization set up
    //       this should be either parallel or serial.
    serial_babybear_icrt(evaluations)
}

/// Same as `babybear_icrt` but performs the inverse CRT in place.
/// Each nonuple from `evaluations[9*i]` to `evaluations[9*i + 8]`,
/// has to be an `Fq9` element. In the order described in `babybear_icrt`'s docstring.
///
/// # Panics
///
/// Panics if `evaluations.len() != 24`.
#[inline(always)]
pub fn babybear_icrt_in_place(evaluations: &mut [Fq]) {
    // TODO: once we have parallelization set up
    //       this should be either parallel or serial.
    serial_babybear_icrt_in_place(evaluations)
}

fn serial_babybear_crt(mut coefficients: Vec<Fq>) -> Vec<Fq9> {
    assert_eq!(coefficients.len(), D);
    serial_babybear_crt_in_place(&mut coefficients);
    fq_vec_to_fq9_vec(coefficients)
}

fn serial_babybear_icrt(evaluations: Vec<Fq9>) -> Vec<Fq> {
    assert_eq!(evaluations.len(), N);
    let mut evaluations = fq9_vec_to_fq_vec(evaluations);
    serial_babybear_icrt_in_place(&mut evaluations);
    evaluations
}

/// (2 * ROOT_OF_UNITY_24[4] - 1)^-1
const KAPPA: Fq = Fq::new(BigInt([1807872479]));
/// 1 / 8
const EIGHT_INV: Fq = Fq::new(BigInt([1761607681]));
/// 1 / 4
const FOUR_INV: Fq = Fq::new(BigInt([1509949441]));

fn serial_babybear_crt_in_place(coefficients: &mut [Fq]) {
    assert_eq!(coefficients.len(), D);

    // Compute f mod X^36-zeta and f mod X^36 - zeta^5
    // Use the technique from https://eprint.iacr.org/2019/040.
    // We take zeta = ROOT_OF_UNITY[4] i.e. the sixth primitive root of unity.
    // Then we can do the reduction using the formulas
    // f_0_i = f_i + zeta f_{i + n/2},
    // f_1_i = f_i + f_{i+n/2} - zeta f_{i + n/2},
    // Since X^72-X^36+1 = (X^36 - zeta) * (X^36 - zeta^5) in Fq
    // and zeta^5 = 1 - zeta.
    for i in 0..(D / 2) {
        let (coeff_i, coeff_d_div_2_plus_i) = (coefficients[i], coefficients[D / 2 + i]);
        let zeta_coeff_d_div_2_plus_i = ROOTS_OF_UNITY_24[4] * coeff_d_div_2_plus_i;

        coefficients[i] = coeff_i + zeta_coeff_d_div_2_plus_i;
        coefficients[i + D / 2] = coeff_i + coeff_d_div_2_plus_i - zeta_coeff_d_div_2_plus_i;
    }
    // After this step we get f_1 in the first half coefficients and f_2 in the second half.

    // From here we can perform radix-2 CRT.

    // Compute f_1 mod X^18-\sigma, f_1 mod X^18-\sigma^7, f_2 mod X^18 -\sigma^5, f_2 mod X^18 -\sigma^11
    for i in 0..(D / 4) {
        // f_1
        {
            let (coeff_i, coeff_d_div_4_plus_i) = (coefficients[i], coefficients[D / 4 + i]);
            let sigma_coeff_d_div_4_plus_i = ROOTS_OF_UNITY_24[2] * coeff_d_div_4_plus_i;

            coefficients[i] = coeff_i + sigma_coeff_d_div_4_plus_i;
            coefficients[i + D / 4] = coeff_i - sigma_coeff_d_div_4_plus_i;
        }

        // f_2
        {
            let (coeff_i, coeff_d_div_4_plus_i) =
                (coefficients[D / 2 + i], coefficients[3 * D / 4 + i]);
            let sigma_coeff_d_div_4_plus_i = ROOTS_OF_UNITY_24[10] * coeff_d_div_4_plus_i;

            coefficients[D / 2 + i] = coeff_i + sigma_coeff_d_div_4_plus_i;
            coefficients[3 * D / 4 + i] = coeff_i - sigma_coeff_d_div_4_plus_i;
        }
    }
    // After this step we have f_1, f_2, f_3, f_4.

    // Compute f_1 mod X^9-NONRESIDUE, f_1 mod X^9-NONRESIDUE^13, f_2 mod X^9-NONRESIDUE^7, f_2 mod X^9-NONRESIDUE^19
    // Compute f_3 mod X^9-NONRESIDUE^5, f_3 mod X^9-NONRESIDUE^17, f_4 mod X^9-NONRESIDUE^11, f_4 mod X^9-NONRESIDUE^23
    for i in 0..(D / 8) {
        // f_1
        {
            let (coeff_i, coeff_d_div_8_plus_i) = (coefficients[i], coefficients[D / 8 + i]);
            let nonresidue_coeff_d_div_8_plus_i = ROOTS_OF_UNITY_24[1] * coeff_d_div_8_plus_i;

            coefficients[i] = coeff_i + nonresidue_coeff_d_div_8_plus_i;
            coefficients[i + D / 8] = coeff_i - nonresidue_coeff_d_div_8_plus_i;
        }

        // f_2
        {
            let (coeff_i, coeff_d_div_8_plus_i) =
                (coefficients[D / 4 + i], coefficients[3 * D / 8 + i]);
            let sigma_coeff_d_div_8_plus_i = ROOTS_OF_UNITY_24[7] * coeff_d_div_8_plus_i;

            coefficients[D / 4 + i] = coeff_i + sigma_coeff_d_div_8_plus_i;
            coefficients[3 * D / 8 + i] = coeff_i - sigma_coeff_d_div_8_plus_i;
        }

        // f_3
        {
            let (coeff_i, coeff_d_div_8_plus_i) =
                (coefficients[D / 2 + i], coefficients[5 * D / 8 + i]);
            let sigma_coeff_d_div_8_plus_i = ROOTS_OF_UNITY_24[5] * coeff_d_div_8_plus_i;

            coefficients[D / 2 + i] = coeff_i + sigma_coeff_d_div_8_plus_i;
            coefficients[5 * D / 8 + i] = coeff_i - sigma_coeff_d_div_8_plus_i;
        }

        // f_4
        {
            let (coeff_i, coeff_d_div_8_plus_i) =
                (coefficients[3 * D / 4 + i], coefficients[7 * D / 8 + i]);
            let sigma_coeff_d_div_8_plus_i = ROOTS_OF_UNITY_24[11] * coeff_d_div_8_plus_i;

            coefficients[3 * D / 4 + i] = coeff_i + sigma_coeff_d_div_8_plus_i;
            coefficients[7 * D / 8 + i] = coeff_i - sigma_coeff_d_div_8_plus_i;
        }
    }

    homogenize_fq9(coefficients);
}

fn serial_babybear_icrt_in_place(evaluations: &mut [Fq]) {
    assert_eq!(evaluations.len(), D);

    dehomogenize_fq9(evaluations);

    // Given f_1 mod X^9-NONRESIDUE, f_1 mod X^9-NONRESIDUE^13, f_2 mod X^9-NONRESIDUE^7, f_2 mod X^9-NONRESIDUE^19,
    //       f_3 mod X^9-NONRESIDUE^5, f_3 mod X^9-NONRESIDUE^17, f_4 mod X^9-NONRESIDUE^11, f_4 mod X^9-NONRESIDUE^23.
    // recreate f_1 mod X^18-\sigma, f_1 mod X^18-\sigma^7, f_2 mod X^18 -\sigma^5, f_2 mod X^18 -\sigma^11, where \sigma=NONRESIDUE ^ 2.
    for i in 0..(D / 8) {
        // f_1
        {
            let (coeff_i, coeff_d_div_8_plus_i) = (evaluations[i], evaluations[D / 8 + i]);

            evaluations[i] = coeff_i + coeff_d_div_8_plus_i;
            evaluations[i + D / 8] = ROOTS_OF_UNITY_24[23] * (coeff_i - coeff_d_div_8_plus_i);
        }

        // f_2
        {
            let (coeff_i, coeff_d_div_8_plus_i) =
                (evaluations[D / 4 + i], evaluations[3 * D / 8 + i]);

            evaluations[D / 4 + i] = coeff_i + coeff_d_div_8_plus_i;
            evaluations[3 * D / 8 + i] = ROOTS_OF_UNITY_24[17] * (coeff_i - coeff_d_div_8_plus_i);
        }

        // f_3
        {
            let (coeff_i, coeff_d_div_8_plus_i) =
                (evaluations[D / 2 + i], evaluations[5 * D / 8 + i]);

            evaluations[D / 2 + i] = coeff_i + coeff_d_div_8_plus_i;
            evaluations[5 * D / 8 + i] = ROOTS_OF_UNITY_24[19] * (coeff_i - coeff_d_div_8_plus_i);
        }

        // f_4
        {
            let (coeff_i, coeff_d_div_8_plus_i) =
                (evaluations[3 * D / 4 + i], evaluations[7 * D / 8 + i]);

            evaluations[3 * D / 4 + i] = coeff_i + coeff_d_div_8_plus_i;
            evaluations[7 * D / 8 + i] = ROOTS_OF_UNITY_24[13] * (coeff_i - coeff_d_div_8_plus_i);
        }
    }

    // Given f_1 mod X^18-\sigma, f_1 mod X^18-\sigma^7, f_2 mod X^18 -\sigma^5,
    //       f_2 mod X^18 -\sigma^11, where \sigma=NONRESIDUE ^ 2,
    // recreate f_1 mod X^18-\sigma, f_1 mod X^18-\sigma^7, f_2 mod X^18 -\sigma^5, f_2 mod X^18 -\sigma^11, where \sigma=NONRESIDUE ^ 2.
    for i in 0..(D / 4) {
        // f_1
        {
            let (coeff_i, coeff_d_div_4_plus_i) = (evaluations[i], evaluations[D / 4 + i]);

            evaluations[i] = coeff_i + coeff_d_div_4_plus_i;
            evaluations[i + D / 4] = ROOTS_OF_UNITY_24[22] * (coeff_i - coeff_d_div_4_plus_i);
        }

        // f_2
        {
            let (coeff_i, coeff_d_div_4_plus_i) =
                (evaluations[D / 2 + i], evaluations[3 * D / 4 + i]);

            evaluations[D / 2 + i] = coeff_i + coeff_d_div_4_plus_i;
            evaluations[3 * D / 4 + i] = ROOTS_OF_UNITY_24[14] * (coeff_i - coeff_d_div_4_plus_i);
        }
    }

    // Rewind the first step of the CRT.
    for i in 0..(D / 2) {
        let (coeff_i, coeff_d_div_2_plus_i) = (evaluations[i], evaluations[D / 2 + i]);
        let kappa_diff = KAPPA * (coeff_i - coeff_d_div_2_plus_i);

        evaluations[i] = EIGHT_INV * (coeff_i + coeff_d_div_2_plus_i - kappa_diff);
        evaluations[i + D / 2] = FOUR_INV * kappa_diff;
    }
}

/// At the end of CRT we get elements in different (although isomorphic)
/// degree-9 extensions of Fq,
/// this function converts each triple `c[9 * i]` to `c[9 * {i + 1) - 1]`
/// into coefficients of an element from Fq[X]/(X^9-NONRESIDUE)=Fq9.
#[inline(always)]
fn homogenize_fq9(c: &mut [Fq]) {
    permute_to_fq9_of_fq3(&mut c[0..9]);
    nonresidue_to_13_to_nonresidue(&mut c[9..18]);
    nonresidue_to_7_to_nonresidue(&mut c[18..27]);
    nonresidue_to_19_to_nonresidue(&mut c[27..36]);
    nonresidue_to_5_to_nonresidue(&mut c[36..45]);
    nonresidue_to_17_to_nonresidue(&mut c[45..54]);
    nonresidue_to_11_to_nonresidue(&mut c[54..63]);
    nonresidue_to_23_to_nonresidue(&mut c[63..72]);
}

/// The inverse of the above.
#[inline(always)]
fn dehomogenize_fq9(c: &mut [Fq]) {
    inverse_permute_to_fq9_of_fq3(&mut c[0..9]);
    nonresidue_to_nonresidue_to_13(&mut c[9..18]);
    nonresidue_to_nonresidue_to_7(&mut c[18..27]);
    nonresidue_to_nonresidue_to_19(&mut c[27..36]);
    nonresidue_to_nonresidue_to_5(&mut c[36..45]);
    nonresidue_to_nonresidue_to_17(&mut c[45..54]);
    nonresidue_to_nonresidue_to_11(&mut c[54..63]);
    nonresidue_to_nonresidue_to_23(&mut c[63..72]);
}

// Different automorphisms with the target Fp(NONRESIDUE) and their inverses.
// TODO: Double check automorphism
#[inline(always)]
fn nonresidue_to_13_to_nonresidue(c: &mut [Fq]) {
    let c1 = c[1];
    c[1] = c[7] * ROOTS_OF_UNITY_24[10];
    c[7] = c[4] * ROOTS_OF_UNITY_24[5];
    c[4] = c1 * ROOTS_OF_UNITY_24[1];

    let c2 = c[2];
    c[2] = c[5] * ROOTS_OF_UNITY_24[7];
    c[5] = c[8] * ROOTS_OF_UNITY_24[11];
    c[8] = c2 * ROOTS_OF_UNITY_24[2];

    c[3] *= ROOTS_OF_UNITY_24[4];
    c[6] *= ROOTS_OF_UNITY_24[8];
    permute_to_fq9_of_fq3(c);
}

#[inline(always)]
fn nonresidue_to_nonresidue_to_13(c: &mut [Fq]) {
    inverse_permute_to_fq9_of_fq3(c);
    let c1 = c[1];
    c[1] = c[4] * ROOTS_OF_UNITY_24[23];
    c[4] = c[7] * ROOTS_OF_UNITY_24[19];
    c[7] = c1 * ROOTS_OF_UNITY_24[14];

    let c2 = c[2];
    c[2] = c[8] * ROOTS_OF_UNITY_24[22];
    c[8] = c[5] * ROOTS_OF_UNITY_24[13];
    c[5] = c2 * ROOTS_OF_UNITY_24[17];

    c[3] *= ROOTS_OF_UNITY_24[20];
    c[6] *= ROOTS_OF_UNITY_24[16];
}

#[inline(always)]
fn nonresidue_to_7_to_nonresidue(c: &mut [Fq]) {
    let c1 = c[1];
    c[1] = c[4] * ROOTS_OF_UNITY_24[3];
    c[4] = c[7] * ROOTS_OF_UNITY_24[5];
    c[7] = c1;

    let c2 = c[2];
    c[2] = c[8] * ROOTS_OF_UNITY_24[6];
    c[8] = c[5] * ROOTS_OF_UNITY_24[3];
    c[5] = c2 * ROOTS_OF_UNITY_24[1];

    c[3] *= ROOTS_OF_UNITY_24[2];
    c[6] *= ROOTS_OF_UNITY_24[4];
    permute_to_fq9_of_fq3(c);
}

#[inline(always)]
fn nonresidue_to_nonresidue_to_7(c: &mut [Fq]) {
    inverse_permute_to_fq9_of_fq3(c);
    let c1 = c[1];
    c[1] = c[7];
    c[7] = c[4] * ROOTS_OF_UNITY_24[19];
    c[4] = c1 * ROOTS_OF_UNITY_24[21];

    let c2 = c[2];
    c[2] = c[5] * ROOTS_OF_UNITY_24[23];
    c[5] = c[8] * ROOTS_OF_UNITY_24[21];
    c[8] = c2 * ROOTS_OF_UNITY_24[18];

    c[3] *= ROOTS_OF_UNITY_24[22];
    c[6] *= ROOTS_OF_UNITY_24[20];
}

#[inline(always)]
fn nonresidue_to_19_to_nonresidue(c: &mut [Fq]) {
    c[1] *= ROOTS_OF_UNITY_24[2];
    c[2] *= ROOTS_OF_UNITY_24[4];
    c[3] *= ROOTS_OF_UNITY_24[6];
    c[4] *= ROOTS_OF_UNITY_24[8];
    c[5] *= ROOTS_OF_UNITY_24[10];
    c[6] = -c[6];
    c[7] *= ROOTS_OF_UNITY_24[14];
    c[8] *= ROOTS_OF_UNITY_24[16];
    permute_to_fq9_of_fq3(c);
}

#[inline(always)]
fn nonresidue_to_nonresidue_to_19(c: &mut [Fq]) {
    inverse_permute_to_fq9_of_fq3(c);
    c[1] *= ROOTS_OF_UNITY_24[22];
    c[2] *= ROOTS_OF_UNITY_24[20];
    c[3] *= ROOTS_OF_UNITY_24[18];
    c[4] *= ROOTS_OF_UNITY_24[16];
    c[5] *= ROOTS_OF_UNITY_24[14];
    c[6] = -c[6];
    c[7] *= ROOTS_OF_UNITY_24[10];
    c[8] *= ROOTS_OF_UNITY_24[8];
}

#[inline(always)]
fn nonresidue_to_5_to_nonresidue(c: &mut [Fq]) {
    let c1 = c[1];
    c[1] = c[2] * ROOTS_OF_UNITY_24[1];
    c[2] = c[4] * ROOTS_OF_UNITY_24[2];
    c[4] = c[8] * ROOTS_OF_UNITY_24[4];
    c[8] = c[7] * ROOTS_OF_UNITY_24[3];
    c[7] = c[5] * ROOTS_OF_UNITY_24[2];
    c[5] = c1;

    let c3 = c[3];
    c[3] = c[6] * ROOTS_OF_UNITY_24[3];
    c[6] = c3 * ROOTS_OF_UNITY_24[1];
    permute_to_fq9_of_fq3(c);
}

#[inline(always)]
fn nonresidue_to_nonresidue_to_5(c: &mut [Fq]) {
    inverse_permute_to_fq9_of_fq3(c);
    let c1 = c[1];
    c[1] = c[5];
    c[5] = c[7] * ROOTS_OF_UNITY_24[22];
    c[7] = c[8] * ROOTS_OF_UNITY_24[21];
    c[8] = c[4] * ROOTS_OF_UNITY_24[20];
    c[4] = c[2] * ROOTS_OF_UNITY_24[22];
    c[2] = c1 * ROOTS_OF_UNITY_24[23];

    let c3 = c[3];
    c[3] = c[6] * ROOTS_OF_UNITY_24[23];
    c[6] = c3 * ROOTS_OF_UNITY_24[21];
}

#[inline(always)]
fn nonresidue_to_17_to_nonresidue(c: &mut [Fq]) {
    let c1 = c[1];
    c[1] = c[8] * ROOTS_OF_UNITY_24[15];
    c[8] = c1 * ROOTS_OF_UNITY_24[1];

    let c2 = c[2];
    c[2] = c[7] * ROOTS_OF_UNITY_24[13];
    c[7] = c2 * ROOTS_OF_UNITY_24[3];

    let c3 = c[3];
    c[3] = c[6] * ROOTS_OF_UNITY_24[11];
    c[6] = c3 * ROOTS_OF_UNITY_24[5];

    let c4 = c[4];
    c[4] = c[5] * ROOTS_OF_UNITY_24[9];
    c[5] = c4 * ROOTS_OF_UNITY_24[7];
    permute_to_fq9_of_fq3(c);
}

#[inline(always)]
fn nonresidue_to_nonresidue_to_17(c: &mut [Fq]) {
    inverse_permute_to_fq9_of_fq3(c);
    let c1 = c[1];
    c[1] = c[8] * ROOTS_OF_UNITY_24[23];
    c[8] = c1 * ROOTS_OF_UNITY_24[9];

    let c2 = c[2];
    c[2] = c[7] * ROOTS_OF_UNITY_24[21];
    c[7] = c2 * ROOTS_OF_UNITY_24[11];

    let c3 = c[3];
    c[3] = c[6] * ROOTS_OF_UNITY_24[19];
    c[6] = c3 * ROOTS_OF_UNITY_24[13];

    let c4 = c[4];
    c[4] = c[5] * ROOTS_OF_UNITY_24[17];
    c[5] = c4 * ROOTS_OF_UNITY_24[15];
}

#[inline(always)]
fn nonresidue_to_11_to_nonresidue(c: &mut [Fq]) {
    let c1 = c[1];
    c[1] = c[5] * ROOTS_OF_UNITY_24[6];
    c[5] = c[7] * ROOTS_OF_UNITY_24[8];
    c[7] = c[8] * ROOTS_OF_UNITY_24[9];
    c[8] = c[4] * ROOTS_OF_UNITY_24[4];
    c[4] = c[2] * ROOTS_OF_UNITY_24[2];
    c[2] = c1 * ROOTS_OF_UNITY_24[1];

    let c3 = c[3];
    c[3] = c[6] * ROOTS_OF_UNITY_24[7];
    c[6] = c3 * ROOTS_OF_UNITY_24[3];
    permute_to_fq9_of_fq3(c);
}

#[inline(always)]
fn nonresidue_to_nonresidue_to_11(c: &mut [Fq]) {
    inverse_permute_to_fq9_of_fq3(c);
    let c1 = c[1];
    c[1] = c[2] * ROOTS_OF_UNITY_24[23];
    c[2] = c[4] * ROOTS_OF_UNITY_24[22];
    c[4] = c[8] * ROOTS_OF_UNITY_24[20];
    c[8] = c[7] * ROOTS_OF_UNITY_24[15];
    c[7] = c[5] * ROOTS_OF_UNITY_24[16];
    c[5] = c1 * ROOTS_OF_UNITY_24[18];

    let c3 = c[3];
    c[3] = c[6] * ROOTS_OF_UNITY_24[21];
    c[6] = c3 * ROOTS_OF_UNITY_24[17];
}

#[inline(always)]
fn nonresidue_to_23_to_nonresidue(c: &mut [Fq]) {
    let c1 = c[1];
    c[1] = c[2] * ROOTS_OF_UNITY_24[5];
    c[2] = c[4] * ROOTS_OF_UNITY_24[10];
    c[4] = c[8] * ROOTS_OF_UNITY_24[20];
    c[8] = c[7] * ROOTS_OF_UNITY_24[17];
    c[7] = -c[5];
    c[5] = c1 * ROOTS_OF_UNITY_24[2];

    let c3 = c[3];
    c[3] = c[6] * ROOTS_OF_UNITY_24[15];
    c[6] = c3 * ROOTS_OF_UNITY_24[7];
    permute_to_fq9_of_fq3(c);
}

#[inline(always)]
fn nonresidue_to_nonresidue_to_23(c: &mut [Fq]) {
    inverse_permute_to_fq9_of_fq3(c);
    let c1 = c[1];
    c[1] = c[5] * ROOTS_OF_UNITY_24[22];
    c[5] = -c[7];
    c[7] = c[8] * ROOTS_OF_UNITY_24[7];
    c[8] = c[4] * ROOTS_OF_UNITY_24[4];
    c[4] = c[2] * ROOTS_OF_UNITY_24[14];
    c[2] = c1 * ROOTS_OF_UNITY_24[19];

    let c3 = c[3];
    c[3] = c[6] * ROOTS_OF_UNITY_24[17];
    c[6] = c3 * ROOTS_OF_UNITY_24[9];
}

const SWAPS: [(usize, usize); 3] = [(1, 3), (2, 6), (5, 7)];

fn permute_to_fq9_of_fq3(fq9_vec: &mut [Fq]) {
    for &(i, j) in SWAPS.iter() {
        fq9_vec.swap(i, j);
    }
}

fn inverse_permute_to_fq9_of_fq3(fq9_vec: &mut [Fq]) {
    for &(i, j) in SWAPS.iter().rev() {
        fq9_vec.swap(i, j);
    }
}

#[cfg(test)]
mod tests {
    use crate::cyclotomic_ring::models::babybear::BabyBear3ExtConfig;

    use super::*;
    use ark_ff::{Field, Fp3Config, MontFp, UniformRand};
    use ark_poly::{
        univariate::{DenseOrSparsePolynomial, DensePolynomial, SparsePolynomial},
        DenseUVPolynomial,
    };
    use rand::thread_rng;

    #[test]
    fn test_babybear_roots_of_unity() {
        // Check if they all are roots of unity of degree 24.
        for root in ROOTS_OF_UNITY_24 {
            assert_eq!(root.pow([24]), Fq::ONE);
        }

        // Check if they are pairwise distinct.
        for (i, item) in ROOTS_OF_UNITY_24.iter().enumerate().take(24) {
            for item2 in ROOTS_OF_UNITY_24.iter().take(24).skip(i + 1) {
                assert_ne!(*item, *item2);
            }
        }

        // Check if they come in order.
        for i in 0..24u64 {
            assert_eq!(ROOTS_OF_UNITY_24[i as usize], ROOTS_OF_UNITY_24[1].pow([i]));
        }
    }

    macro_rules! test_inverses {
        ($($isomorphism:expr, $inverse:expr),+) => {
            let mut rng = thread_rng();
            $({
                let x = Fq9::rand(&mut rng);
                let mut x_prime = fq9_vec_to_fq_vec(vec![x]);
                $isomorphism(&mut x_prime);
                $inverse(&mut x_prime);
                assert_eq!(
                    fq_vec_to_fq9_vec(x_prime).first().cloned().unwrap(),
                    x
                );
            })+
        };
    }

    #[test]
    fn test_babybear_isomorphism_inverses() {
        test_inverses! {
            nonresidue_to_nonresidue_to_13, nonresidue_to_13_to_nonresidue,
            nonresidue_to_nonresidue_to_7, nonresidue_to_7_to_nonresidue,
            nonresidue_to_nonresidue_to_19, nonresidue_to_19_to_nonresidue,
            nonresidue_to_nonresidue_to_5, nonresidue_to_5_to_nonresidue,
            nonresidue_to_nonresidue_to_17, nonresidue_to_17_to_nonresidue,
            nonresidue_to_nonresidue_to_11, nonresidue_to_11_to_nonresidue,
            nonresidue_to_nonresidue_to_23, nonresidue_to_23_to_nonresidue
        }
    }

    macro_rules! test_x_powers {
        ($($inverse:expr),+; $($i:expr),*) => {
            $({
                let mut x_prime: Vec<Fq> =
                    vec![Fq::ZERO, Fq::ONE, Fq::ZERO, Fq::ZERO, Fq::ZERO, Fq::ZERO, Fq::ZERO, Fq::ZERO, Fq::ZERO];
                let mut x_prime_2: Vec<Fq> =
                    vec![Fq::ZERO; 9];
                x_prime_2[$i] = Fq::ONE;
                $inverse(&mut x_prime);
                $inverse(&mut x_prime_2);

                let x_prime = fq_vec_to_fq9_vec(x_prime).first().cloned().unwrap();
                let x_prime_2 = fq_vec_to_fq9_vec(x_prime_2).first().cloned().unwrap();

                let mut result = x_prime;
                for _ in 1..$i {
                    result *= x_prime;
                }

                assert_eq!(result, x_prime_2);
            })+
        };
    }

    #[test]
    fn test_babybear_powers() {
        test_x_powers! {
            nonresidue_to_13_to_nonresidue,
            nonresidue_to_7_to_nonresidue,
            nonresidue_to_19_to_nonresidue,
            nonresidue_to_5_to_nonresidue,
            nonresidue_to_17_to_nonresidue,
            nonresidue_to_11_to_nonresidue,
            nonresidue_to_23_to_nonresidue;
            2,3,4,5,6,7,8
        }
    }

    macro_rules! test_extension_equations {
        ($($root_of_unity:expr, $inverse:expr),+) => {
            $({
                let mut x_prime: Vec<Fq> = vec![Fq::ZERO, Fq::ONE, Fq::ZERO, Fq::ZERO, Fq::ZERO, Fq::ZERO, Fq::ZERO, Fq::ZERO, Fq::ZERO];
                $inverse(&mut x_prime);

                let x_prime = crate::cyclotomic_ring::models::babybear::fq9::fq_vec_to_fq9_vec(x_prime)[0];

                assert_eq!(x_prime * x_prime * x_prime * x_prime * x_prime * x_prime * x_prime * x_prime * x_prime - Fq9::from_base_prime_field($root_of_unity), Fq9::ZERO);
            })+
        };
    }

    #[test]
    fn test_equations() {
        test_extension_equations! {
            ROOTS_OF_UNITY_24[13], nonresidue_to_13_to_nonresidue,
            ROOTS_OF_UNITY_24[7], nonresidue_to_7_to_nonresidue,
            ROOTS_OF_UNITY_24[19], nonresidue_to_19_to_nonresidue,
            ROOTS_OF_UNITY_24[5], nonresidue_to_5_to_nonresidue,
            ROOTS_OF_UNITY_24[17], nonresidue_to_17_to_nonresidue,
            ROOTS_OF_UNITY_24[11], nonresidue_to_11_to_nonresidue,
            ROOTS_OF_UNITY_24[23], nonresidue_to_23_to_nonresidue
        };
    }

    fn test_fq9_multiplication() {
        let mut rng = thread_rng();
        let f1 = Fq9::rand(&mut rng);
        let f2 = Fq9::rand(&mut rng);
        let f_mul = f1 * f2;
        let mut f1_vec = fq9_vec_to_fq_vec(vec![f1]);
        let mut f2_vec = fq9_vec_to_fq_vec(vec![f2]);
        f1_vec.resize(9, Fq::ZERO);
        f2_vec.resize(9, Fq::ZERO);
        inverse_permute_to_fq9_of_fq3(&mut f1_vec);
        inverse_permute_to_fq9_of_fq3(&mut f2_vec);
        let poly1 = DensePolynomial::from_coefficients_slice(&f1_vec);
        let poly2 = DensePolynomial::from_coefficients_slice(&f2_vec);
        let poly_mul = &poly1 * &poly2;
        let minimal_poly = SparsePolynomial::from_coefficients_slice(&[
            (9, Fq::ONE),
            (0, -BabyBear3ExtConfig::NONRESIDUE),
        ]);

        let (_, r) =
            DenseOrSparsePolynomial::divide_with_q_and_r(&poly_mul.into(), &minimal_poly.into())
                .unwrap();
        let mut f_mul_vec = r.coeffs().to_vec();
        permute_to_fq9_of_fq3(&mut f_mul_vec);
        let f_mul_vec = fq_vec_to_fq9_vec(f_mul_vec);
        assert_eq!(f_mul, f_mul_vec.first().cloned().unwrap());
    }

    #[test]
    fn test_fq9_multiplication_1000000_times() {
        for _ in 0..1000000 {
            test_fq9_multiplication()
        }
    }

    #[test]
    fn test_normalize_denormalize() {
        let mut rng = thread_rng();
        let mut x: Vec<Fq> = (0..72).map(|_| Fq::rand(&mut rng)).collect();

        let expected = x.clone();

        homogenize_fq9(&mut x);
        dehomogenize_fq9(&mut x);

        assert_eq!(x, expected);
    }

    fn test_babybear_crt() {
        let mut rng = thread_rng();
        let mut test_poly = (0..D).map(|_| Fq::rand(&mut rng)).collect::<Vec<_>>();

        let poly = DensePolynomial::from_coefficients_slice(&test_poly);

        let r1_poly = SparsePolynomial::from_coefficients_slice(&[
            (9, <Fq as Field>::ONE),
            (0, -ROOTS_OF_UNITY_24[1]),
        ]);
        let r13_poly = SparsePolynomial::from_coefficients_slice(&[
            (9, <Fq as Field>::ONE),
            (0, -ROOTS_OF_UNITY_24[13]),
        ]);
        let r7_poly = SparsePolynomial::from_coefficients_slice(&[
            (9, <Fq as Field>::ONE),
            (0, -ROOTS_OF_UNITY_24[7]),
        ]);
        let r19_poly = SparsePolynomial::from_coefficients_slice(&[
            (9, <Fq as Field>::ONE),
            (0, -ROOTS_OF_UNITY_24[19]),
        ]);
        let r5_poly = SparsePolynomial::from_coefficients_slice(&[
            (9, <Fq as Field>::ONE),
            (0, -ROOTS_OF_UNITY_24[5]),
        ]);
        let r17_poly = SparsePolynomial::from_coefficients_slice(&[
            (9, <Fq as Field>::ONE),
            (0, -ROOTS_OF_UNITY_24[17]),
        ]);
        let r11_poly = SparsePolynomial::from_coefficients_slice(&[
            (9, <Fq as Field>::ONE),
            (0, -ROOTS_OF_UNITY_24[11]),
        ]);
        let r23_poly = SparsePolynomial::from_coefficients_slice(&[
            (9, <Fq as Field>::ONE),
            (0, -ROOTS_OF_UNITY_24[23]),
        ]);

        let (_, r1) =
            DenseOrSparsePolynomial::divide_with_q_and_r(&poly.clone().into(), &r1_poly.into())
                .unwrap();
        let (_, r13) =
            DenseOrSparsePolynomial::divide_with_q_and_r(&poly.clone().into(), &r13_poly.into())
                .unwrap();
        let (_, r7) =
            DenseOrSparsePolynomial::divide_with_q_and_r(&poly.clone().into(), &r7_poly.into())
                .unwrap();
        let (_, r19) =
            DenseOrSparsePolynomial::divide_with_q_and_r(&poly.clone().into(), &r19_poly.into())
                .unwrap();
        let (_, r5) =
            DenseOrSparsePolynomial::divide_with_q_and_r(&poly.clone().into(), &r5_poly.into())
                .unwrap();
        let (_, r17) =
            DenseOrSparsePolynomial::divide_with_q_and_r(&poly.clone().into(), &r17_poly.into())
                .unwrap();
        let (_, r11) =
            DenseOrSparsePolynomial::divide_with_q_and_r(&poly.clone().into(), &r11_poly.into())
                .unwrap();
        let (_, r23) =
            DenseOrSparsePolynomial::divide_with_q_and_r(&poly.clone().into(), &r23_poly.into())
                .unwrap();

        let mut expected: Vec<Fq> = Vec::with_capacity(D);
        let mut r1 = r1.coeffs;
        r1.resize(9, Fq::ZERO);
        let mut r13 = r13.coeffs;
        r13.resize(9, Fq::ZERO);
        let mut r7 = r7.coeffs;
        r7.resize(9, Fq::ZERO);
        let mut r19 = r19.coeffs;
        r19.resize(9, Fq::ZERO);
        let mut r5 = r5.coeffs;
        r5.resize(9, Fq::ZERO);
        let mut r17 = r17.coeffs;
        r17.resize(9, Fq::ZERO);
        let mut r11 = r11.coeffs;
        r11.resize(9, Fq::ZERO);
        let mut r23 = r23.coeffs;
        r23.resize(9, Fq::ZERO);
        expected.extend(r1);
        expected.extend(r13);
        expected.extend(r7);
        expected.extend(r19);
        expected.extend(r5);
        expected.extend(r17);
        expected.extend(r11);
        expected.extend(r23);

        babybear_crt_in_place(&mut test_poly);
        dehomogenize_fq9(&mut test_poly);
        assert_eq!(expected, test_poly);
    }

    #[test]
    fn test_babybear_crt_1000000_times() {
        for _ in 0..1000000 {
            test_babybear_crt()
        }
    }

    #[test]
    fn test_babybear_icrt_hardcoded() {
        let mut initial_ntt = vec![Fq::ZERO; D];
        initial_ntt[0] = MontFp!("1900625136");
        initial_ntt[1] = MontFp!("112939065");
        initial_ntt[2] = MontFp!("80310056");
        initial_ntt[3] = MontFp!("1982426205");
        initial_ntt[4] = MontFp!("1861467068");
        initial_ntt[5] = MontFp!("49211297");
        initial_ntt[6] = MontFp!("1445971202");
        initial_ntt[7] = MontFp!("1729640160");
        initial_ntt[8] = MontFp!("1963708676");
        initial_ntt[9] = MontFp!("954083992");
        initial_ntt[10] = MontFp!("829076309");
        initial_ntt[11] = MontFp!("438872109");
        initial_ntt[12] = MontFp!("1993465853");
        initial_ntt[13] = MontFp!("308690395");
        initial_ntt[14] = MontFp!("488768419");
        initial_ntt[15] = MontFp!("1901524572");
        initial_ntt[16] = MontFp!("1034036719");
        initial_ntt[17] = MontFp!("1241375270");
        initial_ntt[18] = MontFp!("466225138");
        initial_ntt[19] = MontFp!("140725638");
        initial_ntt[20] = MontFp!("1386731196");
        initial_ntt[21] = MontFp!("1170787115");
        initial_ntt[22] = MontFp!("837660087");
        initial_ntt[23] = MontFp!("416513957");
        initial_ntt[24] = MontFp!("1889201657");
        initial_ntt[25] = MontFp!("894740305");
        initial_ntt[26] = MontFp!("1840816808");
        initial_ntt[27] = MontFp!("1326838466");
        initial_ntt[28] = MontFp!("1166698972");
        initial_ntt[29] = MontFp!("712318050");
        initial_ntt[30] = MontFp!("921000781");
        initial_ntt[31] = MontFp!("1023319258");
        initial_ntt[32] = MontFp!("1379077144");
        initial_ntt[33] = MontFp!("501336344");
        initial_ntt[34] = MontFp!("1221812213");
        initial_ntt[35] = MontFp!("128415714");
        initial_ntt[36] = MontFp!("1320317132");
        initial_ntt[37] = MontFp!("289989777");
        initial_ntt[38] = MontFp!("323829215");
        initial_ntt[39] = MontFp!("1532388335");
        initial_ntt[40] = MontFp!("672024586");
        initial_ntt[41] = MontFp!("1004795423");
        initial_ntt[42] = MontFp!("801782387");
        initial_ntt[43] = MontFp!("630843091");
        initial_ntt[44] = MontFp!("152043073");
        initial_ntt[45] = MontFp!("1168499229");
        initial_ntt[46] = MontFp!("261644413");
        initial_ntt[47] = MontFp!("202256778");
        initial_ntt[48] = MontFp!("1733651679");
        initial_ntt[49] = MontFp!("802644602");
        initial_ntt[50] = MontFp!("1547582674");
        initial_ntt[51] = MontFp!("1544500919");
        initial_ntt[52] = MontFp!("58827862");
        initial_ntt[53] = MontFp!("1357070166");
        initial_ntt[54] = MontFp!("1162233649");
        initial_ntt[55] = MontFp!("1001050512");
        initial_ntt[56] = MontFp!("752258130");
        initial_ntt[57] = MontFp!("97056798");
        initial_ntt[58] = MontFp!("87706193");
        initial_ntt[59] = MontFp!("1894023776");
        initial_ntt[60] = MontFp!("1292875010");
        initial_ntt[61] = MontFp!("929142525");
        initial_ntt[62] = MontFp!("1478368962");
        initial_ntt[63] = MontFp!("1981103952");
        initial_ntt[64] = MontFp!("1045062179");
        initial_ntt[65] = MontFp!("1425840265");
        initial_ntt[66] = MontFp!("2027265");
        initial_ntt[67] = MontFp!("1194516240");
        initial_ntt[68] = MontFp!("1379980743");
        initial_ntt[69] = MontFp!("1206462311");
        initial_ntt[70] = MontFp!("1498258989");
        initial_ntt[71] = MontFp!("1811570669");

        let mut expected = vec![Fq::ZERO; D];
        expected[0] = MontFp!("1065674974");
        expected[1] = MontFp!("1170569399");
        expected[2] = MontFp!("170751506");
        expected[3] = MontFp!("265022980");
        expected[4] = MontFp!("1945207175");
        expected[5] = MontFp!("458345263");
        expected[6] = MontFp!("2011655826");
        expected[7] = MontFp!("1046550861");
        expected[8] = MontFp!("264795716");
        expected[9] = MontFp!("1804913559");
        expected[10] = MontFp!("843380477");
        expected[11] = MontFp!("1398172716");
        expected[12] = MontFp!("851789181");
        expected[13] = MontFp!("1613109865");
        expected[14] = MontFp!("702106862");
        expected[15] = MontFp!("341684672");
        expected[16] = MontFp!("1577614606");
        expected[17] = MontFp!("307655228");
        expected[18] = MontFp!("1421181641");
        expected[19] = MontFp!("742137641");
        expected[20] = MontFp!("923616603");
        expected[21] = MontFp!("934523206");
        expected[22] = MontFp!("1207266670");
        expected[23] = MontFp!("487352988");
        expected[24] = MontFp!("958533374");
        expected[25] = MontFp!("997570189");
        expected[26] = MontFp!("746375437");
        expected[27] = MontFp!("449713270");
        expected[28] = MontFp!("1293462949");
        expected[29] = MontFp!("1967479755");
        expected[30] = MontFp!("1128550923");
        expected[31] = MontFp!("78875160");
        expected[32] = MontFp!("557134787");
        expected[33] = MontFp!("489984819");
        expected[34] = MontFp!("1473627119");
        expected[35] = MontFp!("1357428011");
        expected[36] = MontFp!("1445264686");
        expected[37] = MontFp!("380607359");
        expected[38] = MontFp!("1492417418");
        expected[39] = MontFp!("318205607");
        expected[40] = MontFp!("1329808119");
        expected[41] = MontFp!("619981352");
        expected[42] = MontFp!("1642500830");
        expected[43] = MontFp!("1919489665");
        expected[44] = MontFp!("957117942");
        expected[45] = MontFp!("1002845364");
        expected[46] = MontFp!("1844882309");
        expected[47] = MontFp!("153270753");
        expected[48] = MontFp!("1840769424");
        expected[49] = MontFp!("1401353601");
        expected[50] = MontFp!("1250603351");
        expected[51] = MontFp!("20939975");
        expected[52] = MontFp!("424643571");
        expected[53] = MontFp!("1025939175");
        expected[54] = MontFp!("1129236551");
        expected[55] = MontFp!("494827957");
        expected[56] = MontFp!("559340745");
        expected[57] = MontFp!("1220615690");
        expected[58] = MontFp!("96579813");
        expected[59] = MontFp!("1879163772");
        expected[60] = MontFp!("1730330419");
        expected[61] = MontFp!("1082601059");
        expected[62] = MontFp!("949990547");
        expected[63] = MontFp!("129594047");
        expected[64] = MontFp!("940074644");
        expected[65] = MontFp!("1825399223");
        expected[66] = MontFp!("476133872");
        expected[67] = MontFp!("1734778779");
        expected[68] = MontFp!("1594364605");
        expected[69] = MontFp!("1725670109");
        expected[70] = MontFp!("581029317");
        expected[71] = MontFp!("1343349559");

        homogenize_fq9(&mut initial_ntt);
        babybear_icrt_in_place(&mut initial_ntt);
        assert_eq!(initial_ntt, expected);
    }

    fn test_crt_icrt() {
        let mut rng = thread_rng();
        let coefficients: Vec<Fq> = (0..72).map(|_| Fq::rand(&mut rng)).collect();

        let mut input = coefficients.clone();

        serial_babybear_crt_in_place(&mut input);
        serial_babybear_icrt_in_place(&mut input);

        assert_eq!(coefficients, input);
    }
    fn test_crt_icrt_not_in_place() {
        let mut rng = thread_rng();
        let coefficients: Vec<Fq> = (0..72).map(|_| Fq::rand(&mut rng)).collect();

        let input = coefficients.clone();

        babybear_icrt(babybear_crt(input.clone()));

        assert_eq!(coefficients, input);
    }

    #[test]
    fn test_crt_icrt_not_in_place_1000000_times() {
        for _ in 0..1000000 {
            test_crt_icrt_not_in_place()
        }
    }

    #[test]
    fn test_crt_icrt_1000000_times() {
        for _ in 0..1000000 {
            test_crt_icrt()
        }
    }
}
