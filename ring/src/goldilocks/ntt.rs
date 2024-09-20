#![allow(dead_code)]
use ark_ff::{Fp3Config, MontFp};
use ark_std::Zero;

use super::{Fq, Fq3, Goldilocks3Config};

/// The degree of the cyclotomic polynomial.
const D: usize = 24;

/// The number of splits of the cyclotomic ring in the CRT-form.
const N: usize = 8;

/// zeta
const PRIMITIVE_SIXTH_ROOT_OF_UNITY: Fq = MontFp!("18446744065119617026");
/// zeta ^ 5
const PRIMITIVE_SIXTH_ROOT_OF_UNITY_TO_FIVE: Fq = MontFp!("4294967296");
// sigma
const PRIMITIVE_TWELFTH_ROOT_OF_UNITY: Fq = MontFp!("281474976645120");
const PRIMITIVE_TWELFTH_ROOT_OF_UNITY_TO_FIVE: Fq = MontFp!("65536");

const PRIMITIVE_TWENTYFORTH_ROOT_OF_UNITY: Fq = Goldilocks3Config::NONRESIDUE;
const PRIMITIVE_TWENTYFORTH_ROOT_OF_UNITY_TO_7: Fq = MontFp!("72057594021150720");
const PRIMITIVE_TWENTYFORTH_ROOT_OF_UNITY_TO_5: Fq = MontFp!("256");
const PRIMITIVE_TWENTYFORTH_ROOT_OF_UNITY_TO_11: Fq = MontFp!("72057594037927936");

const ROOTS_OF_UNITY_24: &'static [Fq] = &[
    MontFp!("1"),
    MontFp!("1099511627776"),
    MontFp!("281474976645120"),
    MontFp!("18446744069397807105"),
    MontFp!("18446744065119617026"),
    MontFp!("256"),
    MontFp!("281474976710656"),
    MontFp!("72057594021150720"),
    MontFp!("18446744065119617025"),
    MontFp!("18446742969902956801"),
    MontFp!("65536"),
    MontFp!("72057594037927936"),
    MontFp!("18446744069414584320"),
    MontFp!("18446742969902956545"),
    MontFp!("18446462594437939201"),
    MontFp!("16777216"),
    MontFp!("4294967295"),
    MontFp!("18446744069414584065"),
    MontFp!("18446462594437873665"),
    MontFp!("18374686475393433601"),
    MontFp!("4294967296"),
    MontFp!("1099511627520"),
    MontFp!("18446744069414518785"),
    MontFp!("18374686475376656385"),
];

fn serial_goldilock_crt(coefficients: &[Fq]) -> Vec<Fq3> {
    let mut coefficients = coefficients.to_vec();

    serial_goldilock_crt_in_place(&mut coefficients);

    let mut result = vec![Fq3::zero(); N];

    for i in 0..N {
        result[i] = Fq3::new(
            coefficients[i * 3],
            coefficients[i * 3 + 1],
            coefficients[i * 3 + 2],
        );
    }

    result
}

fn serial_goldilock_crt_in_place(coefficients: &mut [Fq]) {
    assert!(coefficients.len() == D);

    // Compute f mod X^12-zeta and f mod X^12 - zeta^5
    for i in 0..(D / 2) {
        let (coeff_i, coeff_d_div_2_plus_i) = (coefficients[i], coefficients[D / 2 + i]);
        let zeta_coeff_d_div_2_plus_i = PRIMITIVE_SIXTH_ROOT_OF_UNITY * coeff_d_div_2_plus_i;

        coefficients[i] = coeff_i + zeta_coeff_d_div_2_plus_i;
        coefficients[i + D / 2] = coeff_i + coeff_d_div_2_plus_i - zeta_coeff_d_div_2_plus_i;
    }
    // After this step we get f_1 in the first half coefficients and f_2 in the second half.

    // Compute f_1 mod X^6-\sigma, f_1 mod X^6-\sigma^7, f_2 mod X^6 -\sigma^5, f_2 mod X^6 -\sigma^11
    for i in 0..(D / 4) {
        // f_1
        {
            let (coeff_i, coeff_d_div_4_plus_i) = (coefficients[i], coefficients[D / 4 + i]);
            let sigma_coeff_d_div_4_plus_i = PRIMITIVE_TWELFTH_ROOT_OF_UNITY * coeff_d_div_4_plus_i;

            coefficients[i] = coeff_i + sigma_coeff_d_div_4_plus_i;
            coefficients[i + D / 4] = coeff_i - sigma_coeff_d_div_4_plus_i;
        }

        // f_2
        {
            let (coeff_i, coeff_d_div_4_plus_i) =
                (coefficients[D / 2 + i], coefficients[3 * D / 4 + i]);
            let sigma_coeff_d_div_4_plus_i =
                PRIMITIVE_TWELFTH_ROOT_OF_UNITY_TO_FIVE * coeff_d_div_4_plus_i;

            coefficients[D / 2 + i] = coeff_i + sigma_coeff_d_div_4_plus_i;
            coefficients[3 * D / 4 + i] = coeff_i - sigma_coeff_d_div_4_plus_i;
        }
    }
    // After this step we have f_1, f_2, f_3, f_4.

    // Compute f_1 mod X^3-NONRESIDUE, f_1 mod X^3-NONRESIDUE^13, f_2 mod X^3-NONRESIDUE^7, f_2 mod X^3-NONRESIDUE^19
    // Compute f_3 mod X^3-NONRESIDUE^5, f_3 mod X^3-NONRESIDUE^17, f_4 mod X^3-NONRESIDUE^11, f_4 mod X^3-NONRESIDUE^23
    for i in 0..(D / 8) {
        // f_1
        {
            let (coeff_i, coeff_d_div_8_plus_i) = (coefficients[i], coefficients[D / 8 + i]);
            let nonresidue_coeff_d_div_8_plus_i =
                PRIMITIVE_TWENTYFORTH_ROOT_OF_UNITY * coeff_d_div_8_plus_i;

            coefficients[i] = coeff_i + nonresidue_coeff_d_div_8_plus_i;
            coefficients[i + D / 8] = coeff_i - nonresidue_coeff_d_div_8_plus_i;
        }

        // f_2
        {
            let (coeff_i, coeff_d_div_8_plus_i) =
                (coefficients[D / 4 + i], coefficients[3 * D / 8 + i]);
            let sigma_coeff_d_div_8_plus_i =
                PRIMITIVE_TWENTYFORTH_ROOT_OF_UNITY_TO_7 * coeff_d_div_8_plus_i;

            coefficients[D / 4 + i] = coeff_i + sigma_coeff_d_div_8_plus_i;
            coefficients[3 * D / 8 + i] = coeff_i - sigma_coeff_d_div_8_plus_i;
        }

        // f_3
        {
            let (coeff_i, coeff_d_div_8_plus_i) =
                (coefficients[D / 2 + i], coefficients[5 * D / 8 + i]);
            let sigma_coeff_d_div_8_plus_i =
                PRIMITIVE_TWENTYFORTH_ROOT_OF_UNITY_TO_5 * coeff_d_div_8_plus_i;

            coefficients[D / 2 + i] = coeff_i + sigma_coeff_d_div_8_plus_i;
            coefficients[5 * D / 8 + i] = coeff_i - sigma_coeff_d_div_8_plus_i;
        }

        // f_4
        {
            let (coeff_i, coeff_d_div_8_plus_i) =
                (coefficients[3 * D / 4 + i], coefficients[7 * D / 8 + i]);
            let sigma_coeff_d_div_8_plus_i =
                PRIMITIVE_TWENTYFORTH_ROOT_OF_UNITY_TO_11 * coeff_d_div_8_plus_i;

            coefficients[3 * D / 4 + i] = coeff_i + sigma_coeff_d_div_8_plus_i;
            coefficients[7 * D / 8 + i] = coeff_i - sigma_coeff_d_div_8_plus_i;
        }
    }

    normalize_fq3(coefficients);
}

fn serial_goldilock_icrt(evaluations: &[Fq3]) {
    let mut evaluations_field: Vec<Fq> = vec![Fq::zero(); N * 3];

    assert_eq!(evaluations.len(), N);

    for i in 0..N {
        evaluations_field[i * 3] = evaluations[i].c0;
        evaluations_field[i * 3 + 1] = evaluations[i].c1;
        evaluations_field[i * 3 + 2] = evaluations[i].c2;
    }

    serial_goldilock_icrt_in_place(&mut evaluations_field);
}

fn serial_goldilock_icrt_in_place(evaluations: &mut [Fq]) {
    assert_eq!(evaluations.len(), D);

    denormalize_fq3(evaluations);

    // After this step we have f_1, f_2, f_3, f_4.

    // Compute f_1 mod X^3-NONRESIDUE, f_1 mod X^3-NONRESIDUE^13, f_2 mod X^3-NONRESIDUE^7, f_2 mod X^3-NONRESIDUE^19
    // Compute f_3 mod X^3-NONRESIDUE^5, f_3 mod X^3-NONRESIDUE^17, f_4 mod X^3-NONRESIDUE^11, f_4 mod X^3-NONRESIDUE^23
    // for i in 0..(D / 8) {
    //     // f_1
    //     {
    //         let (coeff_i, coeff_d_div_8_plus_i) = (coefficients[i], coefficients[D / 8 + i]);
    //         let nonresidue_coeff_d_div_8_plus_i =
    //             PRIMITIVE_TWENTYFORTH_ROOT_OF_UNITY * coeff_d_div_8_plus_i;

    //         coefficients[i] = coeff_i + nonresidue_coeff_d_div_8_plus_i;
    //         coefficients[i + D / 8] = coeff_i - nonresidue_coeff_d_div_8_plus_i;
    //     }

    //     // f_2
    //     {
    //         let (coeff_i, coeff_d_div_8_plus_i) =
    //             (coefficients[D / 4 + i], coefficients[3 * D / 8 + i]);
    //         let sigma_coeff_d_div_8_plus_i =
    //             PRIMITIVE_TWENTYFORTH_ROOT_OF_UNITY_TO_7 * coeff_d_div_8_plus_i;

    //         coefficients[D / 4 + i] = coeff_i + sigma_coeff_d_div_8_plus_i;
    //         coefficients[3 * D / 8 + i] = coeff_i - sigma_coeff_d_div_8_plus_i;
    //     }

    //     // f_3
    //     {
    //         let (coeff_i, coeff_d_div_8_plus_i) =
    //             (coefficients[D / 2 + i], coefficients[5 * D / 8 + i]);
    //         let sigma_coeff_d_div_8_plus_i =
    //             PRIMITIVE_TWENTYFORTH_ROOT_OF_UNITY_TO_5 * coeff_d_div_8_plus_i;

    //         coefficients[D / 2 + i] = coeff_i + sigma_coeff_d_div_8_plus_i;
    //         coefficients[5 * D / 8 + i] = coeff_i - sigma_coeff_d_div_8_plus_i;
    //     }

    //     // f_4
    //     {
    //         let (coeff_i, coeff_d_div_8_plus_i) =
    //             (coefficients[3 * D / 4 + i], coefficients[7 * D / 8 + i]);
    //         let sigma_coeff_d_div_8_plus_i =
    //             PRIMITIVE_TWENTYFORTH_ROOT_OF_UNITY_TO_5 * coeff_d_div_8_plus_i;

    //         coefficients[3 * D / 4 + i] = coeff_i + sigma_coeff_d_div_8_plus_i;
    //         coefficients[7 * D / 8 + i] = coeff_i - sigma_coeff_d_div_8_plus_i;
    //     }
    // }
}

/// Converts each triple `c[3 * i]`, `c[3 * i + 1]`, `c[3 * i + 2]`
/// into coefficients of an element from Fq[X]/(X^3-NONRESIDUE).
#[inline(always)]
fn normalize_fq3(c: &mut [Fq]) {
    nonresidue_to_13_to_nonresidue(&mut c[3..6]);
    nonresidue_to_7_to_nonresidue(&mut c[6..9]);
    nonresidue_to_19_to_nonresidue(&mut c[9..12]);
    nonresidue_to_5_to_nonresidue(&mut c[12..15]);
    nonresidue_to_17_to_nonresidue(&mut c[15..18]);
    nonresidue_to_11_to_nonresidue(&mut c[18..21]);
    nonresidue_to_23_to_nonresidue(&mut c[21..24]);
}

/// The inverse of the above.
#[inline(always)]
fn denormalize_fq3(c: &mut [Fq]) {
    nonresidue_to_nonresidue_to_13(&mut c[3..6]);
    nonresidue_to_nonresidue_to_7(&mut c[6..9]);
    nonresidue_to_nonresidue_to_19(&mut c[9..12]);
    nonresidue_to_nonresidue_to_5(&mut c[12..15]);
    nonresidue_to_nonresidue_to_17(&mut c[15..18]);
    nonresidue_to_nonresidue_to_11(&mut c[18..21]);
    nonresidue_to_nonresidue_to_23(&mut c[21..24]);
}

// Different automorphisms with the target Fp(NONRESIDUE) and their inverses.
#[inline(always)]
fn nonresidue_to_13_to_nonresidue(c: &mut [Fq]) {
    c[1] = -c[1];
}

#[inline(always)]
fn nonresidue_to_nonresidue_to_13(c: &mut [Fq]) {
    c[1] = -c[1];
}

#[inline(always)]
fn nonresidue_to_7_to_nonresidue(c: &mut [Fq]) {
    c[1] *= ROOTS_OF_UNITY_24[2];
    c[2] *= ROOTS_OF_UNITY_24[4];
}

#[inline(always)]
fn nonresidue_to_nonresidue_to_7(c: &mut [Fq]) {
    c[1] *= ROOTS_OF_UNITY_24[22];
    c[2] *= ROOTS_OF_UNITY_24[20];
}

#[inline(always)]
fn nonresidue_to_19_to_nonresidue(c: &mut [Fq]) {
    c[1] *= ROOTS_OF_UNITY_24[6];
    c[2] *= ROOTS_OF_UNITY_24[12];
}

#[inline(always)]
fn nonresidue_to_nonresidue_to_19(c: &mut [Fq]) {
    c[1] *= ROOTS_OF_UNITY_24[18];
    c[2] *= ROOTS_OF_UNITY_24[12];
}

#[inline(always)]
fn nonresidue_to_5_to_nonresidue(c: &mut [Fq]) {
    let c1 = c[1];
    c[1] = c[2] * ROOTS_OF_UNITY_24[3];
    c[2] = c1 * ROOTS_OF_UNITY_24[1];
}

#[inline(always)]
fn nonresidue_to_nonresidue_to_5(c: &mut [Fq]) {
    let c1 = c[1];
    c[1] = c[2] * ROOTS_OF_UNITY_24[23];
    c[2] = c1 * ROOTS_OF_UNITY_24[21];
}

#[inline(always)]
fn nonresidue_to_17_to_nonresidue(c: &mut [Fq]) {
    let c1 = c[1];
    c[1] = c[2] * ROOTS_OF_UNITY_24[11];
    c[2] = c1 * ROOTS_OF_UNITY_24[5];
}

#[inline(always)]
fn nonresidue_to_nonresidue_to_17(c: &mut [Fq]) {
    let c1 = c[1];
    c[1] = c[2] * ROOTS_OF_UNITY_24[19];
    c[2] = c1 * ROOTS_OF_UNITY_24[13];
}

#[inline(always)]
fn nonresidue_to_11_to_nonresidue(c: &mut [Fq]) {
    let c1 = c[1];
    c[1] = c[2] * ROOTS_OF_UNITY_24[7];
    c[2] = c1 * ROOTS_OF_UNITY_24[3];
}

#[inline(always)]
fn nonresidue_to_nonresidue_to_11(c: &mut [Fq]) {
    let c1 = c[1];
    c[1] = c[2] * ROOTS_OF_UNITY_24[21];
    c[2] = c1 * ROOTS_OF_UNITY_24[17];
}

#[inline(always)]
fn nonresidue_to_23_to_nonresidue(c: &mut [Fq]) {
    let c1 = c[1];
    c[1] = c[2] * ROOTS_OF_UNITY_24[15];
    c[2] = c1 * ROOTS_OF_UNITY_24[7];
}

#[inline(always)]
fn nonresidue_to_nonresidue_to_23(c: &mut [Fq]) {
    let c1 = c[1];
    c[1] = c[2] * ROOTS_OF_UNITY_24[17];
    c[2] = c1 * ROOTS_OF_UNITY_24[9];
}

#[cfg(test)]
mod tests {
    use ark_ff::{Field, UniformRand};
    use ark_std::One;
    use rand::thread_rng;

    use crate::goldilocks::Fq3;

    use super::*;

    #[test]
    fn test_roots_of_unity() {
        // Check if they all are roots of unity of degree 24.
        for root in ROOTS_OF_UNITY_24 {
            assert_eq!(root.pow([24]), Fq::ONE);
        }

        // Check if they are pairwise distinct.
        for i in 0..24 {
            for j in (i + 1)..24 {
                assert_ne!(ROOTS_OF_UNITY_24[i], ROOTS_OF_UNITY_24[j]);
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
                let x = Fq3::rand(&mut rng);
                let mut x_prime: Vec<Fq> = vec![x.c0, x.c1, x.c2];
                $isomorphism(&mut x_prime);
                $inverse(&mut x_prime);
                assert_eq!(Fq3::new(x_prime[0], x_prime[1], x_prime[2]), x);
            })+
        };
    }

    #[test]
    fn test_isomorphisms_inverses() {
        test_inverses! {
            nonresidue_to_nonresidue_to_13, nonresidue_to_13_to_nonresidue,
            nonresidue_to_nonresidue_to_7, nonresidue_to_7_to_nonresidue,
            nonresidue_to_nonresidue_to_19, nonresidue_to_19_to_nonresidue,
            nonresidue_to_nonresidue_to_5, nonresidue_to_5_to_nonresidue,
            nonresidue_to_nonresidue_to_17, nonresidue_to_17_to_nonresidue,
            nonresidue_to_nonresidue_to_11, nonresidue_to_11_to_nonresidue,
            nonresidue_to_nonresidue_to_23, nonresidue_to_23_to_nonresidue
        };
    }

    macro_rules! test_x_square {
        ($($inverse:expr),+) => {
            $({
                let mut x_prime: Vec<Fq> = vec![Fq::ZERO, Fq::ONE, Fq::ZERO];
                let mut x_prime_2: Vec<Fq> = vec![Fq::ZERO, Fq::ZERO, Fq::ONE];
                $inverse(&mut x_prime);
                $inverse(&mut x_prime_2);

                let x_prime = Fq3::new(x_prime[0], x_prime[1], x_prime[2]);
                let x_prime_2 = Fq3::new(x_prime_2[0], x_prime_2[1], x_prime_2[2]);

                assert_eq!(x_prime * x_prime, x_prime_2);
            })+
        };
    }

    #[test]
    fn test_squares() {
        test_x_square! {
            nonresidue_to_13_to_nonresidue,
            nonresidue_to_7_to_nonresidue,
            nonresidue_to_19_to_nonresidue,
            nonresidue_to_5_to_nonresidue,
            nonresidue_to_17_to_nonresidue,
            nonresidue_to_11_to_nonresidue,
            nonresidue_to_23_to_nonresidue
        };
    }

    #[test]
    fn test_normalize_denormalize() {
        let mut rng = thread_rng();
        let mut x: Vec<Fq> = (0..24).map(|_| Fq::rand(&mut rng)).collect();

        let expected = x.clone();

        normalize_fq3(&mut x);
        denormalize_fq3(&mut x);

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
            Fq::zero(),
            Fq::zero(),
            Fq::zero(),
            Fq::zero(),
            Fq::zero(),
            Fq::zero(),
            Fq::zero(),
            Fq::one(),
        ];

        let expected: Vec<Fq> = vec![
            MontFp!("3841"),
            MontFp!("2"),
            MontFp!("72057594021150723"),
            MontFp!("18446744069414580482"),
            MontFp!("2"),
            MontFp!("18374686475393433604"),
            MontFp!("1080863910568919041"),
            MontFp!("2"),
            MontFp!("1099511627779"),
            MontFp!("17365880158845665282"),
            MontFp!("2"),
            MontFp!("18446742969902956548"),
            MontFp!("16492674416641"),
            MontFp!("2"),
            MontFp!("72057594037927939"),
            MontFp!("18446727576740167682"),
            MontFp!("2"),
            MontFp!("18374686475376656388"),
            MontFp!("1080863910317260801"),
            MontFp!("2"),
            MontFp!("259"),
            MontFp!("17365880159097323522"),
            MontFp!("2"),
            MontFp!("18446744069414584068"),
        ];

        serial_goldilock_crt_in_place(&mut test_poly);

        denormalize_fq3(&mut test_poly);

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
            MontFp!("11977680547482164101"),
            MontFp!("543543"),
            MontFp!("488514175862046709"),
            MontFp!("11976965864924109701"),
            MontFp!("543543"),
            MontFp!("17958229893552537618"),
            MontFp!("11441394850670851783"),
            MontFp!("543543"),
            MontFp!("10160120756981332284"),
            MontFp!("1497446875752052425"),
            MontFp!("543543"),
            MontFp!("8286623312433252043"),
            MontFp!("50172301757990"),
            MontFp!("543543"),
            MontFp!("60243283203"),
            MontFp!("50172301591590"),
            MontFp!("543543"),
            MontFp!("18446744009171301124"),
            MontFp!("4971923820610324773"),
            MontFp!("543543"),
            MontFp!("10164068860789127484"),
            MontFp!("13474719904200919336"),
            MontFp!("543543"),
            MontFp!("8282675208625456843"),
        ];

        serial_goldilock_crt_in_place(&mut test_poly);

        denormalize_fq3(&mut test_poly);

        assert_eq!(test_poly, expected);
    }
}
