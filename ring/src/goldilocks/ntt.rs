use ark_ff::{Fp3Config, MontFp};

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

    serial_goldilock_crt_in_place(&mut coefficients)
}

fn serial_goldilock_crt_in_place(coefficients: &mut [Fq]) -> Vec<Fq3> {
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
                PRIMITIVE_TWENTYFORTH_ROOT_OF_UNITY_TO_5 * coeff_d_div_8_plus_i;

            coefficients[3 * D / 4 + i] = coeff_i + sigma_coeff_d_div_8_plus_i;
            coefficients[7 * D / 8 + i] = coeff_i - sigma_coeff_d_div_8_plus_i;
        }
    }

    vec![
        Fq3::new(coefficients[0], coefficients[1], coefficients[2]),
        nonresidue_to_13_to_nonresidue(coefficients[3], coefficients[4], coefficients[5]),
        nonresidue_to_7_to_nonresidue(coefficients[6], coefficients[7], coefficients[8]),
        nonresidue_to_19_to_nonresidue(coefficients[9], coefficients[10], coefficients[11]),
        nonresidue_to_5_to_nonresidue(coefficients[12], coefficients[13], coefficients[14]),
        nonresidue_to_17_to_nonresidue(coefficients[15], coefficients[16], coefficients[17]),
        nonresidue_to_11_to_nonresidue(coefficients[18], coefficients[19], coefficients[20]),
        nonresidue_to_23_to_nonresidue(coefficients[21], coefficients[22], coefficients[23]),
    ]
}

// Different automorphisms with the target Fp(NONRESIDUE) and their inverses.
#[inline(always)]
const fn nonresidue_to_13_to_nonresidue(c0: Fq, c1: Fq, c2: Fq) -> Fq3 {
    Fq3::new(c0, -c1, c2)
}

#[inline(always)]
const fn nonresidue_to_nonresidue_to_13(c: Fq3) -> (Fq, Fq, Fq) {
    (c.c0, -c.c1, c.c2)
}

#[inline(always)]
const fn nonresidue_to_7_to_nonresidue(c0: Fq, c1: Fq, c2: Fq) -> Fq3 {
    Fq3::new(c0, ROOTS_OF_UNITY_24[2] * c1, ROOTS_OF_UNITY_24[4] * c2)
}

#[inline(always)]
const fn nonresidue_to_nonresidue_to_7(c: Fq3) -> (Fq, Fq, Fq) {
    (c.c0, ROOTS_OF_UNITY_24[22] * c.c1, ROOTS_OF_UNITY_24[20] * c.c2)
}

#[inline(always)]
const fn nonresidue_to_19_to_nonresidue(c0: Fq, c1: Fq, c2: Fq) -> Fq3 {
    Fq3::new(c0, ROOTS_OF_UNITY_24[6] * c1, ROOTS_OF_UNITY_24[12] * c2)
}

#[inline(always)]
const fn nonresidue_to_nonresidue_to_19(c: Fq3) -> (Fq, Fq, Fq) {
    (c.c0, ROOTS_OF_UNITY_24[18] * c.c1, ROOTS_OF_UNITY_24[12] * c.c2)
}

#[inline(always)]
const fn nonresidue_to_5_to_nonresidue(c0: Fq, c1: Fq, c2: Fq) -> Fq3 {
    Fq3::new(c0, ROOTS_OF_UNITY_24[3] * c2, ROOTS_OF_UNITY_24[1] * c1)
}

#[inline(always)]
const fn nonresidue_to_nonresidue_to_5(c: Fq3) -> (Fq, Fq, Fq) {
    (c.c0, ROOTS_OF_UNITY_24[23] * c.c2, ROOTS_OF_UNITY_24[21] * c.c1)
}

#[inline(always)]
const fn nonresidue_to_17_to_nonresidue(c0: Fq, c1: Fq, c2: Fq) -> Fq3 {
    Fq3::new(c0, ROOTS_OF_UNITY_24[11] * c2, ROOTS_OF_UNITY_24[5] * c1)
}

#[inline(always)]
const fn nonresidue_to_nonresidue_to_17(c: Fq3) -> (Fq, Fq, Fq) {
    (c.c0, ROOTS_OF_UNITY_24[19] * c.c2, ROOTS_OF_UNITY_24[13] * c.c1)
}

#[inline(always)]
const fn nonresidue_to_11_to_nonresidue(c0: Fq, c1: Fq, c2: Fq) -> Fq3 {
    Fq3::new(c0, ROOTS_OF_UNITY_24[7] * c2, ROOTS_OF_UNITY_24[3] * c1)
}

#[inline(always)]
const fn nonresidue_to_nonresidue_to_11(c: Fq3) -> (Fq, Fq, Fq) {
    (c.c0, ROOTS_OF_UNITY_24[21] * c.c2, ROOTS_OF_UNITY_24[17] * c.c1)
}

#[inline(always)]
const fn nonresidue_to_23_to_nonresidue(c0: Fq, c1: Fq, c2: Fq) -> Fq3 {
    Fq3::new(c0, ROOTS_OF_UNITY_24[15] * c2, ROOTS_OF_UNITY_24[7] * c1)
}

#[inline(always)]
const fn nonresidue_to_nonresidue_to_23(c: Fq3) -> (Fq, Fq, Fq) {
    (c.c0, ROOTS_OF_UNITY_24[17] * c.c2, ROOTS_OF_UNITY_24[9] * c.c1)
}

#[cfg(test)]
mod tests{
    use ark_ff::{Field, UniformRand};
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
            for j in (i+1)..24 {
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
                let x_prime = $isomorphism(x);
                assert_eq!($inverse(x_prime.0, x_prime.1, x_prime.2), x);
            })+
        };
    }

    #[test]
    fn test_isomorphisms_inverses() {
        test_inverses!{
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
                let x_prime = $inverse(Fq::ZERO, Fq::ONE, Fq::ZERO);
                let x_prime_2 = $inverse(Fq::ZERO, Fq::ZERO, Fq::ONE);

                assert_eq!(x_prime * x_prime, x_prime_2);
            })+
        };
    }

    #[test]
    fn test_squares() {
        test_x_square!{
            nonresidue_to_13_to_nonresidue,
            nonresidue_to_7_to_nonresidue,
            nonresidue_to_19_to_nonresidue,
            nonresidue_to_5_to_nonresidue,
            nonresidue_to_17_to_nonresidue,
            nonresidue_to_11_to_nonresidue,
            nonresidue_to_23_to_nonresidue
        };
    }
}