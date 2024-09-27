//!
//! A CRT implementation for the ring Fq[X]/(X^16 + 1), where q is the Starknet prime.
//!
use ark_ff::MontFp;

use super::Fq;

/// The degree of the cyclotomic polynomial.
pub(super) const D: usize = 16;

/// The number of splits of the cyclotomic ring in the CRT-form.
pub(super) const N: usize = 16;

/// All 32 roots of unity of degree 32.
const ROOTS_OF_UNITY_32: &[Fq] = &[
    MontFp!("1"),
    MontFp!("3409443867035641044245057348756544640549407421541289951053907001322227935403"),
    MontFp!("2679026602897868112349604024891625875968950767352485125058791696935099163961"),
    MontFp!("1247662575873211570659477042654408208661347962178767032684320489684194658007"),
    MontFp!("2804690217475462062143361339624939640984649667966511418446363596075299761851"),
    MontFp!("936046279675834030954630398480954073813587618256291413305593916375183853324"),
    MontFp!("20018235416501333678500397115203031959003915510248567384130230248350883666"),
    MontFp!("1517735079997376536627264160940840285257771804302276892720884016887753722993"),
    MontFp!("2779265777490745486907812647164038750808419677710890077123925563301162072035"),
    MontFp!("1500963833878796838526042408495665868007267193736142409970969381386334285397"),
    MontFp!("2625973900227797461931995796145047439335352296623810184744144340654095677649"),
    MontFp!("1839320071806632312009665796873469510898389465043833405723636246914115759733"),
    MontFp!("1683934744455360247964225572400674938023097276580318475317001647084465960397"),
    MontFp!("2950642793828688006577837425001810606298859709516175614105408635229182523157"),
    MontFp!("1010767318980875148721624974497647105601682168692758069501464334811644188044"),
    MontFp!("1947442199807235156009795270208565947955700433016980456436986602271286759616"),
    MontFp!("3618502788666131213697322783095070105623107215331596699973092056135872020480"),
    MontFp!("209058921630490169452265434338525465073699793790306748919185054813644085078"),
    MontFp!("939476185768263101347718758203444229654156447979111574914300359200772856520"),
    MontFp!("2370840212792919643037845740440661896961759253152829667288771566451677362474"),
    MontFp!("813812571190669151553961443470130464638457547365085281526728460060572258630"),
    MontFp!("2682456508990297182742692384614116031809519597075305286667498139760688167157"),
    MontFp!("3598484553249629880018822385979867073664103299821348132588961825887521136815"),
    MontFp!("2100767708668754677070058622154229820365335411029319807252208039248118297488"),
    MontFp!("839237011175385726789510135931031354814687537620706622849166492834709948446"),
    MontFp!("2117538954787334375171280374599404237615840021595454290002122674749537735084"),
    MontFp!("992528888438333751765326986950022666287754918707786515228947715481776342832"),
    MontFp!("1779182716859498901687656986221600594724717750287763294249455809221756260748"),
    MontFp!("1934568044210770965733097210694395167600009938751278224656090409051406060084"),
    MontFp!("667859994837443207119485358093259499324247505815421085867683420906689497324"),
    MontFp!("2607735469685256064975697808597423000021425046638838630471627721324227832437"),
    MontFp!("1671060588858896057687527512886504157667406782314616243536105453864585260865"),
];

const SIXTEEN_INV: Fq =
    MontFp!("3392346364374498012841240109151628224021663014373371906224773802627380019201");
/// 1 / 16 * ROOTS_OF_UNITY_32[24]
const SIXTEEN_INV_TIMES_ROOT_OF_UNITY_32_24: Fq =
    MontFp!("504765161781728009636509731382573222878806373017743751424709412819153374338");

/// Given `coefficients` of a polynoimial `f mod X^16 + 1`
/// returns its evaluations on the sequence:
/// `ROOTS_OF_UNITY_32[1], ROOTS_OF_UNITY_32[17], ROOTS_OF_UNITY_32[9], ROOTS_OF_UNITY_32[25], ROOTS_OF_UNITY_32[5], ROOTS_OF_UNITY_32[21],
/// ROOTS_OF_UNITY_32[13], ROOTS_OF_UNITY_32[29], ROOTS_OF_UNITY_32[3], ROOTS_OF_UNITY_32[19], ROOTS_OF_UNITY_32[11], ROOTS_OF_UNITY_32[27],
/// ROOTS_OF_UNITY_32[7], ROOTS_OF_UNITY_32[23], ROOTS_OF_UNITY_32[15], ROOTS_OF_UNITY_32[31]`. In this exact order.
///
/// # Panics
///
/// Panics if `coefficients.len() != 16`.
#[inline(always)]
pub fn stark_prime_crt(coefficients: Vec<Fq>) -> Vec<Fq> {
    // TODO: once we have parallelization set up
    //       this should be either parallel or serial.
    serial_stark_prime_crt(coefficients)
}

/// Same as `stark_prime_crt` but performs the CRT in place.
///
/// # Panics
///
/// Panics if `coefficients.len() != 16`.
#[inline(always)]
pub fn stark_prime_crt_in_place(coefficients: &mut [Fq]) {
    // TODO: once we have parallelization set up
    //       this should be either parallel or serial.
    serial_stark_prime_crt_in_place(coefficients)
}

/// The inverse CRT. The order of evaluations as in the return of the direct CRT.
///
/// Returns the coefficients of the polynomial encoded by this CRT form.
///
/// # Panics
///
/// Panics if `evaluations.len() != 16`.
#[inline(always)]
pub fn stark_prime_icrt(evaluations: Vec<Fq>) -> Vec<Fq> {
    // TODO: once we have parallelization set up
    //       this should be either parallel or serial.
    serial_stark_prime_icrt(evaluations)
}

/// Same as `stark_prime_icrt` but performs the inverse CRT in place.
///
/// # Panics
///
/// Panics if `evaluations.len() != 16`.
#[inline(always)]
pub fn stark_prime_icrt_in_place(evaluations: &mut [Fq]) {
    // TODO: once we have parallelization set up
    //       this should be either parallel or serial.
    serial_stark_prime_icrt_in_place(evaluations)
}

fn serial_stark_prime_crt(mut coefficients: Vec<Fq>) -> Vec<Fq> {
    serial_stark_prime_crt_in_place(&mut coefficients);

    coefficients
}

fn serial_stark_prime_crt_in_place(coefficients: &mut [Fq]) {
    assert_eq!(coefficients.len(), D);

    for i in 0..(D / 2) {
        let (coeff_i, coeff_d_div_2_plus_i) = (coefficients[i], coefficients[D / 2 + i]);
        let zeta_coeff_d_div_2_plus_i = ROOTS_OF_UNITY_32[8] * coeff_d_div_2_plus_i;

        coefficients[i] = coeff_i + zeta_coeff_d_div_2_plus_i;
        coefficients[i + D / 2] = coeff_i - zeta_coeff_d_div_2_plus_i;
    }

    for i in 0..(D / 4) {
        // f_1
        {
            let (coeff_i, coeff_d_div_4_plus_i) = (coefficients[i], coefficients[D / 4 + i]);
            let sigma_coeff_d_div_4_plus_i = ROOTS_OF_UNITY_32[4] * coeff_d_div_4_plus_i;

            coefficients[i] = coeff_i + sigma_coeff_d_div_4_plus_i;
            coefficients[i + D / 4] = coeff_i - sigma_coeff_d_div_4_plus_i;
        }

        // f_2
        {
            let (coeff_i, coeff_d_div_4_plus_i) =
                (coefficients[D / 2 + i], coefficients[3 * D / 4 + i]);
            let sigma_coeff_d_div_4_plus_i = ROOTS_OF_UNITY_32[12] * coeff_d_div_4_plus_i;

            coefficients[D / 2 + i] = coeff_i + sigma_coeff_d_div_4_plus_i;
            coefficients[3 * D / 4 + i] = coeff_i - sigma_coeff_d_div_4_plus_i;
        }
    }

    for i in 0..(D / 8) {
        // f_1
        {
            let (coeff_i, coeff_d_div_8_plus_i) = (coefficients[i], coefficients[D / 8 + i]);
            let nonresidue_coeff_d_div_8_plus_i = ROOTS_OF_UNITY_32[2] * coeff_d_div_8_plus_i;

            coefficients[i] = coeff_i + nonresidue_coeff_d_div_8_plus_i;
            coefficients[i + D / 8] = coeff_i - nonresidue_coeff_d_div_8_plus_i;
        }

        // f_2
        {
            let (coeff_i, coeff_d_div_8_plus_i) =
                (coefficients[D / 4 + i], coefficients[3 * D / 8 + i]);
            let sigma_coeff_d_div_8_plus_i = ROOTS_OF_UNITY_32[10] * coeff_d_div_8_plus_i;

            coefficients[D / 4 + i] = coeff_i + sigma_coeff_d_div_8_plus_i;
            coefficients[3 * D / 8 + i] = coeff_i - sigma_coeff_d_div_8_plus_i;
        }

        // f_3
        {
            let (coeff_i, coeff_d_div_8_plus_i) =
                (coefficients[D / 2 + i], coefficients[5 * D / 8 + i]);
            let sigma_coeff_d_div_8_plus_i = ROOTS_OF_UNITY_32[6] * coeff_d_div_8_plus_i;

            coefficients[D / 2 + i] = coeff_i + sigma_coeff_d_div_8_plus_i;
            coefficients[5 * D / 8 + i] = coeff_i - sigma_coeff_d_div_8_plus_i;
        }

        // f_4
        {
            let (coeff_i, coeff_d_div_8_plus_i) =
                (coefficients[3 * D / 4 + i], coefficients[7 * D / 8 + i]);
            let sigma_coeff_d_div_8_plus_i = ROOTS_OF_UNITY_32[14] * coeff_d_div_8_plus_i;

            coefficients[3 * D / 4 + i] = coeff_i + sigma_coeff_d_div_8_plus_i;
            coefficients[7 * D / 8 + i] = coeff_i - sigma_coeff_d_div_8_plus_i;
        }
    }

    {
        let (a_0, mut a_1) = (coefficients[0], coefficients[1]);
        a_1 *= ROOTS_OF_UNITY_32[1];
        coefficients[0] += a_1;
        coefficients[1] = a_0 - a_1;

        let (a_0, mut a_1) = (coefficients[2], coefficients[3]);
        a_1 *= ROOTS_OF_UNITY_32[9];
        coefficients[2] += a_1;
        coefficients[3] = a_0 - a_1;

        let (a_0, mut a_1) = (coefficients[4], coefficients[5]);
        a_1 *= ROOTS_OF_UNITY_32[5];
        coefficients[4] += a_1;
        coefficients[5] = a_0 - a_1;

        let (a_0, mut a_1) = (coefficients[6], coefficients[7]);
        a_1 *= ROOTS_OF_UNITY_32[13];
        coefficients[6] += a_1;
        coefficients[7] = a_0 - a_1;

        let (a_0, mut a_1) = (coefficients[8], coefficients[9]);
        a_1 *= ROOTS_OF_UNITY_32[3];
        coefficients[8] += a_1;
        coefficients[9] = a_0 - a_1;

        let (a_0, mut a_1) = (coefficients[10], coefficients[11]);
        a_1 *= ROOTS_OF_UNITY_32[11];
        coefficients[10] += a_1;
        coefficients[11] = a_0 - a_1;

        let (a_0, mut a_1) = (coefficients[12], coefficients[13]);
        a_1 *= ROOTS_OF_UNITY_32[7];
        coefficients[12] += a_1;
        coefficients[13] = a_0 - a_1;

        let (a_0, mut a_1) = (coefficients[14], coefficients[15]);
        a_1 *= ROOTS_OF_UNITY_32[15];
        coefficients[14] += a_1;
        coefficients[15] = a_0 - a_1;
    }
}

fn serial_stark_prime_icrt(mut evaluations: Vec<Fq>) -> Vec<Fq> {
    assert_eq!(evaluations.len(), N);

    serial_stark_prime_icrt_in_place(&mut evaluations);

    evaluations
}

fn serial_stark_prime_icrt_in_place(evaluations: &mut [Fq]) {
    assert_eq!(evaluations.len(), D);

    {
        let (a_0, a_1) = (evaluations[0], evaluations[1]);
        evaluations[0] = a_0 + a_1;
        evaluations[1] = ROOTS_OF_UNITY_32[31] * (a_0 - a_1);

        let (a_0, a_1) = (evaluations[2], evaluations[3]);
        evaluations[2] = a_0 + a_1;
        evaluations[3] = ROOTS_OF_UNITY_32[23] * (a_0 - a_1);

        let (a_0, a_1) = (evaluations[4], evaluations[5]);
        evaluations[4] = a_0 + a_1;
        evaluations[5] = ROOTS_OF_UNITY_32[27] * (a_0 - a_1);

        let (a_0, a_1) = (evaluations[6], evaluations[7]);
        evaluations[6] = a_0 + a_1;
        evaluations[7] = ROOTS_OF_UNITY_32[19] * (a_0 - a_1);

        let (a_0, a_1) = (evaluations[8], evaluations[9]);
        evaluations[8] = a_0 + a_1;
        evaluations[9] = ROOTS_OF_UNITY_32[29] * (a_0 - a_1);

        let (a_0, a_1) = (evaluations[10], evaluations[11]);
        evaluations[10] = a_0 + a_1;
        evaluations[11] = ROOTS_OF_UNITY_32[21] * (a_0 - a_1);

        let (a_0, a_1) = (evaluations[12], evaluations[13]);
        evaluations[12] = a_0 + a_1;
        evaluations[13] = ROOTS_OF_UNITY_32[25] * (a_0 - a_1);

        let (a_0, a_1) = (evaluations[14], evaluations[15]);
        evaluations[14] = a_0 + a_1;
        evaluations[15] = ROOTS_OF_UNITY_32[17] * (a_0 - a_1);
    }

    for i in 0..(D / 8) {
        // f_1
        {
            let (coeff_i, coeff_d_div_8_plus_i) = (evaluations[i], evaluations[D / 8 + i]);

            evaluations[i] = coeff_i + coeff_d_div_8_plus_i;
            evaluations[i + D / 8] = ROOTS_OF_UNITY_32[30] * (coeff_i - coeff_d_div_8_plus_i);
        }

        // f_2
        {
            let (coeff_i, coeff_d_div_8_plus_i) =
                (evaluations[D / 4 + i], evaluations[3 * D / 8 + i]);

            evaluations[D / 4 + i] = coeff_i + coeff_d_div_8_plus_i;
            evaluations[3 * D / 8 + i] = ROOTS_OF_UNITY_32[22] * (coeff_i - coeff_d_div_8_plus_i);
        }

        // f_3
        {
            let (coeff_i, coeff_d_div_8_plus_i) =
                (evaluations[D / 2 + i], evaluations[5 * D / 8 + i]);

            evaluations[D / 2 + i] = coeff_i + coeff_d_div_8_plus_i;
            evaluations[5 * D / 8 + i] = ROOTS_OF_UNITY_32[26] * (coeff_i - coeff_d_div_8_plus_i);
        }

        // f_4
        {
            let (coeff_i, coeff_d_div_8_plus_i) =
                (evaluations[3 * D / 4 + i], evaluations[7 * D / 8 + i]);

            evaluations[3 * D / 4 + i] = coeff_i + coeff_d_div_8_plus_i;
            evaluations[7 * D / 8 + i] = ROOTS_OF_UNITY_32[18] * (coeff_i - coeff_d_div_8_plus_i);
        }
    }

    for i in 0..(D / 4) {
        // f_1
        {
            let (coeff_i, coeff_d_div_4_plus_i) = (evaluations[i], evaluations[D / 4 + i]);

            evaluations[i] = coeff_i + coeff_d_div_4_plus_i;
            evaluations[i + D / 4] = ROOTS_OF_UNITY_32[28] * (coeff_i - coeff_d_div_4_plus_i);
        }

        // f_2
        {
            let (coeff_i, coeff_d_div_4_plus_i) =
                (evaluations[D / 2 + i], evaluations[3 * D / 4 + i]);

            evaluations[D / 2 + i] = coeff_i + coeff_d_div_4_plus_i;
            evaluations[3 * D / 4 + i] = ROOTS_OF_UNITY_32[20] * (coeff_i - coeff_d_div_4_plus_i);
        }
    }

    // Rewind the first step of the CRT.
    for i in 0..(D / 2) {
        let (coeff_i, coeff_d_div_2_plus_i) = (evaluations[i], evaluations[D / 2 + i]);

        evaluations[i] = SIXTEEN_INV * (coeff_i + coeff_d_div_2_plus_i);
        evaluations[i + D / 2] =
            SIXTEEN_INV_TIMES_ROOT_OF_UNITY_32_24 * (coeff_i - coeff_d_div_2_plus_i);
    }
}

#[cfg(test)]
mod tests {
    use ark_ff::{Field, MontFp, UniformRand};
    use ark_std::Zero;
    use rand::thread_rng;

    use super::*;

    #[test]
    fn test_roots_of_unity() {
        // Check if they all are roots of unity of degree 32.
        for root in ROOTS_OF_UNITY_32 {
            assert_eq!(root.pow([32]), Fq::ONE);
        }

        // Check if they are pairwise distinct.
        for (i, item) in ROOTS_OF_UNITY_32.iter().enumerate().take(24) {
            for item2 in ROOTS_OF_UNITY_32.iter().take(24).skip(i + 1) {
                assert_ne!(*item, *item2);
            }
        }

        // Check if they come in order.
        for i in 0..32u64 {
            assert_eq!(ROOTS_OF_UNITY_32[i as usize], ROOTS_OF_UNITY_32[1].pow([i]));
        }
    }

    #[test]
    fn test_crt() {
        // 1 + 2 * X + 3 * X^2 + 15 * X^15
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
            MontFp!("645567075879839201317982428175614860863887056415581724161883453888590511158"),
            MontFp!("954581386842844617990350588793859972457388686372942226298498503178516390686"),
            MontFp!("3256862904360112534622944024048777964481830288622596035686455700740188841883"),
            MontFp!("2379994210249466073463368525171887413443108399252073413799346454464448297239"),
            MontFp!("1154787450473232990877007307548392072878001948725547466967243997065670930781"),
            MontFp!("127044796229026925925676336941612140641682969690926841605253822315415053191"),
            MontFp!("2264291370407632332462551969004321602019453100516039364350823369215509801549"),
            MontFp!("72379171556238964432087169600744290083969196399083027049770867539276234964"),
            MontFp!("3392988301432672676654143254540743340205316729651753579425056524920715568596"),
            MontFp!("345623899732466539114181911245544957171813978741334524852816912705261753883"),
            MontFp!("3185268477836746232513995247235360610631152596074905508986422324167448586405"),
            MontFp!("313124898330376979112325153168491303237931126195199786681888350478318132082"),
            MontFp!("3252710104456412601428779891144119168067847513403787243699057366054350147707"),
            MontFp!("2811893809428838290900969955841763465542245498752761173309728642815514980559"),
            MontFp!("1345201406880623537768056884941378708623964668432682170554896434038851699412"),
            MontFp!("3445703045232519210994161617357948974635263965405559512355593725498899233769"),
        ];

        serial_stark_prime_crt_in_place(&mut test_poly);

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
            MontFp!("3342128707467438650582796293028326366188841359173385550320192593983949514781"),
            MontFp!("145415729351446532104221590506099754017149096611578597467764922845836744024"),
            MontFp!("28344076774352150081365919378257754937378686670188713656703512828535506046"),
            MontFp!("102614275072893880928938980182386230479738072876443838528431026477550264998"),
            MontFp!("3033695382336505385241782670919441783483213948688819507086473908690567613286"),
            MontFp!("939956913328774515564849172117096646347233743420544100796619190984456389268"),
            MontFp!("1637409282267998787972433807386232485827432500117030096420080018113277841431"),
            MontFp!("1625943999398983738615579915767369295588334238436799695643010994483442206345"),
            MontFp!("3594911410986703872812844193278353776989462412403241940436425945741340009951"),
            MontFp!("75507220001790948222892434565881002017963472150232301947099166071873513613"),
            MontFp!("1219346410107661887858150731282067078930201967594791163094249908344089354239"),
            MontFp!("2347240536236105718500758207063838353308586578514927994468409092114441172527"),
            MontFp!("3374920525187094988736089450950872833267621485172465575429173175774334119406"),
            MontFp!("1711209537345756975647900555311037454790763165006988461096559542358841336444"),
            MontFp!("1505297521723099557839333110912200164930605895060277511805849845124000458838"),
            MontFp!("645577993076310905171322449016029758257223885423461851614601549014568135642"),
        ];

        serial_stark_prime_crt_in_place(&mut test_poly);

        assert_eq!(test_poly, expected);
    }

    #[test]
    fn test_icrt() {
        // 1 + 2 * X + 3 * X^2 + 15 * X^15
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
            MontFp!("645567075879839201317982428175614860863887056415581724161883453888590511158"),
            MontFp!("954581386842844617990350588793859972457388686372942226298498503178516390686"),
            MontFp!("3256862904360112534622944024048777964481830288622596035686455700740188841883"),
            MontFp!("2379994210249466073463368525171887413443108399252073413799346454464448297239"),
            MontFp!("1154787450473232990877007307548392072878001948725547466967243997065670930781"),
            MontFp!("127044796229026925925676336941612140641682969690926841605253822315415053191"),
            MontFp!("2264291370407632332462551969004321602019453100516039364350823369215509801549"),
            MontFp!("72379171556238964432087169600744290083969196399083027049770867539276234964"),
            MontFp!("3392988301432672676654143254540743340205316729651753579425056524920715568596"),
            MontFp!("345623899732466539114181911245544957171813978741334524852816912705261753883"),
            MontFp!("3185268477836746232513995247235360610631152596074905508986422324167448586405"),
            MontFp!("313124898330376979112325153168491303237931126195199786681888350478318132082"),
            MontFp!("3252710104456412601428779891144119168067847513403787243699057366054350147707"),
            MontFp!("2811893809428838290900969955841763465542245498752761173309728642815514980559"),
            MontFp!("1345201406880623537768056884941378708623964668432682170554896434038851699412"),
            MontFp!("3445703045232519210994161617357948974635263965405559512355593725498899233769"),
        ];

        serial_stark_prime_icrt_in_place(&mut evaluations);

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
            MontFp!("3342128707467438650582796293028326366188841359173385550320192593983949514781"),
            MontFp!("145415729351446532104221590506099754017149096611578597467764922845836744024"),
            MontFp!("28344076774352150081365919378257754937378686670188713656703512828535506046"),
            MontFp!("102614275072893880928938980182386230479738072876443838528431026477550264998"),
            MontFp!("3033695382336505385241782670919441783483213948688819507086473908690567613286"),
            MontFp!("939956913328774515564849172117096646347233743420544100796619190984456389268"),
            MontFp!("1637409282267998787972433807386232485827432500117030096420080018113277841431"),
            MontFp!("1625943999398983738615579915767369295588334238436799695643010994483442206345"),
            MontFp!("3594911410986703872812844193278353776989462412403241940436425945741340009951"),
            MontFp!("75507220001790948222892434565881002017963472150232301947099166071873513613"),
            MontFp!("1219346410107661887858150731282067078930201967594791163094249908344089354239"),
            MontFp!("2347240536236105718500758207063838353308586578514927994468409092114441172527"),
            MontFp!("3374920525187094988736089450950872833267621485172465575429173175774334119406"),
            MontFp!("1711209537345756975647900555311037454790763165006988461096559542358841336444"),
            MontFp!("1505297521723099557839333110912200164930605895060277511805849845124000458838"),
            MontFp!("645577993076310905171322449016029758257223885423461851614601549014568135642"),
        ];

        serial_stark_prime_icrt_in_place(&mut evaluations);

        assert_eq!(evaluations, expected);
    }

    fn test_crt_icrt() {
        let mut rng = thread_rng();
        let coefficients: Vec<Fq> = (0..16).map(|_| Fq::rand(&mut rng)).collect();

        let mut input = coefficients.clone();

        serial_stark_prime_crt_in_place(&mut input);
        serial_stark_prime_icrt_in_place(&mut input);

        assert_eq!(coefficients, input);
    }

    #[test]
    fn test_crt_icrt_1000000_times() {
        for _i in 0..1000000 {
            test_crt_icrt();
        }
    }
}
