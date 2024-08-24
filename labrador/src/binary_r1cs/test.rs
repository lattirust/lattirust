use ark_r1cs_std::alloc::AllocVar;
use ark_r1cs_std::boolean::Boolean;
use ark_relations::r1cs::ConstraintSystemRef;
use indicatif::ProgressBar;
use log::debug;
use nimue::IOPattern;
use tracing_subscriber::fmt::format;
use tracing_subscriber::fmt::format::FmtSpan;

use lattirust_arithmetic::ntt::ntt_modulus;
use lattirust_arithmetic::ring::Pow2CyclotomicPolyRingNTT;

use crate::binary_r1cs::prover::prove_reduction_binaryr1cs_labradorpr;
use crate::binary_r1cs::util::{BinaryR1CSCRS, Z2};
use crate::binary_r1cs::verifier::verify_reduction_binaryr1cs_labradorpr;
use crate::iopattern::LabradorIOPattern;

const Q: u64 = ntt_modulus::<64>(32);
const D: usize = 64;

type R = Pow2CyclotomicPolyRingNTT<Q, D>;

fn init() {
    let _ = tracing_subscriber::fmt::fmt()
        .with_span_events(FmtSpan::ENTER | FmtSpan::CLOSE)
        .event_format(format().compact())
        .with_env_filter("none,labrador=trace")
        .try_init();
}

fn valid_binr1cs_instance_witness(num_r1cs_constraints: usize) -> ConstraintSystemRef<Z2> {
    let cs = ark_relations::r1cs::ConstraintSystem::<Z2>::new_ref();
    let pb = ProgressBar::new(num_r1cs_constraints as u64);
    let v = Boolean::new_witness(cs.clone(), || Ok(true)).unwrap();
    cs.enforce_constraint(v.lc(), v.lc(), v.lc()).unwrap();
    for i in 2..num_r1cs_constraints {
        let _v = Boolean::new_witness(cs.clone(), || Ok(i % 2 == 0)).unwrap();
        // cs.enforce_constraint(v.lc(), v.lc(), v.lc()).unwrap();
        // v.enforce_equal(&v).unwrap();
        pb.inc(1);
    }
    cs.finalize();
    let mats = cs.to_matrices().unwrap();
    debug!(
        "\t{} ≈ 2^{} constraints, {}+{} ≈ 2^{} variables, ({}, {}, {}) non-zero entries",
        cs.num_constraints(),
        cs.num_constraints().next_power_of_two().ilog2(),
        cs.num_instance_variables(),
        cs.num_witness_variables(),
        (cs.num_instance_variables() + cs.num_witness_variables())
            .next_power_of_two()
            .ilog2(),
        mats.a_num_non_zero,
        mats.b_num_non_zero,
        mats.c_num_non_zero
    );
    debug_assert_eq!(cs.num_constraints(), num_r1cs_constraints);
    debug_assert_eq!(cs.num_instance_variables(), 1);
    debug_assert_eq!(cs.num_witness_variables(), num_r1cs_constraints - 1);
    // debug_assert!(mats.a_num_non_zero > 0);
    // debug_assert!(mats.b_num_non_zero > 0);
    // debug_assert!(mats.c_num_non_zero > 0);
    cs
}

fn invalid_binr1cs_instance_witness(num_r1cs_constraints: usize) -> ConstraintSystemRef<Z2> {
    let cs = ark_relations::r1cs::ConstraintSystem::<Z2>::new_ref();
    let pb = ProgressBar::new(num_r1cs_constraints as u64);
    let v = Boolean::new_witness(cs.clone(), || Ok(true)).unwrap();
    cs.enforce_constraint(
        v.lc(),
        Boolean::constant(false).lc(),
        Boolean::constant(true).lc(),
    )
    .unwrap(); // v * 0 == 1 cannot be satisfied
    for i in 2..num_r1cs_constraints {
        let _v = Boolean::new_witness(cs.clone(), || Ok(i % 2 == 0)).unwrap();
        // cs.enforce_constraint(v.lc(), v.lc(), v.lc()).unwrap();
        // v.enforce_equal(&v).unwrap();
        pb.inc(1);
    }
    cs.finalize();
    let mats = cs.to_matrices().unwrap();
    debug!(
        "\t{} ≈ 2^{} constraints, {}+{} ≈ 2^{} variables, ({}, {}, {}) non-zero entries",
        cs.num_constraints(),
        cs.num_constraints().next_power_of_two().ilog2(),
        cs.num_instance_variables(),
        cs.num_witness_variables(),
        (cs.num_instance_variables() + cs.num_witness_variables())
            .next_power_of_two()
            .ilog2(),
        mats.a_num_non_zero,
        mats.b_num_non_zero,
        mats.c_num_non_zero
    );
    debug_assert_eq!(cs.num_constraints(), num_r1cs_constraints);
    debug_assert_eq!(cs.num_instance_variables(), 1);
    debug_assert_eq!(cs.num_witness_variables(), num_r1cs_constraints - 1);
    // debug_assert!(mats.a_num_non_zero > 0);
    // debug_assert!(mats.b_num_non_zero > 0);
    // debug_assert!(mats.c_num_non_zero > 0);
    cs
}

#[test]
fn test_reduction_binaryr1cs_principalrelation_completeness() {
    init();

    let mut num_r1cs_constraints = 1 << 12;
    num_r1cs_constraints = (num_r1cs_constraints / D) * D;

    let cs = valid_binr1cs_instance_witness(num_r1cs_constraints);
    debug_assert_eq!(
        cs.is_satisfied(),
        Ok(true),
        "generated R1CS is not satisfied, aborting test"
    );

    let crs = BinaryR1CSCRS::<R>::new(num_r1cs_constraints, num_r1cs_constraints);

    let io = IOPattern::new("labrador_binaryr1cs").labrador_binaryr1cs_io(&cs, &crs);

    let mut merlin = io.to_merlin();
    let (pr_instance, pr_witness) = prove_reduction_binaryr1cs_labradorpr(&crs, &mut merlin, &cs);

    debug_assert!(pr_instance.is_valid_witness(&pr_witness), "reduction is not complete; the output PrincipalRelation witness is not a valid witness for the output instance");

    let proof = merlin.transcript();
    debug!("Finished proving, proof size = {} bytes", proof.len());

    let mut arthur = io.to_arthur(proof);
    let verifier_res = verify_reduction_binaryr1cs_labradorpr(&mut arthur, &cs, &crs);

    debug_assert!(
        verifier_res.is_ok(),
        "reduction is not complete; the verifier rejected the proof"
    );
    debug_assert_eq!(
        verifier_res.unwrap(),
        pr_instance,
        "reduction is not complete; the prover and verifier output different instances"
    );
}

#[test]
fn test_reduction_binaryr1cs_principalrelation_soundness() {
    init();

    let mut num_r1cs_constraints = 1 << 12;
    num_r1cs_constraints = (num_r1cs_constraints / D) * D;

    let cs = invalid_binr1cs_instance_witness(num_r1cs_constraints);
    debug_assert_eq!(
        cs.is_satisfied(),
        Ok(false),
        "generated R1CS instance is satisfied, aborting test"
    );

    let crs = BinaryR1CSCRS::<R>::new(num_r1cs_constraints, num_r1cs_constraints);

    let io = IOPattern::new("labrador_binaryr1cs").labrador_binaryr1cs_io(&cs, &crs);

    let mut merlin = io.to_merlin();
    let (pr_instance, pr_witness) = prove_reduction_binaryr1cs_labradorpr(&crs, &mut merlin, &cs);

    debug_assert!(!pr_instance.is_valid_witness(&pr_witness), "reduction is not sound; the output PrincipalRelation witness is a valid witness for the output instance");

    let proof = merlin.transcript();
    debug!("Finished proving, proof size = {} bytes", proof.len());

    let mut arthur = io.to_arthur(proof);
    let verifier_res = verify_reduction_binaryr1cs_labradorpr(&mut arthur, &cs, &crs);

    debug_assert!(
        verifier_res.is_err(),
        "reduction is not sound; the verifier accepted the proof"
    );
}
