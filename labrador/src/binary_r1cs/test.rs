use ark_r1cs_std::alloc::AllocVar;
use ark_r1cs_std::boolean::Boolean;
use indicatif::ProgressBar;
use log::debug;
use nimue::IOPattern;
use tracing_subscriber::fmt::format;
use tracing_subscriber::fmt::format::FmtSpan;

use lattirust_arithmetic::ntt::ntt_modulus;
use lattirust_arithmetic::ring::Pow2CyclotomicPolyRingNTT;

use crate::binary_r1cs::prover::prove_binary_r1cs;
use crate::binary_r1cs::util::{BinaryR1CSCRS, Z2};
use crate::binary_r1cs::verifier::verify_binary_r1cs;
use crate::iopattern::LabradorIOPattern;

const Q: u64 = ntt_modulus::<64>(32);
const D: usize = 64;

type R = Pow2CyclotomicPolyRingNTT<Q, 64>;

fn init() {
    // let _ = env_logger::builder().is_test(true).try_init();
    tracing_subscriber::fmt::fmt()
        .with_span_events(FmtSpan::ENTER | FmtSpan::CLOSE)
        .event_format(format().compact())
        .with_env_filter("none,labrador=trace")
        .init();
}

#[test]
#[allow(non_snake_case)]
fn test_prove_binary_r1cs() {
    init();

    let num_r1cs_constraints = 1 << 12;
    let k = (num_r1cs_constraints / D) * D;

    debug!("Constructing constraint system...");
    // Generate a dummy constraint system with num_r1cs_constraints constraints and variables.
    let cs = ark_relations::r1cs::ConstraintSystem::<Z2>::new_ref();
    let pb = ProgressBar::new(k as u64);
    let v = Boolean::new_witness(cs.clone(), || Ok(true)).unwrap();
    cs.enforce_constraint(v.lc(), v.lc(), v.lc()).unwrap();
    for i in 2..k {
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
    debug_assert_eq!(cs.is_satisfied(), Ok(true));

    debug!("Constructing common reference string...");
    let crs = BinaryR1CSCRS::<R>::new(k, k);

    debug!("Setting IOPattern...");
    let io = IOPattern::new("labrador_binaryr1cs")
        .labrador_binaryr1cs_io(&cs, &crs)
        .ratchet()
        // .labrador_crs(&crs.pr_crs())
        // .ratchet()
        // .labrador_instance(&PrincipalRelation::<R>::new_empty(&crs.pr_crs()))
        // .ratchet()
        .labrador_io(&crs.core_crs);
    let mut arthur = io.to_arthur();

    debug!("Proving...");
    let proof = prove_binary_r1cs(&crs, &mut arthur, &cs).unwrap();
    debug!("Finished proving, proof size = {} bytes", proof.len());

    let mut merlin = io.to_merlin(proof);
    verify_binary_r1cs(&mut merlin, &cs, &crs).unwrap();
}
