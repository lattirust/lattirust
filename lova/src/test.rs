use lattirust_arithmetic::linear_algebra::Matrix;
use lattirust_arithmetic::ring::Fq;
use nimue::IOPattern;
use pretty_env_logger::env_logger;
use relations::traits::Relation;

use crate::prover::prove_folding;
use crate::util::{BaseRelation, Instance, LovaIOPattern, OptimizationMode, PublicParameters};
use crate::verifier::verify_folding;

const Q: u64 = ((1u128 << 64) - 1) as u64;

type F = Fq<Q>;

const N: usize = 1 << 4;

fn init() {
    let _ = env_logger::builder().is_test(true).try_init();
}

#[test]
fn test() {
    init();

    let pp = PublicParameters::new(N, OptimizationMode::OptimizeForSpeed);

    let witness_1 = Matrix::<F>::zeros(N, pp.security_parameter);
    let instance_1 = Instance::new(&pp, &witness_1);
    debug_assert!(BaseRelation::is_satisfied(&pp, &instance_1, &witness_1));

    let witness_2 = Matrix::<F>::identity(N, pp.security_parameter);
    let instance_2 = Instance::new(&pp, &witness_2);
    debug_assert!(BaseRelation::is_satisfied(&pp, &instance_2, &witness_2));

    let io = IOPattern::new("lova").folding_round(&pp);

    // Prove folding
    let mut arthur = io.to_arthur();
    let new_witness = prove_folding(&mut arthur, &pp, &instance_1, &witness_1, &instance_2, &witness_2).unwrap();
    let folding_proof = arthur.transcript();

    // Verify folding
    let mut merlin = io.to_merlin(folding_proof);
    let new_instance = verify_folding(&mut merlin, &pp, &instance_1, &instance_2).unwrap();

    // Check that the folded instance and witness are in the relation
    assert!(BaseRelation::is_satisfied(&pp, &new_instance, &new_witness));
}