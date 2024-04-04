use nimue::IOPattern;
use pretty_env_logger::env_logger;

use crate::lattice_arithmetic::matrix::Matrix;
use crate::lattice_arithmetic::ring::Fq;
use crate::lattice_arithmetic::traits::Modulus;
use crate::lova::prover::prove_folding;
use crate::lova::util::{BaseRelation, Instance, LovaIOPattern, PublicParameters, SECPARAM};
use crate::lova::verifier::verify_folding;
use crate::relations::traits::Relation;

const Q: u64 = (1 << 7) - 1;

type F = Fq<Q>;

const N: usize = (1 << 6);

fn init() {
    let _ = env_logger::builder().is_test(true).try_init();
}
#[test]
fn test() {
    init();

    let pp = PublicParameters::new(N, (F::modulus() as f64).sqrt());

    let witness_1 = Matrix::<F>::zeros(N, SECPARAM);
    let instance_1 = Instance::new(&pp, &witness_1);
    debug_assert!(BaseRelation::is_satisfied(&pp, &instance_1, &witness_1));

    let witness_2 = Matrix::<F>::identity(N, SECPARAM);
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