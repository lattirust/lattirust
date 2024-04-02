use ark_serialize::{CanonicalSerialize, Compress};
use nimue::IOPattern;
use num_traits::Zero;

use crate::labrador::util::inner_products_mat;
use crate::lattice_arithmetic::challenge_set::ternary::{TernaryChallengeSet, Trit};
use crate::lattice_arithmetic::matrix::Matrix;
use crate::lattice_arithmetic::ring::Fq;
use crate::lattice_arithmetic::traits::Modulus;
use crate::lova::prover::prove_folding;
use crate::lova::util::{BaseRelation, Instance, PublicParameters, SECPARAM};
use crate::lova::verifier::verify_folding;
use crate::nimue::iopattern::SqueezeFromRandomBytes;
use crate::relations::traits::Relation;

const Q: u64 = (1 << 8) + 1;

type F = Fq<Q>;

const N: usize = (1 << 9);

#[test]
fn test() {
    let pp = PublicParameters::new(N, (F::modulus() as f64).sqrt());
    let b = pp.norm_bound.sqrt().round_ties_even() as u128;
    let k = pp.norm_bound.log(b as f64).ceil() as usize;
    let m = pp.commitment_mat.nrows();

    let s = F::zero().serialized_size(Compress::Yes);

    let io = IOPattern::new("lova")
        .absorb((2 * SECPARAM) * k * m * s, "t")
        .absorb(((2 * SECPARAM) * (2 * SECPARAM + 1)).div_ceil(2) * k * k, "g")
        .squeeze_matrices::<Trit, TernaryChallengeSet<Trit>>(SECPARAM, 2 * SECPARAM, k, "c");

    let mut arthur = io.to_arthur();

    let witness = Matrix::<F>::zeros(N, 2 * SECPARAM);
    let instance = Instance { commitment: &pp.commitment_mat * &witness, inner_products: inner_products_mat(&witness).into() };

    // Fold
    let new_witness = prove_folding(&mut arthur, &pp, &instance, &witness).unwrap();
    let folding_proof = arthur.transcript();
    let mut merlin = io.to_merlin(folding_proof);
    let new_instance = verify_folding(&mut merlin, &pp, &instance).unwrap();

    // Check that the folded instance and witness are in the relation
    assert!(BaseRelation::is_satisfied(&pp, &new_instance, &new_witness));
}