use ark_serialize::{CanonicalSerialize, Compress};
use nimue::IOPattern;
use num_traits::Zero;

use crate::lattice_arithmetic::challenge_set::ternary::{TernaryChallengeSet, Trit};
use crate::lattice_arithmetic::matrix::Matrix;
use crate::lattice_arithmetic::ring::Fq;
use crate::lattice_arithmetic::traits::{Modulus, WithL2Norm};
use crate::lova::prover::prove_folding;
use crate::lova::util::{CRS, SECPARAM};
use crate::lova::verifier::{decide, verify_folding};
use crate::nimue::iopattern::SqueezeFromRandomBytes;

const Q: u64 = (1 << 8) + 1;

type F = Fq<Q>;

const N: usize = (1 << 10);

#[test]
fn test_prover() {
    let crs = CRS::new(N, (F::modulus() as f64).sqrt());
    let b = crs.norm_bound.sqrt().round_ties_even() as u128;
    let k = crs.norm_bound.log(b as f64).ceil() as usize;
    let m = crs.commitment_mat.nrows();

    let s = F::zero().serialized_size(Compress::Yes);

    let mut io = IOPattern::new("lova")
        .absorb((2 * SECPARAM) * k * m * s, "t")
        .absorb(((2 * SECPARAM) * (2 * SECPARAM + 1)).div_ceil(2) * k * k, "g")
        .squeeze_matrices::<Trit, TernaryChallengeSet<Trit>>(SECPARAM, 2 * SECPARAM, k, "c");

    let mut arthur = io.to_arthur();

    let s = Matrix::<F>::zeros(N, 2 * SECPARAM);
    let t = &crs.commitment_mat * &s;

    // Fold
    let s_new = prove_folding(&mut arthur, &crs, &t, &s).unwrap();
    let folding_proof = arthur.transcript();
    let mut merlin = io.to_merlin(folding_proof);
    let t_new = verify_folding(&mut merlin, &crs, &t).unwrap();

    // Check that (t_new, s_new) is in the relation
    assert_eq!(t_new, &crs.commitment_mat * &s_new);
    for s_i in s_new.iter() {
        assert!(s_i.l2_norm_squared() as f64 <= crs.norm_bound);

    }

    // Check that decider accepts
    let verifies = decide(&mut merlin, &crs, &t_new, &s_new);
    assert!(verifies.is_ok());
}