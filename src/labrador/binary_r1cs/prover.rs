#![allow(non_snake_case)]

use anyhow::Error;

use crate::labrador::binary_r1cs::util::*;
use crate::labrador::prover::prove_principal_relation;
use crate::labrador::util::concat;
use crate::labrador::witness::Witness;
use crate::lattice_arithmetic::challenge_set::labrador_challenge_set::LabradorChallengeSet;
use crate::lattice_arithmetic::challenge_set::weighted_ternary::WeightedTernaryChallengeSet;
use crate::lattice_arithmetic::poly_ring::PolyRing;
use crate::lattice_arithmetic::traits::FromRandomBytes;
use crate::nimue::arthur::LatticeArthur;

pub fn prove_binary_r1cs<'a, R: PolyRing>(crs: &BinaryR1CSCRS<R>, arthur: &'a mut LatticeArthur<R>, instance: &BinaryR1CSInstance, witness: &BinaryR1CSWitness) -> Result<&'a [u8], Error>
    where LabradorChallengeSet<R>: FromRandomBytes<R>, WeightedTernaryChallengeSet<R>: FromRandomBytes<R>
{
    debug_assert!(is_satisfied(&crs, &instance, &witness));
    let (A, B, C) = (&instance.A, &instance.B, &instance.C);
    let w = &witness.w;
    let (k, n) = (A.nrows(), A.ncols());
    assert_eq!(k, n, "the current implementation only support k = n"); // TODO: remove this restriction by splitting a,b,c or w into multiple vectors

    // TODO: add statement

    let a = A * w;
    let b = B * w;
    let c = C * w;

    let a_R = lift::<R>(&a);
    let b_R = lift::<R>(&b);
    let c_R = lift::<R>(&c);
    let w_R = lift::<R>(w);

    let v = concat(&[&a_R, &b_R, &c_R, &w_R]);
    let t = &crs.A * &v;

    arthur.absorb_vector(&t).unwrap();

    let alpha = arthur.squeeze_binary_matrix(SECPARAM, k).unwrap();
    let beta = arthur.squeeze_binary_matrix(SECPARAM, k).unwrap();
    let gamma = arthur.squeeze_binary_matrix(SECPARAM, k).unwrap();

    // delta_i is computed mod 2, i.e., over Z2
    let delta = &alpha * A + &beta * B + &gamma * C;

    // g_i is computed over Zq
    let (a_BR, b_BR, c_BR) = (a.map(embed::<R::BaseRing>), b.map(embed::<R::BaseRing>), c.map(embed::<R::BaseRing>));
    let (alpha_BR, beta_BR, gamma_BR) = (alpha.map(embed::<R::BaseRing>), beta.map(embed::<R::BaseRing>), gamma.map(embed::<R::BaseRing>));
    let (delta_BR, w_BR) = (delta.map(embed::<R::BaseRing>), w.map(embed::<R::BaseRing>));

    let g = &alpha_BR * &a_BR + &beta_BR * &b_BR + &gamma_BR * &c_BR - &delta_BR * &w_BR;
    arthur.absorb_vector(&g).unwrap();

    let a_tilde = R::sigma_vec(&a_R);
    let b_tilde = R::sigma_vec(&b_R);
    let c_tilde = R::sigma_vec(&c_R);
    let w_tilde = R::sigma_vec(&w_R);

    let transcript = BinaryR1CSTranscript { t, alpha, beta, gamma, g, delta };

    let instance_pr = reduce(&crs, &instance, &transcript);

    let witness_pr = Witness::<R> {
        s: vec![a_R, b_R, c_R, w_R, a_tilde, b_tilde, c_tilde, w_tilde] // see definition of indices above
    };

    arthur.ratchet()?;
    prove_principal_relation(arthur, &instance_pr, &witness_pr, &crs.pr_crs())
}

