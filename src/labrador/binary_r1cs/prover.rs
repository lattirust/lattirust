#![allow(non_snake_case)]

use anyhow::Error;
use log::debug;

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
    debug!("labrador::binary_r1cs starting BinR1CS -> PrincipalRelation reduction");
    print!("__labrador::binary_r1cs starting BinR1CS -> PrincipalRelation reduction");
    debug_assert!(is_satisfied(&crs, &instance, &witness));
    let (A, B, C) = (&instance.A, &instance.B, &instance.C);
    let w = &witness.w;
    let (k, n) = (A.nrows(), A.ncols());
    assert_eq!(k, n, "the current implementation only support k = n"); // TODO: remove this restriction by splitting a,b,c or w into multiple vectors

    // TODO: add statement

    let a = A * w;
    let b = B * w;
    let c = C * w;

    debug!("labrador::binary_r1cs computed A*w, B*w, C*w");

    let a_R = lift::<R>(&a);
    let b_R = lift::<R>(&b);
    let c_R = lift::<R>(&c);
    let w_R = lift::<R>(w);

    let v = concat(&[&a_R, &b_R, &c_R, &w_R]);
    let t = &crs.A * &v;

    arthur.absorb_vector(&t).unwrap();

    debug!("labrador::binary_r1cs squeeze alpha in {{0,1}}^{SECPARAM}x{k}");
    let alpha = arthur.squeeze_binary_matrix(SECPARAM, k).unwrap();
    debug!("labrador::binary_r1cs squeeze beta in {{0,1}}^{SECPARAM}x{k}");
    let beta = arthur.squeeze_binary_matrix(SECPARAM, k).unwrap();
    debug!("labrador::binary_r1cs squeeze gamma in {{0,1}}^{SECPARAM}x{k}");
    let gamma = arthur.squeeze_binary_matrix(SECPARAM, k).unwrap();

    // delta_i is computed mod 2, i.e., over Z2
    let delta = &alpha * A + &beta * B + &gamma * C;

    // g_i is computed over Zq
    let (a_BR, b_BR, c_BR) = (a.map(embed::<R::BaseRing>), b.map(embed::<R::BaseRing>), c.map(embed::<R::BaseRing>));
    let (alpha_BR, beta_BR, gamma_BR) = (alpha.map(embed::<R::BaseRing>), beta.map(embed::<R::BaseRing>), gamma.map(embed::<R::BaseRing>));
    let (delta_BR, w_BR) = (delta.map(embed::<R::BaseRing>), w.map(embed::<R::BaseRing>));

    let g = &alpha_BR * &a_BR + &beta_BR * &b_BR + &gamma_BR * &c_BR - &delta_BR * &w_BR;

    debug_assert_eq!(g.len(), SECPARAM, "g has length {} but should have length {}", g.len(), SECPARAM);
    debug!("labrador::binary_r1cs absorb g in R_q^{SECPARAM}");
    arthur.absorb_vector_baseringlem(&g).unwrap();

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

    debug!("labrador::binary_r1cs finished BinR1CS -> PrincipalRelation reduction");
    prove_principal_relation(arthur, &instance_pr, &witness_pr, &crs.pr_crs())
}

