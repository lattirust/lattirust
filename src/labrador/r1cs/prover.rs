#![allow(non_snake_case)]

use anyhow::Error;
use ark_ff::{Field, Fp, MontBackend};
use ark_ff::fields::MontConfig;

use crate::labrador::binary_r1cs::util::R1CSWitness;
use crate::labrador::r1cs::util::{R1CSCRS, R1CSInstance};
use crate::labrador::util::{concat, flatten_vec_vector};
use crate::lattice_arithmetic::challenge_set::labrador_challenge_set::LabradorChallengeSet;
use crate::lattice_arithmetic::challenge_set::weighted_ternary::WeightedTernaryChallengeSet;
use crate::lattice_arithmetic::matrix::Vector;
use crate::lattice_arithmetic::poly_ring::PolyRing;
use crate::lattice_arithmetic::ring::Fq;
use crate::lattice_arithmetic::traits::FromRandomBytes;
use crate::nimue::arthur::LatticeArthur;

#[derive(MontConfig)]
#[modulus = "18446744073709551617"] // 2^64+1
#[generator = "3"]
pub struct FqConfig;

pub type F64b = Fp<MontBackend<FqConfig, 2>, 2>;
// pub type Z64 = Fq<18446744073709551617>;
pub type Z64 = Fq<3>; // TODO

fn enc<R: PolyRing>(vec: &Vector<Z64>) -> Vector<R> {
    todo!()
}

struct R1CSTranscript<R: PolyRing> {
    t: Vec<Vector<R>>,
    phi: Vec<Vector<Z64>>,
    // t_d: _,
    // c: OMatrix<Z64, Dyn, U1>,
    // g: _,
    // r: _,
}

pub fn prove_r1cs<'a, R: PolyRing>(crs: &R1CSCRS<R>, arthur: &'a mut LatticeArthur<R>, instance: &R1CSInstance<Z64>, witness: &R1CSWitness<Z64>) -> Result<&'a [u8], Error>
    where LabradorChallengeSet<R>: FromRandomBytes<R>, WeightedTernaryChallengeSet<R>: FromRandomBytes<R>
{
    let (A, B, C) = (&instance.A, &instance.B, &instance.C);
    let w = &witness.w;
    let (k, n) = (A.nrows(), A.ncols());
    assert_eq!(k, n, "the current implementation only support k = n"); // TODO: remove this restriction by splitting a,b,c or w into multiple vectors
    let l = crs.l;

    let a = A * w;
    let b = B * w;
    let c = C * w;

    let a_R = enc::<R>(&a);
    let b_R = enc::<R>(&b);
    let c_R = enc::<R>(&c);
    let w_R = enc::<R>(w);

    let v = concat(&[&a_R, &b_R, &c_R, &w_R]);
    let t = &crs.A * &v;

    arthur.absorb_vector(&t).unwrap();

    let phi = arthur.squeeze_vectors::<Z64, Z64>(l, k)?;

    let mut d_R = Vec::<Vector::<R>>::with_capacity(l);
    for i in 0..l {
        d_R.push(enc::<R>(&phi[i].component_mul(&a)));
    }
    let d_R_flat = flatten_vec_vector(&d_R);
    let t_d = &crs.B * &d_R_flat;

    arthur.absorb_vector(&t_d)?;

    // let alpha = arthur.squeeze_vector::<Z64, Z64>()
    let c = arthur.squeeze_vectors::<Z64, Z64>(l, k * (l + 3) + l);

    let g = (0..l).map(|i|
        // f(i, a_R, b_R, c_R, w_R, d_R_flat)
        i
    ).collect::<Vec<_>>();

    // arthur.absorb_vec(&g)?;
    todo!();

    // let transcript = R1CSTranscript { t, phi, t_d, c, g, r };

    // let instance_pr = reduce(&crs, &instance, &transcript);

    // let witness_pr = Witness::<R> {
    //     s: vec![a_R, b_R, c_R, w_R, d_R] // see definition of indices above
    // };

    // arthur.ratchet()?;
    // prove_principal_relation(arthur, &instance_pr, &witness_pr, &crs.pr_crs())
    Ok(arthur.transcript())
}