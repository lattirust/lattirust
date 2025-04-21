#![allow(non_snake_case)]

use std::fmt::Debug;

use nimue::{Merlin, ProofResult};
use rayon::prelude::*;

use lattirust_arithmetic::challenge_set::labrador_challenge_set::LabradorChallengeSet;
use lattirust_arithmetic::challenge_set::weighted_ternary::WeightedTernaryChallengeSet;
use lattirust_arithmetic::decomposition::balanced_decomposition::decompose_vec_vector_dimfirst;
use lattirust_arithmetic::decomposition::decomposition::decompose_vec_polyring;
use lattirust_arithmetic::decomposition::DecompositionFriendlySignedRepresentative;
use lattirust_arithmetic::linear_algebra::inner_products::{inner_products, inner_products2};
use lattirust_arithmetic::linear_algebra::Vector;
use lattirust_arithmetic::nimue::merlin::SerMerlin;
use lattirust_arithmetic::nimue::traits::ChallengeFromRandomBytes;
use lattirust_arithmetic::ring::PolyRing;
use lattirust_arithmetic::ring::representatives::WithSignedRepresentative;
use lattirust_arithmetic::traits::FromRandomBytes;
use relations::principal_relation::{Index, Instance, Witness};

use crate::common_reference_string::CommonReferenceString;
use crate::shared::{
    compute_a__, compute_phi, compute_phi__, fold_instance, Layouter, TranscriptView,
};
use crate::util::*;

pub fn prove_principal_relation_oneround<'a, R: PolyRing>(
    merlin: &'a mut Merlin,
    crs: &CommonReferenceString<R>,
    index: &Index<R>,
    instance: &Instance<R>,
    witness: &Witness<R>,
) -> ProofResult<(Index<R>, Instance<R>, Witness<R>)>
where
    LabradorChallengeSet<R>: FromRandomBytes<R>,
    WeightedTernaryChallengeSet<R>: FromRandomBytes<R>,
    <R as PolyRing>::BaseRing: WithSignedRepresentative,
    <R::BaseRing as WithSignedRepresentative>::SignedRepresentative:
        DecompositionFriendlySignedRepresentative,
    <R as TryFrom<u128>>::Error: Debug,
{
    // Check dimensions and norms
    debug_assert!(index.is_wellformed_instance(instance).is_ok());
    debug_assert!(index.is_wellformed_witness(witness).is_ok());

    // Initialize Fiat-Shamir transcript
    // TODO: add the following, but at the moment this causes a stack overflow somewhere in nimue/keccak
    // merlin.absorb_crs(&crs)?;
    // merlin.ratchet()?;
    // merlin.absorb_instance(&instance)?;
    // merlin.ratchet()?;

    // Prove
    let num_constraints = index.num_constraints;
    let num_ct_constraints = index.num_constant_constraints;

    // Message 1
    // let s_mat = Matrix::<R>::from_columns(&witness.s); // s_mat = [s_1 | ... | s_r] in Rq^{n x r}
    // let t = &crs.A * &s_mat; // t = A * s_mat in Rq^{k x r}
    // let t_decomp = decompose_balanced_polyring(&t, crs.b1, Some(crs.t1)); // t_decomp = [t_1 | ... | t_r] in Rq^{k x r}
    let t: Vec<Vector<R>> = witness
        .s
        .par_iter()
        .map(|s_i| commit(&crs.A, s_i))
        .collect(); // r x k, t[i, j] = (A * s_i)[j]

    let t_decomp = decompose_vec_vector_dimfirst(&t, crs.b1, Some(crs.t1)); // t1 x r x k
    let t_flat = flatten_vec_vec_vector(&t_decomp); // t1 * r * k

    let G = inner_products(&witness.s); // upper triangular matrix in Rq^{r x r}
    let G_decomp = decompose_symmetric_matrix(&G, crs.b2, Some(crs.t2)); // t2 x r x r
    let G_flat = flatten_vec_symmetric_matrix(&G_decomp); // t2 * r * r

    let u_1 = commit(&crs.B, &t_flat) + commit(&crs.C, &G_flat);
    merlin
        .absorb_vector(&u_1)
        .expect("error absorbing prover message 1");

    // Challenge 1
    let num_projections = 256; // TODO: set in CRS
    let Pi = merlin
        .challenge_matrices::<R, WeightedTernaryChallengeSet<R>>(num_projections, crs.n, crs.r)
        .expect("error squeezing verifier message 1"); // r matrices in R^{num_projections x n}

    // Message 2
    let mut p = Vector::<R::BaseRing>::zeros(num_projections);
    for i in 0..crs.r {
        let s_i_vec = R::flattened(&witness.s[i]); // in R::BaseRing^{n*d}
        let pi_i = &Pi[i];
        for j in 0..num_projections {
            let pi_ij = &pi_i.row(j).transpose();
            let pi_ij_vec = R::flattened(pi_ij); // in R::BaseRing^{n*d}
            p[j] += pi_ij_vec.dot(&s_i_vec);
        }
    }
    merlin
        .absorb_vector_canonical::<R::BaseRing>(&p)
        .expect("error absorbing prover message 2");

    // Challenge 2
    let psi = merlin
        .challenge_vectors::<R::BaseRing, R::BaseRing>(num_ct_constraints, crs.num_aggregs)
        .expect("error squeezing verifier message 2 (psi)");
    let omega = merlin
        .challenge_vectors::<R::BaseRing, R::BaseRing>(256, crs.num_aggregs)
        .expect("error squeezing verifier message 2 (omega)");

    // Message 3
    let phi__ = compute_phi__(crs, index, instance, &Pi, &psi, &omega);
    let mut b__ = vec![R::zero(); crs.num_aggregs];
    let a__ = compute_a__(crs, instance, &psi);

    for k in 0..crs.num_aggregs {
        for i in 0..crs.r {
            for j in 0..crs.r {
                b__[k] += a__[k][(i, j)] * G[(i, j)];
            }
            b__[k] += &phi__[k][i].dot(&witness.s[i]);
        }
    }

    merlin
        .absorb_vec(&b__)
        .expect("error absorbing prover message 3");

    // Challenge 3
    let alpha = merlin
        .challenge_vector::<R, R>(num_constraints)
        .expect("error squeezing verifier message 3 (alpha)");
    let beta = merlin
        .challenge_vector::<R, R>(crs.num_aggregs)
        .expect("error squeezing verifier message 3 (beta)");

    // Message 4
    let phi = compute_phi(crs, instance, &alpha, &beta, &phi__);

    let two_inv = R::inverse(&R::try_from(2u64).unwrap()).unwrap();
    let mut H = inner_products2(&phi, &witness.s);
    let H_2 = inner_products2(&witness.s, &phi);
    for i in 0..crs.r {
        for j in 0..i + 1 {
            H[(i, j)] += H_2[(i, j)];
            H[(i, j)] *= two_inv;
        }
    }
    let H_decomp = decompose_symmetric_matrix(&H, crs.b1, Some(crs.t1)); // t1 x r x r
    let H_flat = flatten_vec_symmetric_matrix(&H_decomp); // t1 * r * r

    let u_2 = commit(&crs.D, &H_flat);
    merlin
        .absorb_vector(&u_2)
        .expect("error absorbing prover message 4");

    // Challenge 4
    let c = merlin
        .challenge_vec::<R, LabradorChallengeSet<R>>(crs.r)
        .expect("error squeezing verifier message 4");

    // Compute next instance
    let transcript = TranscriptView {
        u_1,
        b__,
        alpha,
        beta,
        u_2,
        c: c.clone(),
        phi,
        a__,
    };
    let (index_next, instance_next) = fold_instance(&crs, &instance, &transcript);
    let next_size = crs.next_size();

    // Compute next witnes
    let z: Vector<R> = witness
        .clone()
        .s
        .into_iter()
        .zip(c.into_iter())
        .map(|(s_i, c_i)| (s_i * c_i).into())
        .sum();

    let z_decomp = decompose_vec_polyring(&z.as_slice(), crs.b, Some(2usize));

    let mut layouter = Layouter::<R>::new(next_size);

    layouter.set_z0(z_decomp[0].as_slice());
    layouter.set_z1(z_decomp[1].as_slice());
    layouter.set_t(t_flat.as_slice());
    layouter.set_g(G_flat.as_slice());
    layouter.set_h(H_flat.as_slice());

    let witness_next = Witness::<R>::new(layouter.split());

    Ok((index_next, instance_next, witness_next))
}

pub fn prove_principal_relation<'a, R: PolyRing>(
    merlin: &'a mut Merlin,
    mut crs: &CommonReferenceString<R>,
    index: &Index<R>,
    instance: &Instance<R>,
    witness: &Witness<R>,
) -> ProofResult<&'a [u8]>
where
    LabradorChallengeSet<R>: FromRandomBytes<R>,
    WeightedTernaryChallengeSet<R>: FromRandomBytes<R>,
    <R as PolyRing>::BaseRing: WithSignedRepresentative,
    <R::BaseRing as WithSignedRepresentative>::SignedRepresentative:
        DecompositionFriendlySignedRepresentative,
    <R as TryFrom<u128>>::Error: Debug,
{
    let mut index_curr = index.clone();
    let mut instance_curr = instance.clone();
    let mut witness_curr = witness.clone();

    while crs.next_crs.is_some() {
        (index_curr, instance_curr, witness_curr) = prove_principal_relation_oneround(
            merlin,
            crs,
            &index_curr,
            &instance_curr,
            &witness_curr,
        )?;
        crs = crs.next_crs.as_ref().unwrap();
    }
    // TODO: add index/instance to the transcript
    Ok(merlin.transcript())
}
