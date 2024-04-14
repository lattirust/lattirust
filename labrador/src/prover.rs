#![allow(non_snake_case)]

use nimue::{Arthur, ProofResult};
use rayon::prelude::*;

use lattirust_arithmetic::balanced_decomposition::{decompose_balanced_polyring, decompose_balanced_vec_polyring};
use lattirust_arithmetic::challenge_set::labrador_challenge_set::LabradorChallengeSet;
use lattirust_arithmetic::challenge_set::weighted_ternary::WeightedTernaryChallengeSet;
use lattirust_arithmetic::linear_algebra::inner_products::{inner_products, inner_products2};
use lattirust_arithmetic::nimue::arthur::SerArthur;
use lattirust_arithmetic::nimue::traits::ChallengeFromRandomBytes;
use lattirust_arithmetic::poly_ring::{PolyRing, SignedRepresentative};
use lattirust_arithmetic::traits::FromRandomBytes;
use lattirust_arithmetic::linear_algebra::Vector;
use relations::principal_relation::{PrincipalRelation, Witness};

use crate::common_reference_string::{CommonReferenceString, fold_instance};
use crate::shared::BaseTranscript;
use crate::util::*;

pub struct Prover<'a, R: PolyRing> {
    pub transcript: BaseTranscript<'a, R>,
    pub witness: &'a Witness<R>,
    pub(crate) phi__: Option<Vec<Vec<Vector<R>>>>,
}

impl<'a, R: PolyRing> Prover<'a, R> {
    fn new(instance: &'a PrincipalRelation<R>, witness: &'a Witness<R>, crs: &'a CommonReferenceString<R>) -> Self {
        Self {
            transcript: BaseTranscript::<R>::init(&crs, &instance),
            witness: &witness,
            phi__: None,
        }
    }

    fn instance(&self) -> &PrincipalRelation<R> { &self.transcript.instance }
    fn crs(&self) -> &CommonReferenceString<R> { &self.transcript.crs }
    fn witness(&self) -> &Witness<R> { &self.witness }

    fn prove_1(&mut self) {
        let (witness, crs) = (self.witness, self.crs());
        let t: Vec<Vector<R>> = witness.s.par_iter().map(
            |s_i| commit(&crs.A, s_i)
        ).collect();
        let t_decomp: Vec<Vec<Vector<R>>> = t.par_iter().map(|t_i| decompose_balanced_vec_polyring(t_i, crs.b1, Some(crs.t1))).collect();

        let G = inner_products(&witness.s);

        let G_decomp: Vec<Vec<Vec<R>>> = G.rows().par_iter().map(
            |G_i| G_i.par_iter().map(
                |G_ij| decompose_balanced_polyring(G_ij, crs.b2, Some(crs.t2))
            ).collect()
        ).collect();

        let mut u_1 = Vector::<R>::zeros(crs.k1);
        for i in 0..crs.r {
            assert_eq!(t_decomp[i].len(), crs.t1, "decomposition of t has the wrong number of elements");
            for k in 0..crs.t1 {
                u_1 += &crs.B[i][k] * &t_decomp[i][k];
            }

            for j in 0..i + 1 {
                assert_eq!(G_decomp[i][j].len(), crs.t2, "decomposition of G has the wrong number of elements");
                for k in 0..crs.t2 {
                    u_1 += &crs.C[i][j][k] * G_decomp[i][j][k];
                }
            }
        }

        self.transcript.u_1.replace(u_1);
        self.transcript.t.replace(t);
        self.transcript.G.replace(G);
    }

    fn prove_2(&mut self) {
        let (witness, crs) = (self.witness(), self.crs());
        let Pi = self.transcript.Pi.as_ref().expect("Pi not available");
        let mut p = Vector::<R::BaseRing>::zeros(256);
        for j in 0..256 {
            for i in 0..crs.r {
                let pi_ij = &Pi[i].row(j).transpose();
                p[j] += R::flattened(&pi_ij).dot(&R::flattened(&witness.s[i]));
            }
        }
        self.transcript.p.replace(p);
    }

    fn prove_3(&mut self) {
        let (instance, witness, crs) = (self.instance(), self.witness(), &self.crs());
        let Pi = self.transcript.Pi.as_ref().expect("Pi not available");
        let psi = self.transcript.psi.as_ref().expect("psi not available");
        let omega = self.transcript.omega.as_ref().expect("omega not available");
        debug_assert_eq!(Pi.len(), crs.r);
        debug_assert_eq!(psi.len(), crs.num_aggregs);
        debug_assert_eq!(omega.len(), crs.num_aggregs);
        for Pi_i in Pi {
            debug_assert_eq!(Pi_i.nrows(), 256);
            debug_assert_eq!(Pi_i.ncols(), crs.n);
        }
        for psi_i in psi { debug_assert_eq!(psi_i.len(), instance.ct_quad_dot_prod_funcs.len()); }
        for omega_i in omega { debug_assert_eq!(omega_i.nrows(), 256); }

        let mut phi__ = vec![vec![Vector::<R>::zeros(crs.n); crs.r]; crs.num_aggregs];
        let mut b__ = vec![R::zero(); crs.num_aggregs];

        for k in 0..crs.num_aggregs {
            for i in 0..crs.r {
                for j in 0..crs.r {
                    // Compute a''_{ij}^{(k)}
                    let mut a_ijk = R::zero();
                    for l in 0..instance.ct_quad_dot_prod_funcs.len() {
                        if let Some(A_l) = &instance.ct_quad_dot_prod_funcs[l].A {
                            a_ijk += A_l[(i, j)] * psi[k][l];
                        }
                    }

                    b__[k] += a_ijk * &witness.s[i].dot(&witness.s[j]);
                }
                // Compute vec{phi}''_i^{(k)}
                for l in 0..instance.ct_quad_dot_prod_funcs.len() {
                    phi__[k][i] += mul_basescalar_vector(psi[k][l], &instance.ct_quad_dot_prod_funcs[l].phi[i]);
                }
                for j in 0..256 {
                    phi__[k][i] += mul_basescalar_vector(omega[k][j], &R::sigma_vec(&Pi[i].row(j).transpose()));
                }
                b__[k] += &phi__[k][i].dot(&witness.s[i]);
            }
        }
        self.phi__.replace(phi__);
        self.transcript.b__.replace(b__);
    }

    fn prove_4(&mut self) {
        let (instance, witness, crs) = (self.instance(), self.witness(), self.crs());
        let alpha = self.transcript.alpha.as_ref().expect("alpha not available");
        let beta = self.transcript.beta.as_ref().expect("beta not available");
        let phi__ = self.phi__.as_ref().expect("phi'' not available");

        // Compute phi (in parallel)
        let phi = (0..crs.r).into_par_iter().map(|i| {
            let mut phi_i = Vector::<R>::zeros(crs.n);
            for k in 0..instance.quad_dot_prod_funcs.len() {
                phi_i += &instance.ct_quad_dot_prod_funcs[k].phi[i] * alpha[k];
            }
            for k in 0..crs.num_aggregs {
                phi_i += &phi__[k][i] * beta[k];
            }
            phi_i
        }).collect::<Vec<_>>();

        // Compute H (in parallel)
        let mut H = inner_products2(&phi, &witness.s);
        let H_2 = inner_products2(&witness.s, &phi);
        for i in 0..crs.r {
            for j in 0..i + 1 {
                H[(i, j)] += H_2[(i, j)];
                H[(i, j)] = R::from(
                    H[(i, j)].coeffs().into_iter().map(|h|
                        Into::<R::BaseRing>::into(SignedRepresentative(Into::<SignedRepresentative>::into(h).0 / 2))
                    ).collect::<Vec<_>>()
                );
            }
        }

        let mut u_2 = Vector::<R>::zeros(crs.k2);
        for i in 0..crs.r {
            for j in 0..i + 1 {
                let h_ij_decomp = decompose_balanced_polyring(&H[(i, j)], crs.b, Some(crs.t1));
                for k in 0..crs.t1 {
                    u_2 += &crs.D[i][j][k] * h_ij_decomp[k];
                }
            }
        }

        self.transcript.u_2.replace(u_2);
        self.transcript.phi.replace(phi);
        self.transcript.H.replace(H);
    }

    fn prove_5(&mut self) {
        let (witness, crs) = (self.witness(), self.crs());
        let c = self.transcript.c.as_ref().expect("c not available");
        let mut z = Vector::<R>::zeros(crs.n);
        debug_assert_eq!(c.len(), crs.r);
        for i in 0..crs.r {
            z += &witness.s[i] * c[i];
        }
        self.transcript.z.replace(z);
    }
}


pub fn prove_principal_relation<'a, R: PolyRing>(arthur: &'a mut Arthur, instance: &PrincipalRelation<R>, witness: &Witness<R>, crs: &CommonReferenceString<R>) -> ProofResult<&'a [u8]>
    where LabradorChallengeSet<R>: FromRandomBytes<R>, WeightedTernaryChallengeSet<R>: FromRandomBytes<R>
{
    // Check dimensions and norms
    debug_assert!(crs.is_wellformed_instance(&instance));
    debug_assert!(crs.is_wellformed_witness(&witness));

    // Check that the witness is valid for this instance
    debug_assert!(instance.is_valid_witness(&witness));

    // Initialize Fiat-Shamir transcript
    // TODO: add the following, but at the moment this causes a stack overflow somewhere in nimue/keccak
    // arthur.absorb_crs(&crs)?;
    // arthur.ratchet()?;
    // arthur.absorb_instance(&instance)?;
    // arthur.ratchet()?;

    // Prove
    let num_constraints = instance.quad_dot_prod_funcs.len();
    let num_ct_constraints = instance.ct_quad_dot_prod_funcs.len();
    let mut prover = Prover::new(&instance, &witness, &crs);

    prover.prove_1();
    let u_1 = prover.transcript.u_1.as_ref().unwrap();
    arthur.absorb_vector(&u_1).expect("error absorbing prover message 1");

    let pis = arthur.challenge_matrices::<R, WeightedTernaryChallengeSet<R>>(256, crs.n, crs.r).expect("error squeezing verifier message 1");
    prover.transcript.Pi.replace(pis);

    prover.prove_2();
    let p = prover.transcript.p.as_ref().unwrap();
    arthur.absorb_vector_canonical::<R::BaseRing>(&p).expect("error absorbing prover message 2");

    let psi = arthur.challenge_vectors::<R::BaseRing, R::BaseRing>(num_ct_constraints, crs.num_aggregs).expect("error squeezing verifier message 2 (psi)");
    let omega = arthur.challenge_vectors::<R::BaseRing, R::BaseRing>(256, crs.num_aggregs).expect("error squeezing verifier message 2 (omega)");
    prover.transcript.psi.replace(psi);
    prover.transcript.omega.replace(omega);

    prover.prove_3();
    let b__ = prover.transcript.b__.as_ref().unwrap();
    arthur.absorb_vec(&b__).expect("error absorbing prover message 3");

    let alpha = arthur.challenge_vector::<R, R>(num_constraints).expect("error squeezing verifier message 3 (alpha)");
    let beta = arthur.challenge_vector::<R, R>(crs.num_aggregs).expect("error squeezing verifier message 3 (beta)");
    prover.transcript.alpha.replace(alpha);
    prover.transcript.beta.replace(beta);

    prover.prove_4();
    let u_2 = prover.transcript.u_2.as_ref().unwrap();
    arthur.absorb_vector(&u_2).expect("error absorbing prover message 4");

    let c = arthur.challenge_vec::<R, LabradorChallengeSet<R>>(crs.r).expect("error squeezing verifier message 4");
    prover.transcript.c.replace(c);

    prover.prove_5();

    let recurse = crs.next_crs.is_some();
    if recurse {
        // Fold relation and recurse
        let (instance_next, witness_next) = fold_instance(&prover.transcript, true);
        prove_principal_relation(arthur, &instance_next, &witness_next.unwrap(), &crs.next_crs.as_ref().unwrap())
    } else {
        // Append final messages to Fiat-Shamir transcript and finish
        // TODO: implement optimization of sending G, H first
        let z = prover.transcript.z.as_ref().unwrap();
        let t = prover.transcript.t.as_ref().unwrap();
        let G = prover.transcript.G.as_ref().unwrap();
        let H = prover.transcript.H.as_ref().unwrap();

        arthur.absorb_vector(&z).expect("error absorbing prover message 5 (z)");
        arthur.absorb_vectors(&t).expect("error absorbing prover message 5 (t)");
        arthur.absorb_symmetric_matrix(&G).expect("error absorbing prover message 5 (G)");
        arthur.absorb_symmetric_matrix(&H).expect("error absorbing prover message 5 (H)");

        Ok(arthur.transcript())
    }
}