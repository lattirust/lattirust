use nimue::InvalidTag;
use num_traits::Zero;
use crate::lattice_arithmetic::balanced_decomposition::decompose_balanced_vec;
use crate::lattice_arithmetic::challenge_set::labrador_challenge_set::LabradorChallengeSet;
use crate::lattice_arithmetic::challenge_set::weighted_ternary::WeightedTernaryChallengeSet;
use crate::lattice_arithmetic::matrix::{norm_sq_ringelem, norm_vec_basering, Vector};
use crate::lattice_arithmetic::poly_ring::PolyRing;
use crate::lattice_arithmetic::traits::FromRandomBytes;
use crate::nimue::merlin::LatticeMerlin;
use crate::relations::labrador::principal_relation::PrincipalRelation;
use crate::labrador::common_reference_string::CommonReferenceString;

#[macro_export]
macro_rules! check {
    ($ cond: expr) => {
        {
            if !($cond) {
                return Err(InvalidTag::from("invalid proof"));
            }
        }
    };
    ( $ cond : expr , $ ( $ arg : tt ) + ) => {
        {
            if !($cond) {
                return Err(InvalidTag::from(format!("invalid proof: {}", format!($($arg)+))));
            }
        }
    };
}

#[macro_export]
macro_rules! check_eq {
    ( $ a : expr , $ b : expr ) => { check!($a == $b) };
    ( $ a : expr , $ b : expr , $ ( $ arg : tt ) + ) => { check!($a == $b, $($arg)+) };
}



pub fn verify_principal_relation<R: PolyRing>(merlin: &mut LatticeMerlin, instance: PrincipalRelation<R>, crs: CommonReferenceString<R>) -> Result<(), InvalidTag>
    where LabradorChallengeSet<R>: FromRandomBytes<R>, WeightedTernaryChallengeSet<R>: FromRandomBytes<R>, i128: From<<R as PolyRing>::BaseRing>, i64: From<<R as PolyRing>::BaseRing>
{
    let (n, r) = (crs.n, crs.r);
    let num_constraints = instance.quad_dot_prod_funcs.len();
    let num_ct_constraints = instance.ct_quad_dot_prod_funcs.len();

    let u_1 = merlin.next_vector::<R>(crs.k1).expect("error extracting prover message 1 from transcript");

    let pis = merlin.challenge_matrices::<R, WeightedTernaryChallengeSet<R>>(256, n, r).expect("error extracting verifier message 1 from transcript");

    let p = merlin.next_vector::<R::BaseRing>(256).expect("error extracting prover message 2 from transcript");
    let norm_p = norm_vec_basering::<R>(&p);
    let p_norm_bound = 128f64.sqrt() * instance.norm_bound;
    check!(norm_p <=p_norm_bound, "||p||_2 = {} is not <= sqrt(128)*beta = {}", norm_p, p_norm_bound);

    let psi = merlin.challenge_vectors::<R::BaseRing, R::BaseRing>(num_ct_constraints, crs.num_aggregs).expect("error extracting verifier message 2 (psi) from transcript");
    let omega = merlin.challenge_vectors::<R::BaseRing, R::BaseRing>(256, crs.num_aggregs).expect("error extracting verifier message 2 (omega) from transcript");

    let b__ = merlin.next_vec::<R>(crs.num_aggregs).expect("error extracting prover message 3 from transcript");

    for k in 0..crs.num_aggregs {
        let mut rhs_k = omega[k].dot(&p);
        for l in 0..num_ct_constraints {
            rhs_k += psi[k][l] * instance.ct_quad_dot_prod_funcs[l].b.coeffs()[0];
        }
        check_eq!(b__[k].coeffs()[0], rhs_k, "constant coeff of b''^(k) = {:?} is not equal to rhs = {:?}", b__[k].coeffs()[0], rhs_k);
    }

    let alpha = merlin.challenge_vector::<R, R>(num_constraints).expect("error extracting verifier message 3 (alpha) from transcript");
    let beta = merlin.challenge_vector::<R, R>(crs.num_aggregs).expect("error extracting verifier message 3 (beta) from transcript");

    let u_2 = merlin.next_vector::<R>(crs.k2).expect("error extracting prover message 4 from transcript");

    let c = merlin.challenge_vec::<R, LabradorChallengeSet<R>>(crs.r).expect("error extracting verifier message 4 from transcript");

    let z = merlin.next_vector::<R>(crs.n).expect("error extracting prover message 5 (z) from transcript");
    let t = merlin.next_vectors::<R>(crs.k, crs.t1).expect("error extracting prover message 5 (t) from transcript");
    let G = merlin.next_symmetric_matrix::<R>(crs.r).expect("error extracting prover message 5 (G) from transcript");
    let H = merlin.next_symmetric_matrix::<R>(crs.r).expect("error extracting prover message 5 (H) from transcript");

    // Bulk of the verification checks
    // TODO: compute a__ijk, etc.

    // We don't need to check if G and H are symmetric, since we only send the upper triangular part of the matrices

    // Decompose z
    let mut z_decomp = decompose_balanced_vec(&z, crs.decomposition_basis, Some(2usize));
    assert_eq!(z_decomp.len(), 2);
    check!(l_inf_norm(&z_decomp[0]) * 2 <= i64::from(crs.decomposition_basis) as u64);

    // Decompose t_i
    for t_i in t {
        let mut t_decomp = decompose_balanced_vec(&t_i, crs.decomposition_basis, Some(crs.t1));
        assert_eq!(t_decomp.len(), crs.t1);
        // check!();
    }

    Ok(())
}

fn l_inf_norm<R: PolyRing>(v: &Vector<R>) -> u64
    where i64: From<R::BaseRing>
{
    R::flattened_coeffs(v).into_iter().map(|x| i64::from(x).abs() as u64).max().unwrap()
}