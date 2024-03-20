use ark_std::iterable::Iterable;
use nimue::{Arthur, IOPatternError};
use rayon::prelude::*;

use crate::lattice_arithmetic::balanced_decomposition::decompose_balanced_vec;
use crate::lattice_arithmetic::challenge_set::ternary::{TernaryChallengeSet, Trit};
use crate::lattice_arithmetic::matrix::Vector;
use crate::lattice_arithmetic::poly_ring::ConvertibleField;
use crate::lova::util::{CRS, SECPARAM};
use crate::nimue::arthur::SerArthur;
use crate::nimue::traits::ChallengeFromRandomBytes;

pub fn prove_folding<F: ConvertibleField>(arthur: &mut Arthur, crs: &CRS<F>, s: &Vec<Vector<F>>) -> Result<(Vec<Vector<F>>), IOPatternError> {
    debug_assert_eq!(s.len(), 2 * SECPARAM);

    // s_decomp in F^(2*SEC x k x n)
    let mut s_decomp: Vec<Vec<Vector<F>>> = s.par_iter().map(|s_i|
        decompose_balanced_vec(s_i.as_slice(), crs.decomposition_basis, Some(crs.decomposition_length)).into_iter().map(|v| Vector::<F>::from_vec(v)).collect()
    ).collect();

    // debug_assert_eq!(s_decomp.len(), 2*SECPARAM);
    // debug_assert!(s_decomp.iter().all(|v| v.len() == k));
    // debug_assert!(s_decomp.iter().all(|v| v.iter().all(|w| w.len() == n)));

    // t_decomp in F^(2*SEC x k x m)
    let t_decomp = s_decomp.par_iter().map(|s_i|
        s_i.iter().map(|s_ij|
            &crs.commitment_mat * s_ij
        ).collect::<Vec<_>>()
    ).collect::<Vec<_>>();

    // Add t to FS transcript
    for t_i in t_decomp.iter() {
        for t_ij in t_i.iter() {
            for t_ijk in t_ij.iter() {
                arthur.absorb_canonical_serializable(t_ijk)?;
            }
        }
    }

    // g in F^(2*SEC x 2*SEC x k x k)
    let g: Vec<Vec<Vec<Vec<F>>>> = (0..(2 * SECPARAM)).into_par_iter().map(|i| {
        (0..i + 1).into_par_iter().map(|j| {
            (0..crs.decomposition_length).map(|k1|
                (0..crs.decomposition_length).map(|k2|
                    s_decomp[i][k1].dot(&s_decomp[j][k2])
                ).collect()
            ).collect()
        }).collect()
    }).collect();

    // Add g to FS transcript
    for g_i in g.iter() {
        for g_ij in g_i.iter() {
            for g_ijk in g_ij.iter() {
                for g_ijkl in g_ijk.iter() {
                    arthur.absorb_canonical_serializable(g_ijkl)?;
                }
            }
        }
    }

    // Get challenge
    let c = arthur.challenge_matrices::<Trit, TernaryChallengeSet<Trit>>(SECPARAM, 2 * SECPARAM, crs.decomposition_length)?;

    // Compute z
    let z: Vec<Vector<F>> = c.par_iter().map(
        |c_i|
            c_i.row_iter()
                .into_iter()
                .zip(s_decomp.iter())
                .map(
                    |(c_ij, s_i)|
                        c_ij.into_iter()
                            .zip(s_i)
                            .filter(|(c_ijk, _)| **c_ijk != Trit::Zero)
                            .map(|(c_ijk, s_ij)|
                                // s_ij * c_ijk
                                match c_ijk {
                                    Trit::MinusOne => -s_ij,
                                    Trit::One => s_ij.clone(),
                                    _ => unreachable!()
                                }
                            )
                            .fold(Vector::<F>::zeros(crs.n), |acc, x| acc + x) // we can't use sum(), as the iterator might be empty
                )
                .fold(Vector::<F>::zeros(crs.n), |acc, x| acc + x) // here it could be empty with negligible probability, but let's make sure tests don't panic
    ).collect();

    Ok(z)
}