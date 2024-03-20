use ark_serialize::Compress;
use nimue::{IOPatternError, Merlin, ProofError};
use rayon::iter::{ParallelIterator, IntoParallelRefIterator};

use crate::{check, check_eq};

use crate::lattice_arithmetic::challenge_set::ternary::{TernaryChallengeSet, Trit};
use crate::lattice_arithmetic::matrix::Vector;
use crate::lattice_arithmetic::poly_ring::ConvertibleField;
use crate::lattice_arithmetic::traits::WithL2Norm;
use crate::lova::util::{CRS, SECPARAM};
use crate::nimue::merlin::SerMerlin;
use crate::nimue::traits::ChallengeFromRandomBytes;

pub fn verify_folding<F: ConvertibleField>(merlin: &mut Merlin, crs: &CRS<F>, z: &Vec<Vector<F>>) -> Result<(), ProofError> {
    let f_size = F::zero().serialized_size(Compress::Yes);

    let t_decomp: Vec<Vec<Vector<F>>> =
        (0..2 * SECPARAM).map(|i|
            (0..crs.decomposition_length).map(|j|
                Vector::<F>::from(
                    (0..crs.commitment_mat.nrows()).map(|k|
                        merlin.next_canonical_serializable_size::<F>(f_size).unwrap()
                    ).collect::<Vec<_>>()
                )
            ).collect()
        ).collect();

    let g: Vec<Vec<Vec<Vec<F>>>> =
        (0..(2 * SECPARAM)).map(|i|
            (0..i + 1).map(|j|
                (0..crs.decomposition_length).map(|k1|
                    (0..crs.decomposition_length).map(|k2|
                        merlin.next_canonical_serializable_size::<F>(f_size).unwrap()
                    ).collect()
                ).collect()
            ).collect()
        ).collect();

    let c = merlin.challenge_matrices::<Trit, TernaryChallengeSet<Trit>>(SECPARAM, 2 * SECPARAM, crs.decomposition_length)?;

    // Check Az = sum_{i in [2*SEC]} sum_{j in [k]} c_{i,j} * t_{i,j}
    let z_com: Vec<Vector<F>> = z.iter().map(|z_i|
        &crs.commitment_mat * z_i
    ).collect();

    let ct_lincomb: Vec<Vector<F>> = c.par_iter().map(
        |c_i|
            c_i.row_iter()
                .into_iter()
                .zip(t_decomp.iter())
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
    check_eq!(z_com, ct_lincomb);

    // Check ||z_i|| <= norm_bound_z
    let norm_bound_z = crs.norm_bound * (2 * SECPARAM) as f64; // TODO
    let norm_bound_z_sq = norm_bound_z.powi(2);
    for z_i in z {
        check!(z_i.l2_norm_squared() <= norm_bound_z_sq as u64);
    }

    // Check <z_i, z_i> = sum_{j,k in [2*SEC]} sum_{l,m in [k]} c_{i,j,l} * c_{i,k,m} * g_{j,k,l,m}
    for (i, z_i) in z.iter().enumerate() {
        let z_i_dotprod = z_i.dot(z_i);
        let mut rhs = F::zero();
        for j in 0..2 * SECPARAM {
            for k in 0..2 * SECPARAM {
                for l in 0..crs.decomposition_length {
                    for m in 0..crs.decomposition_length {
                        let c1 = c[i][(j, l)];
                        let c2 = c[i][(k, m)];
                        if c1 != Trit::Zero && c2 != Trit::Zero {
                            if c1 == c2 {
                                rhs += g[j][k][l][m];
                            } else {
                                rhs -= g[j][k][l][m];
                            }
                        }
                    }
                }
            }
        }
        check_eq!(z_i_dotprod, rhs);
    }

    Ok(())
}