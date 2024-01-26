#![allow(non_snake_case)]
use ark_ff::Field;
use num_traits::{One, Zero};
use crate::labrador::binary_r1cs::prover_binary_r1cs::{F2, Z2};
use crate::lattice_arithmetic::matrix::Vector;
use crate::lattice_arithmetic::poly_ring::PolyRing;

fn inner_prod<F: Field>(row: &[(F, usize)], witness: &[F]) -> F {
    let mut acc = F::zero();
    for &(ref coeff, i) in row {
        acc += &(if coeff.is_one() { witness[i] } else { witness[i] * coeff });
    }
    acc
}

pub fn mul_sparse_mat_vec<R: Field>(matrix: &Vec<Vec<(R, usize)>>, vec: &Vec<R>) -> Vec<R> {
    ark_std::cfg_iter!(matrix)
        .map(|row| inner_prod(row, vec))
        .collect::<Vec<_>>()
}

pub fn embed<R: PolyRing>(vec: &Vector<Z2>) -> Vector<R> {
    Vector::<R>::from_fn(vec.len(), |i, _| if vec[i].is_zero() { R::zero() } else { R::one() })
}

// pub fn to_sparse_matrix(mat: &ark_relations::r1cs::Matrix<F2>) -> nalgebra_sparse::CsrMatrix<F2> {
//     let mut coo_mat = nalgebra_sparse::CooMatrix::<F2>::new(mat.len(), mat[0].len());
//     for (i, row) in mat.iter().enumerate() {
//         for (coeff, j) in row.into_iter() {
//             if !coeff.is_zero() {
//                 coo_mat.push(i, *j, *coeff);
//             }
//         }
//     }
//     nalgebra_sparse::CsrMatrix::<F2>::from(&coo_mat)
// }

pub fn to_Z2_vec(vec: &Vec<F2>) -> Vector<Z2> {
    Vector::<Z2>::from_vec(vec.iter().map(|x| if x.is_zero() { Z2::zero() } else { Z2::one() }).collect())
}

pub fn to_vector(vec: &Vec<F2>) -> Vector<Z2> {
    Vector::<Z2>::from_vec(vec.iter().map(|x| if x.is_zero() { Z2::zero() } else { Z2::one() }).collect())
}
