use std::iter::Sum;
use std::ops::{Add, Div, Mul, Rem, Sub};

use nalgebra::Scalar;
use num_traits::{Signed, Zero};
use rayon::iter::{
    IndexedParallelIterator, IntoParallelIterator, IntoParallelRefIterator, ParallelIterator,
};
use rounded_div::RoundedDiv;

use crate::linear_algebra::{
    ClosedAddAssign, ClosedMulAssign, Matrix, RowVector, SymmetricMatrix, Vector,
};
use crate::ring::Ring;

pub mod balanced_decomposition;
#[allow(clippy::module_inception)]
pub mod decomposition;

pub trait DecompositionFriendlySignedRepresentative:
    Clone
    + Zero
    + Signed
    + RoundedDiv<Output = Self>
    + Div<Self, Output = Self>
    + Rem<Self, Output = Self>
    + Add<Self, Output = Self>
    + Sub<Self, Output = Self>
    + PartialOrd
{
}

impl<T> DecompositionFriendlySignedRepresentative for T where
    T: Clone
        + Zero
        + Signed
        + RoundedDiv<Output = Self>
        + Div<Self, Output = Self>
        + Rem<Self, Output = Self>
        + Add<Self, Output = Self>
        + Sub<Self, Output = Self>
        + PartialOrd
{
}

pub(crate) fn pad_zeros<R: Zero + Clone>(v: &mut Vec<R>, padding_size: Option<usize>) {
    if let Some(padding_size) = padding_size {
        if v.len() > padding_size {
            // Strip trailing zeros
            while v.last().is_some_and(|last| last.is_zero()) {
                v.pop();
            }
        }
        assert!(
            v.len() <= padding_size,
            "padding_size = {} must be at least v.len() = {}",
            padding_size,
            v.len(),
        );
        v.resize(padding_size, R::zero());
    } else {
        // Strip trailing zeros
        while v.last().is_some_and(|last| last.is_zero()) {
            v.pop();
        }
    }
}

/// Given a vector of vectors `v`, pads each row to the same length $l = \max_i \texttt{v}\[i\]\texttt{.len()}$ and transposes the result. The output is a Vec of Vec of dimensions `l` times `v.len()`.
pub fn pad_and_transpose<F: Copy + Zero>(
    mut v: Vec<Vec<F>>,
    padding_size: Option<usize>,
) -> Vec<Vec<F>> {
    if v.is_empty() {
        return vec![];
    }
    let rows = v.len();
    let max_cols = v.iter().map(|d_i| d_i.len()).max().unwrap();
    let cols = match padding_size {
        None => max_cols,
        Some(padding_length) => {
            assert!(padding_length >= max_cols);
            padding_length
        }
    };

    // Pad each row to the same length `cols`
    for row in &mut v {
        row.resize(cols, F::zero());
    }

    // Reshape as cols x rows
    (0..cols)
        .map(|col| (0..rows).map(|row| v[row][col]).collect::<Vec<F>>())
        .collect()
}

pub fn pad_and_transpose_vec_vec_vector<F: Copy + Scalar + Zero>(
    mut v: Vec<Vec<Vector<F>>>,
    padding_size: Option<usize>,
) -> Vec<Vec<Vector<F>>> {
    if v.is_empty() {
        return vec![];
    }
    let rows = v.len();
    let max_cols = v.iter().map(|d_i| d_i.len()).max().unwrap();
    let vector_len = v[0][0].len();

    let cols = match padding_size {
        None => max_cols,
        Some(padding_length) => {
            assert!(padding_length >= max_cols);
            padding_length
        }
    };

    // Pad each row to the same length `cols`

    for row in &mut v {
        row.resize(cols, Vector::<F>::zeros(vector_len));
    }

    // Reshape as cols x rows
    let res: Vec<Vec<Vector<F>>> = (0..cols)
        .map(|col| {
            (0..rows)
                .map(|row| v[row][col].clone())
                .collect::<Vec<Vector<F>>>()
        })
        .collect();
    debug_assert_eq!(res.len(), cols);
    debug_assert!(res.iter().all(|r| r.len() == rows));
    debug_assert!(res
        .iter()
        .all(|r| { r.iter().all(|v| v.len() == vector_len) }));

    res
}

pub fn recompose<A, B>(vec: &Vec<A>, basis: B) -> A
where
    A: Mul<B, Output = A> + Copy + Sum + Send + Sync,
    B: Ring,
{
    vec.par_iter()
        .enumerate()
        .map(|(i, v_i)| *v_i * basis.pow([i as u64]))
        .sum()
}

/// Given a `m x n*k` matrix `mat` decomposed in basis b and a slice \[1, b, ..., b^(k-1)] `powers_of_basis`, returns the `m x n` recomposed matrix.
pub fn recompose_matrix<F: Scalar + ClosedAddAssign + ClosedMulAssign + Zero + Send + Sync>(
    mat: &Matrix<F>,
    powers_of_basis: &[F],
) -> Matrix<F> {
    let nk = mat.ncols();
    let k = powers_of_basis.len();
    debug_assert!(nk % k == 0, "matrix `mat` to be recomposed should have dimensions m x nk, where k is the length of `powers_of_basis`, but its number of columns is not a multiple of k");
    let n = nk / k;
    let pows = Vector::<F>::from_slice(powers_of_basis).transpose();
    let rows = mat
        .par_row_iter()
        .map(|r_i| {
            RowVector::from(
                (0..n)
                    .into_par_iter()
                    .map(|j| r_i.columns_range(j * k..j * k + k).dot(&pows))
                    .collect::<Vec<_>>(),
            )
        })
        .collect::<Vec<_>>();
    Matrix::<F>::from_rows(rows.as_slice())
}

/// Given a `n*d x n*d` symmetric matrix `mat` and a slice `\[1, b, ..., b^(d-1)\]` `powers_of_basis`, returns the `n x n` symmetric matrix corresponding to $G^T \textt{mat} G$, where $G = I_n \otimes (1, b, ..., b^(\textt{d}-1))$ is the gadget matrix of dimensions `n*d x n`.
pub fn recompose_left_right_symmetric_matrix<F: Clone + Sum + Send + Sync>(
    mat: &SymmetricMatrix<F>,
    powers_of_basis: &[F],
) -> SymmetricMatrix<F>
where
    for<'a> &'a F: Mul<&'a F, Output = F>,
{
    let (nd, d) = (mat.size(), powers_of_basis.len());
    assert_eq!(nd % d, 0);

    let n = nd / d;
    (0..n)
        .into_par_iter()
        .map(|i| {
            (0..=i)
                .into_par_iter()
                .map(|j| {
                    (0..nd)
                        .filter(|k| k / d == i)
                        .flat_map(|k| {
                            (0..nd).filter(|l| l / d == j).map(move |l| {
                                &mat[(k, l)] * &(&powers_of_basis[k % d] * &powers_of_basis[l % d])
                            })
                        })
                        .sum()
                })
                .collect()
        })
        .collect::<Vec<_>>()
        .into()
}
