#![allow(non_snake_case)]

use lattirust_arithmetic::decomposition::balanced_decomposition::decompose_balanced_vec_polyring;
use lattirust_arithmetic::decomposition::DecompositionFriendlySignedRepresentative;
use lattirust_arithmetic::ring::representatives::WithSignedRepresentative;
use num_traits::{One, Zero};
use tracing::debug;

use lattirust_arithmetic::linear_algebra::{Matrix, Scalar, SymmetricMatrix, Vector};
use lattirust_arithmetic::ring::Ring;
use lattirust_arithmetic::ring::{PolyRing, Z2};

pub fn commit<R: Ring>(A: &Matrix<R>, s: &Vector<R>) -> Vector<R> {
    A * s
}

pub fn flatten_vec_vector<R: Ring>(v: &Vec<Vector<R>>) -> Vec<R> {
    let mut res = Vec::<R>::with_capacity(v.len() * v[0].len());
    for v_i in v {
        res.extend(v_i.as_slice());
    }
    res
}

pub fn flatten_vec_vec_vector<R: Ring>(v: &Vec<Vec<Vector<R>>>) -> Vector<R> {
    let res = v
        .iter()
        .flat_map(|v| flatten_vec_vector(&v))
        .collect::<Vec<_>>();
    Vector::<R>::from_vec(res)
}

pub fn flatten_symmetric_matrix<R: Ring>(v: &SymmetricMatrix<R>) -> Vec<R> {
    v.to_vec()
}

pub fn flatten_vec_symmetric_matrix<R: Ring>(v: &Vec<SymmetricMatrix<R>>) -> Vector<R> {
    let res = v
        .iter()
        .flat_map(|v| flatten_symmetric_matrix(&v))
        .collect::<Vec<_>>();
    Vector::<R>::from_vec(res)
}

pub fn unflatten_symmetric_matrix<R: Ring>(v: &[R], n: usize) -> SymmetricMatrix<R> {
    let mut chunks = Vec::<Vec<R>>::with_capacity(n);
    for i in 0..n {
        let start = (i * (i + 1)) / 2;
        let end = ((i + 1) * (i + 2)) / 2;
        chunks.push(v[start..end].to_vec());
    }
    SymmetricMatrix::<R>::from(chunks)
}

pub fn concat<R: Clone + Scalar>(vecs: &[&[R]]) -> Vector<R> {
    let mut vals = Vec::<R>::with_capacity(vecs.iter().map(|v| v.len()).sum());
    for v in vecs {
        vals.extend_from_slice(v);
    }
    Vector::<R>::from_vec(vals)
}

pub fn shift_right<R: Ring>(v: &[R], shift: usize) -> Vector<R> {
    let mut res = vec![R::zero(); shift];
    res.extend(v);
    Vector::<R>::from_vec(res)
}

pub fn mul_matrix_basescalar<R: PolyRing>(A: &Matrix<R>, x: R::BaseRing) -> Matrix<R> {
    A.map(|a_ij| a_ij * x)
}

pub fn mul_basescalar_vector<R: PolyRing>(s: R::BaseRing, A: &Vector<R>) -> Vector<R> {
    A.map(|a_ij| a_ij * s)
}

/// Compute $\sum_{i,j \in \[r\]} A_ij c_i  c_j$
pub fn linear_combination_symmetric_matrix<R: Ring>(A: &SymmetricMatrix<R>, c: &[R]) -> R {
    let n = A.size();
    debug_assert_eq!(c.len(), n);
    let mut lc = R::zero();
    for i in 0..n {
        for j in 0..n {
            lc += A[(i, j)] * c[i] * c[j];
        }
    }
    lc
}

/// Reinterprets a vector of k = k' * d binary coefficients as k' vectors of d binary coefficients, represented as a vector of k' elements of the polynomial ring R with dimension d.
pub fn lift<R: PolyRing>(vec: &Vector<Z2>) -> Vector<R> {
    let d = R::dimension();
    debug_assert_eq!(
        vec.len() % d,
        0,
        "vector length {} must be multiple of dimension {}",
        vec.len(),
        d
    );
    let coeffs = vec
        .as_slice()
        .chunks(d)
        .map(|chunk| R::from(chunk.iter().copied().map(embed).collect::<Vec<_>>()))
        .collect();
    Vector::<R>::from_vec(coeffs)
}

/// Upcast an element in Z2 to an element in R
pub fn embed<R: Zero + One>(x: Z2) -> R {
    if x.is_zero() {
        R::zero()
    } else {
        debug_assert!(x.is_one());
        R::one()
    }
}

pub fn basis_vector<R: PolyRing>(i: usize, n: usize) -> Vector<R> {
    debug_assert!(i < n, "i = {} must be less than n = {}", i, n);
    let mut coeffs = vec![R::zero(); n];
    coeffs[i] = R::one();
    Vector::<R>::from_vec(coeffs)
}

// In: v: n x n; output: k x n x n where k is the decomposition length
pub fn decompose_symmetric_matrix<PR: PolyRing>(
    mat: &SymmetricMatrix<PR>,
    b: u128,
    padding_size: Option<usize>,
) -> Vec<SymmetricMatrix<PR>>
where
    PR::BaseRing: WithSignedRepresentative,
    <PR::BaseRing as WithSignedRepresentative>::SignedRepresentative:
        DecompositionFriendlySignedRepresentative,
{
    let vec = flatten_symmetric_matrix(mat);
    let decomp = decompose_balanced_vec_polyring(&vec, b, padding_size);
    decomp
        .into_iter()
        .map(|v_i| unflatten_symmetric_matrix(&v_i.as_slice(), mat.size()))
        .collect()
}
