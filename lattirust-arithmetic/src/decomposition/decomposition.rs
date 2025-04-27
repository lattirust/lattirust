use num_bigint::BigUint;
use num_traits::{One, ToPrimitive, Zero};
use rayon::prelude::*;

use crate::decomposition::pad_zeros;
use crate::linear_algebra::{Matrix, RowVector, Vector};
use crate::ring::representatives::WithSignedRepresentative;
use crate::ring::{PolyRing, Ring};

use super::{pad_and_transpose, DecompositionFriendlySignedRepresentative};

/// Returns the maximum number of terms in the decomposition in basis `b` of any `x` with $|\textt{x}| \leq \textt{max}$.
pub fn decomposition_max_length(b: u128, max: BigUint) -> usize {
    if max.is_zero() {
        0
    } else {
        max.to_f64().unwrap().log(b as f64).floor() as usize + 1
    }
}

/// Returns the decomposition of a slice as a Vec of Vecs.
///
/// # Arguments
/// * `v`: input element
/// * `b`: basis for the decomposition, must be even
/// * `padding_size`: indicates whether the output should be padded with zeros to a specified length `k` if `padding_size` is `Some(k)`, or if it should be padded to the largest decomposition length required for `v` if `padding_size` is `None`
///
/// # Output
/// Returns `d`, the decomposition in basis `basis` as a Vec of size `decomp_size`, i.e.,
/// $\texttt{v}\[i\] = \sum_{j \in \[k\]} \texttt{basis}^j \texttt{d}\[j\]$ and $|\texttt{d}\[j\]| \le \texttt{basis}$.
pub fn decompose<R>(v: &R, basis: u128, padding_size: Option<usize>) -> Vec<R>
where
    R: Ring + WithSignedRepresentative,
    R::SignedRepresentative: DecompositionFriendlySignedRepresentative,
{
    assert!(
        !basis.is_zero() && !basis.is_one(),
        "cannot decompose in basis 0 or 1"
    );
    // TODO: not sure if this really necessary, but having b be even allow for more efficient divisions/remainders
    assert_eq!(basis % 2, 0, "decomposition basis must be even");

    let mut decomp = Vec::<R>::new();
    let mut curr = Into::<R::SignedRepresentative>::into(*v);
    let b = R::SignedRepresentative::try_from(basis as i128).unwrap();
    loop {
        let rem = curr.clone() % b.clone(); // rem = curr % b is in [-(b-1), (b-1)]
        decomp.push(R::try_from(rem).unwrap());
        curr = curr.clone() / b.clone(); // Rust integer division rounds towards zero
        if curr.is_zero() {
            break;
        }
    }

    pad_zeros(&mut decomp, padding_size);

    decomp
}

/// Returns the decomposition of a slice as a Vec of Vecs.
///
/// # Arguments
/// * `v`: input slice, of length `l`
/// * `b`: basis for the decomposition, must be even
/// * `padding_size`: indicates whether the output should be padded with zeros to a specified length `k` if `padding_size` is `Some(k)`, or if it should be padded to the largest decomposition length required for `v` if `padding_size` is `None`
///
/// # Output
/// Returns `d` the decomposition in basis `b` as a Vec of size `decomp_size`, with each item being a Vec of length `l`, i.e.,
/// for all $i \in \[l\]: \texttt{v}\[i\] = \sum_{j \in \[k\]} \texttt{b}^j \texttt{d}\[i\]\[j\]$ and $|\texttt{d}\[i\]\[j\]| \le \frac{\texttt{b}$.
pub fn decompose_vec<R>(v: &[R], b: u128, padding_size: Option<usize>) -> Vec<Vec<R>>
where
    R: Ring + WithSignedRepresentative,
    R::SignedRepresentative: DecompositionFriendlySignedRepresentative,
{
    let decomp: Vec<Vec<R>> = v
        .par_iter()
        .map(|v_i| decompose(v_i, b, padding_size))
        .collect(); // v.len() x decomp_size
    pad_and_transpose(decomp, padding_size) // decomp_size x v.len()
}

/// Returns the decomposition of a [`PolyRing`] element as a Vec of [`PolyRing`] elements.
///
/// # Arguments
/// * `v`: `PolyRing` element to be decomposed
/// * `b`: basis for the decomposition, must be even
/// * `padding_size`: indicates whether the output should be padded with zeros to a specified length `k` if `padding_size` is `Some(k)`, or if it should be padded to the largest decomposition length required for `v` if `padding_size` is `None`
///
/// # Output
/// Returns `d` the decomposition in basis `b` as a Vec of size `decomp_size`, i.e.,
/// for all $\texttt{v} = \sum_{j \in \[k\]} \texttt{b}^j \texttt{d}\[j\]$ and $|\texttt{d}\[j\]| \le \texttt{b}$.
pub fn decompose_polyring<PR: PolyRing>(v: &PR, b: u128, padding_size: Option<usize>) -> Vec<PR>
where
    PR::BaseRing: WithSignedRepresentative,
    <PR::BaseRing as WithSignedRepresentative>::SignedRepresentative:
        DecompositionFriendlySignedRepresentative,
{
    decompose_vec::<PR::BaseRing>(v.coefficients().as_slice(), b, padding_size)
        .into_par_iter()
        .map(|v_i| PR::from(v_i))
        .collect()
}

/// Returns the decomposition of a slice of [`PolyRing`] elements as a Vec of [`Vector`] of [`PolyRing`] elements.
///
/// # Arguments
/// * `v`: input slice, of length `l`
/// * `b`: basis for the decomposition, must be even
/// * `padding_size`: indicates whether the output should be padded with zeros to a specified length `k` if `padding_size` is `Some(k)`, or if it should be padded to the largest decomposition length required for `v` if `padding_size` is `None`
///
/// # Output
/// Returns `d` the decomposition in basis `b` as a Vec of size `decomp_size`, with each item being a Vec of length `l`, i.e.,
/// for all $i \in \[l\]: \texttt{v}\[i\] = \sum_{j \in \[k\]} \texttt{b}^j \texttt{d}\[i\]\[j\]$ and $|\texttt{d}\[i\]\[j\]| \le \texttt{b}$.
pub fn decompose_vec_polyring<PR: PolyRing>(
    v: &[PR],
    b: u128,
    padding_size: Option<usize>,
) -> Vec<Vector<PR>>
where
    PR::BaseRing: WithSignedRepresentative,
    <PR::BaseRing as WithSignedRepresentative>::SignedRepresentative:
        DecompositionFriendlySignedRepresentative,
{
    let decomp: Vec<Vec<PR>> = v
        .par_iter()
        .map(|ring_elem| decompose_polyring(ring_elem, b, padding_size))
        .collect(); // v.len() x decomp_size
    pad_and_transpose(decomp, padding_size)
        .into_par_iter()
        .map(Vector::from)
        .collect() // decomp_size x v.len()
}

/// Returns the gadget decomposition of a [`Matrix`] of dimensions `m × n` as a matrix of dimensions `m × (k * n)`.
///
/// # Arguments
/// * `mat`: input matrix of dimensions `m × n`
/// * `b`: basis for the decomposition, must be even
/// * `padding_size`: indicates whether the decomposition length is the specified `k` if `padding_size` is `Some(k)`, or if `k` is the largest decomposition length required for `mat` if `padding_size` is `None`
///
/// # Output
/// Returns `d` the decomposition in basis `b` as a Matrix of dimensions `m × (k * n)`, i.e.,
/// $\texttt{mat} = \texttt{d} \times G_\texttt{n}$ where $G_\texttt{n} = I_\texttt{n} \otimes (1, \texttt{b}, \ldots, \texttt{b}^k) \in R^{\texttt{k}\texttt{n} \times \texttt{n}}$ and $|\texttt{d}\[i\]\[j\]| \le \texttt{b}$.
pub fn decompose_matrix<R>(
    mat: &Matrix<R>,
    decomposition_basis: u128,
    decomposition_length: usize,
) -> Matrix<R>
where
    R: Ring + WithSignedRepresentative,
    R::SignedRepresentative: DecompositionFriendlySignedRepresentative,
{
    Matrix::<R>::from_rows(
        mat.par_row_iter()
            .map(|s_i| {
                RowVector::<R>::from(
                    s_i.map(|s_ij| {
                        decompose(&s_ij, decomposition_basis, Some(decomposition_length))
                    })
                    .as_slice()
                    .concat(),
                )
            })
            .collect::<Vec<_>>()
            .as_slice(),
    )
}

#[cfg(test)]
mod tests {
    use crate::decomposition::balanced_decomposition::recompose;
    use crate::decomposition::{recompose_left_right_symmetric_matrix, recompose_matrix};
    use crate::linear_algebra::SymmetricMatrix;
    use ark_std::test_rng;

    use crate::ring;
    use crate::ring::ntt::ntt_prime;
    use crate::ring::pow2_cyclotomic_poly_ring_ntt::Pow2CyclotomicPolyRingNTT;
    use crate::ring::util::powers_of_basis;
    use crate::ring::Zq1;
    use crate::traits::WithLinfNorm;

    use super::*;

    const N: usize = 128;
    const Q: u64 = ntt_prime::<N>(12);
    const VEC_LENGTH: usize = 32;
    const BASIS_TEST_RANGE: [u128; 5] = [2, 4, 8, 16, 32];

    type R = Zq1<Q>;
    type PolyR = Pow2CyclotomicPolyRingNTT<R, N>;

    #[test]
    fn test_decompose() {
        let vs: Vec<R> = (0..Q).map(|v| R::try_from(v).unwrap()).collect();
        for b in BASIS_TEST_RANGE {
            for v in &vs {
                let decomp = decompose(v, b, None);

                // Check that the decomposition length is correct
                let max = v.linf_norm();
                let max_decomp_length = decomposition_max_length(b, max);
                assert!(decomp.len() <= max_decomp_length);

                // Check that all entries are smaller than b/2
                for v_i in &decomp {
                    assert!(v_i.linf_norm() < BigUint::from(b));
                }

                // Check that the decomposition is correct
                assert_eq!(*v, recompose(&decomp, R::try_from(b).unwrap()));
            }
        }
    }

    fn get_test_vec() -> [R; N] {
        core::array::from_fn(|i| R::try_from((i as u64 * Q) / (N as u64)).unwrap())
    }

    #[test]
    fn test_decompose_vec() {
        let v = get_test_vec();
        for b in BASIS_TEST_RANGE {
            let decomp = decompose_vec(&v, b, None);

            // Check that the decomposition length is correct
            let max = v.linf_norm();
            let max_decomp_length = decomposition_max_length(b, max);
            assert!(decomp.len() <= max_decomp_length);

            // Check that all entries are smaller than b in absolute value
            for d_i in &decomp {
                for d_ij in d_i {
                    assert!(d_ij.linf_norm() < b.into());
                }
            }

            for i in 0..decomp.len() {
                // Check that the decomposition is correct
                let decomp_i = decomp.iter().map(|d_j| d_j[i]).collect::<Vec<_>>();
                assert_eq!(v[i], recompose(&decomp_i, R::try_from(b).unwrap()));
            }
        }
    }

    #[test]
    fn test_decompose_polyring() {
        let v = PolyR::from_coefficient_array(get_test_vec());
        for b in BASIS_TEST_RANGE {
            let decomp: Vec<PolyR> = decompose_polyring(&v, b, None);

            for d_i in &decomp {
                for d_ij in d_i.coefficients() {
                    assert!(d_ij.linf_norm() <= b.into());
                }
            }

            assert_eq!(v, recompose(&decomp, R::try_from(b).unwrap()));
        }
    }

    #[test]
    fn test_decompose_vec_polyring() {
        let v = Vector::<PolyR>::from_fn(VEC_LENGTH, |i, _| {
            PolyR::try_from_coefficients(
                get_test_vec()
                    .into_iter()
                    .map(|v| v + R::try_from(i as u64).unwrap())
                    .collect::<Vec<_>>()
                    .as_slice(),
            )
            .unwrap()
        });
        for b in BASIS_TEST_RANGE {
            let decomp = decompose_vec_polyring::<PolyR>(&v.as_slice(), b, None);

            for v_i in &decomp {
                for v_ij in v_i.as_slice() {
                    for v_ijk in v_ij.coefficients() {
                        assert!(v_ijk.linf_norm() < b.into());
                    }
                }
            }

            let mut recomposed = Vector::<PolyR>::zeros(v.len());
            for (i, v_i) in decomp.iter().enumerate() {
                recomposed +=
                    v_i * PolyR::from_scalar(Ring::pow(&R::try_from(b).unwrap(), i as u64));
            }
            assert_eq!(v, recomposed);
        }
    }

    #[test]
    pub fn test_decompose_matrix() {
        use num_traits::Signed;

        const NROWS: usize = 101;
        const NCOLS: usize = 42;

        let rng = &mut test_rng();
        let mat = Matrix::<R>::rand(NROWS, NCOLS, rng);

        for b in BASIS_TEST_RANGE {
            let decomp_size = decomposition_max_length(b, ((Q - 1) as u128).into());

            let decomp = decompose_matrix(&mat, b, decomp_size);
            assert_eq!(decomp.nrows(), mat.nrows());
            assert_eq!(decomp.ncols(), mat.ncols() * decomp_size);
            for d_ij in decomp.iter() {
                assert!(d_ij.as_signed_representative().0.abs() <= b.into());
            }

            let pows_b = ring::util::powers_of_basis(R::try_from(b).unwrap(), decomp_size);
            let recomp = recompose_matrix(&decomp, &pows_b);

            assert_eq!(mat, recomp);

            // Check that recompose is equivalent to right-multiplication by the gadget matrix
            let gadget_vector = Vector::<R>::from_slice(&pows_b);
            let gadget_matrix = Matrix::<R>::identity(NCOLS, NCOLS).kronecker(&gadget_vector);
            let mat_g = &decomp * &gadget_matrix;
            assert_eq!(recomp, mat_g);

            // Check that A * M == recompose(A * decompose(M))
            let mat_l = Matrix::<R>::rand(64, NROWS, rng);
            let lhs = &mat_l * &mat;
            let rhs = recompose_matrix(&(&mat_l * &decomp), &pows_b);
            assert_eq!(lhs, rhs);
        }
    }

    #[test]
    fn test_left_right_decompose_symmetric_matrix() {
        const SIZE: usize = 13;
        let rng = &mut test_rng();

        for b in BASIS_TEST_RANGE {
            let decomp_size = decomposition_max_length(b, (Q - 1).into());

            let pows_b = powers_of_basis(R::try_from(b).unwrap(), decomp_size);
            let mat = SymmetricMatrix::<R>::rand(SIZE * decomp_size, rng);

            let expected: SymmetricMatrix<R> = recompose_matrix(
                &recompose_matrix(&mat.clone().into(), pows_b.as_slice()).transpose(),
                pows_b.as_slice(),
            )
            .into();

            let actual = recompose_left_right_symmetric_matrix(&mat, pows_b.as_slice());
            assert_eq!(expected, actual);
        }
    }
}
