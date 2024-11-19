use ark_ff::{vec::*, UniformRand, Zero};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::rand::Rng;

#[derive(Clone, Debug, Eq, PartialEq, CanonicalSerialize, CanonicalDeserialize)]
pub struct SparseMatrix<
    R1: Clone
        + Send
        + Sync
        + ark_serialize::Valid
        + ark_serialize::CanonicalSerialize
        + ark_serialize::CanonicalDeserialize,
> {
    pub n_rows: usize,
    pub n_cols: usize,
    pub coeffs: Vec<Vec<(R1, usize)>>,
}

impl<
        R: Copy
            + Send
            + Sync
            + Zero
            + UniformRand
            + ark_serialize::Valid
            + ark_serialize::CanonicalSerialize
            + ark_serialize::CanonicalDeserialize,
    > SparseMatrix<R>
{
    pub fn empty() -> Self {
        Self {
            n_rows: 0,
            n_cols: 0,
            coeffs: vec![],
        }
    }

    pub fn rand<RND: Rng>(rng: &mut RND, n_rows: usize, n_cols: usize) -> Self {
        const ZERO_VAL_PROBABILITY: f64 = 0.8f64;

        let dense = (0..n_rows)
            .map(|_| {
                (0..n_cols)
                    .map(|_| {
                        if !rng.gen_bool(ZERO_VAL_PROBABILITY) {
                            return R::rand(rng);
                        }
                        R::zero()
                    })
                    .collect::<Vec<R>>()
            })
            .collect::<Vec<Vec<R>>>();
        dense_matrix_to_sparse(dense)
    }

    pub fn to_dense(&self) -> Vec<Vec<R>> {
        let mut r: Vec<Vec<R>> = vec![vec![R::zero(); self.n_cols]; self.n_rows];
        for (row_i, row) in self.coeffs.iter().enumerate() {
            for &(value, col_i) in row.iter() {
                r[row_i][col_i] = value;
            }
        }
        r
    }

    pub fn nrows(&self) -> usize {
        self.n_rows
    }

    pub fn ncols(&self) -> usize {
        self.n_cols
    }
}

pub fn dense_matrix_to_sparse<
    R: Copy
        + Send
        + Sync
        + Zero
        + UniformRand
        + ark_serialize::Valid
        + ark_serialize::CanonicalSerialize
        + ark_serialize::CanonicalDeserialize,
>(
    m: Vec<Vec<R>>,
) -> SparseMatrix<R> {
    let mut r = SparseMatrix::<R> {
        n_rows: m.len(),
        n_cols: m[0].len(),
        coeffs: Vec::new(),
    };
    for m_row in m.iter() {
        let mut row: Vec<(R, usize)> = Vec::new();
        for (col_i, value) in m_row.iter().enumerate() {
            if !value.is_zero() {
                row.push((*value, col_i));
            }
        }
        r.coeffs.push(row);
    }
    r
}
