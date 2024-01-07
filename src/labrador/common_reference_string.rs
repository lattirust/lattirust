use crate::lattice_arithmetic::matrix::{Matrix, sample_uniform_mat, sample_uniform_vec, Vector};
use crate::lattice_arithmetic::poly_ring::PolyRing;
use crate::lattice_arithmetic::traits::WithLog2;

#[derive(Clone, Debug)]
pub struct CommonReferenceString<R: PolyRing> {
    pub r: usize,
    pub n: usize,
    pub d: usize,
    pub beta: f64,
    pub k: usize,
    pub k1: usize,
    pub k2: usize,
    pub t1: usize,
    pub t2: usize,
    pub num_aggregs: usize,
    pub A: Matrix<R>,
    // k x n
    pub B: Vec<Vec<Matrix<R>>>,
    // (r x t1) x k1 x k
    pub C: Vec<Vec<Vec<Vector<R>>>>,
    // (r x r x t2) x k2 x 1
    pub D: Vec<Vec<Vec<Vector<R>>>>,
    // (r x r x t1) x k2 x 1
    pub decomposition_basis: R::BaseRing,
}

impl<R: PolyRing> CommonReferenceString<R> {
    pub fn new(r: usize, n: usize, d: usize, beta: f64, k: usize, k1: usize, k2: usize, decomposition_basis: R::BaseRing) -> CommonReferenceString<R> {
        let log2_q: f64 = R::BaseRing::log2_q();
        let log2_b: f64 = R::BaseRing::log2(&decomposition_basis);
        let t1 = (log2_q / log2_b).round() as usize;
        let s_std_dev_sq: f64 = (beta * beta) / ((r * n * d) as f64); // standard deviation of s vectors = beta / sqrt(r * n * d)
        let t2 = (f64::log2(f64::sqrt((24 * n * d) as f64) * s_std_dev_sq) / log2_b).round() as usize;
        let num_aggregs = (128. / log2_q).ceil() as usize;

        CommonReferenceString {
            r,
            n,
            d,
            beta,
            k,
            k1,
            k2,
            t1,
            t2,
            num_aggregs,
            A: sample_uniform_mat(k, n),
            B: (0..r).map(
                |_| (0..t1).map(
                    |_| sample_uniform_mat(k1, k)
                ).collect()
            ).collect(),
            C: (0..r).map(
                |i| (i..r).map(
                    |_|
                        (0..t2).map(
                            |_| sample_uniform_vec(k1)
                        ).collect()
                ).collect()
            ).collect(),
            D: (0..r).map(
                |i| (i..r).map(
                    |_|
                        (0..t1).map(
                            |_| sample_uniform_vec(k2)
                        ).collect()
                ).collect()
            ).collect(),
            decomposition_basis,
        }
    }
}