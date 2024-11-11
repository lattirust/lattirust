use ark_ff::{UniformRand, Zero};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::{
    cfg_iter,
    collections::{BTreeMap, HashMap},
    hash::{DefaultHasher, Hash},
    log2,
    ops::{Add, AddAssign, Index, Neg, Sub, SubAssign},
};
use lattirust_linear_algebra::SparseMatrix;
use rand::Rng;
#[cfg(feature = "parallel")]
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};

use super::{swap_bits, MultilinearExtension};
use lattirust_ring::Ring;

#[derive(Debug, Clone, PartialEq, Eq, Hash, Default, CanonicalDeserialize, CanonicalSerialize)]
pub struct SparseMultilinearExtension<Rn: Ring> {
    /// The evaluation over {0,1}^`num_vars`
    pub evaluations: BTreeMap<usize, Rn>,
    /// Number of variables
    pub num_vars: usize,
    zero: Rn,
}
impl<R: Ring> SparseMultilinearExtension<R> {
    pub fn from_evaluations<'a>(
        num_vars: usize,
        evaluations: impl IntoIterator<Item = &'a (usize, R)>,
    ) -> Self {
        let bit_mask = 1 << num_vars;

        let evaluations: Vec<_> = evaluations
            .into_iter()
            .map(|(i, v): &(usize, R)| {
                assert!(*i < bit_mask, "index out of range");
                (*i, *v)
            })
            .collect();
        Self {
            evaluations: tuples_to_treemap(&evaluations),
            num_vars,
            zero: R::zero(),
        }
    }
    pub fn evaluate(&self, point: &[R]) -> R {
        assert!(point.len() == self.num_vars);
        self.fix_variables(point)[0]
    }
    /// Outputs an `l`-variate multilinear extension where value of evaluations
    /// are sampled uniformly at random. The number of nonzero entries is
    /// `num_nonzero_entries` and indices of those nonzero entries are
    /// distributed uniformly at random.
    ///
    /// Note that this function uses rejection sampling. As number of nonzero
    /// entries approach `2 ^ num_vars`, sampling will be very slow due to
    /// large number of collisions.
    pub fn rand_with_config<Rn: Rng>(
        num_vars: usize,
        num_nonzero_entries: usize,
        rng: &mut Rn,
    ) -> Self {
        assert!(num_nonzero_entries <= 1 << num_vars);

        let mut map =
            HashMap::with_hasher(core::hash::BuildHasherDefault::<DefaultHasher>::default());
        for _ in 0..num_nonzero_entries {
            let mut index = usize::rand(rng) & ((1 << num_vars) - 1);
            while map.contains_key(&index) {
                index = usize::rand(rng) & ((1 << num_vars) - 1);
            }
            map.entry(index).or_insert(R::rand(rng));
        }
        let mut buf = Vec::new();
        for (arg, v) in map.iter() {
            if *v != R::zero() {
                buf.push((*arg, *v));
            }
        }
        let evaluations = hashmap_to_treemap(&map);
        Self {
            num_vars,
            evaluations,
            zero: R::zero(),
        }
    }

    /// Returns the sparse MLE from the given matrix, without modifying the original matrix.
    pub fn from_matrix(m: &SparseMatrix<R>) -> Self {
        let n_rows = m.nrows().next_power_of_two();
        let n_cols = m.ncols().next_power_of_two();
        let n_vars: usize = (log2(n_rows * n_cols)) as usize; // n_vars = s + s'

        // build the sparse vec representing the sparse matrix
        let mut v: Vec<(usize, R)> = Vec::with_capacity(m.nnz());

        for (col, row, val) in m.triplet_iter().filter(|(_, _, val)| !val.is_zero()) {
            v.push((row * n_cols + col, *val));
        }

        // convert the sparse vector into a mle
        Self::from_sparse_slice(n_vars, &v)
    }

    /// Takes n_vars and a sparse slice and returns its sparse MLE.
    pub fn from_sparse_slice(n_vars: usize, v: &[(usize, R)]) -> Self {
        SparseMultilinearExtension::<R>::from_evaluations(n_vars, v)
    }

    /// Takes n_vars and a dense slice and returns its sparse MLE.
    pub fn from_slice(n_vars: usize, v: &[R]) -> Self {
        let v_sparse = v
            .iter()
            .enumerate()
            .map(|(i, v_i)| (i, *v_i))
            .collect::<Vec<(usize, R)>>();
        SparseMultilinearExtension::<R>::from_evaluations(n_vars, &v_sparse)
    }
}

impl<R: Ring> MultilinearExtension<R> for SparseMultilinearExtension<R> {
    fn num_vars(&self) -> usize {
        self.num_vars
    }
    /// Outputs an `l`-variate multilinear extension where value of evaluations
    /// are sampled uniformly at random. The number of nonzero entries is
    /// `sqrt(2^num_vars)` and indices of those nonzero entries are distributed
    /// uniformly at random.
    fn rand<Rn: rand::Rng>(num_vars: usize, rng: &mut Rn) -> Self {
        Self::rand_with_config(num_vars, 1 << (num_vars / 2), rng)
    }

    fn relabel(&self, mut a: usize, mut b: usize, k: usize) -> Self {
        if a > b {
            // swap
            core::mem::swap(&mut a, &mut b);
        }
        // sanity check
        assert!(
            a + k < self.num_vars && b + k < self.num_vars,
            "invalid relabel argument"
        );
        if a == b || k == 0 {
            return self.clone();
        }
        assert!(a + k <= b, "overlapped swap window is not allowed");
        let ev: Vec<_> = cfg_iter!(self.evaluations)
            .map(|(&i, &v)| (swap_bits(i, a, b, k), v))
            .collect();
        Self {
            num_vars: self.num_vars,
            evaluations: tuples_to_treemap(&ev),
            zero: R::zero(),
        }
    }

    fn fix_variables(&self, partial_point: &[R]) -> Self {
        let dim = partial_point.len();
        assert!(dim <= self.num_vars, "invalid partial point dimension");

        let mut window = ark_std::log2(self.evaluations.len()) as usize;
        if window == 0 {
            window = 1;
        }
        let mut point = partial_point;
        let mut last = treemap_to_hashmap(&self.evaluations);

        // batch evaluation
        while !point.is_empty() {
            let focus_length = if point.len() > window {
                window
            } else {
                point.len()
            };
            let focus = &point[..focus_length];
            point = &point[focus_length..];
            let pre = precompute_eq(focus);
            let dim = focus.len();
            let mut result =
                HashMap::with_hasher(core::hash::BuildHasherDefault::<DefaultHasher>::default());
            for src_entry in last.iter() {
                let old_idx = *src_entry.0;
                let gz = pre[old_idx & ((1 << dim) - 1)];
                let new_idx = old_idx >> dim;
                let dst_entry = result.entry(new_idx).or_insert(R::zero());
                *dst_entry += gz * src_entry.1;
            }
            last = result;
        }
        let evaluations = hashmap_to_treemap(&last);
        Self {
            num_vars: self.num_vars - dim,
            evaluations,
            zero: R::zero(),
        }
    }
    fn to_evaluations(&self) -> Vec<R> {
        let mut evaluations: Vec<_> = (0..1 << self.num_vars).map(|_| R::zero()).collect();
        self.evaluations
            .iter()
            .map(|(&i, &v)| {
                evaluations[i] = v;
            })
            .last();
        evaluations
    }
}
impl<R: Ring> Zero for SparseMultilinearExtension<R> {
    fn zero() -> Self {
        Self {
            num_vars: 0,
            evaluations: tuples_to_treemap(&Vec::new()),
            zero: R::zero(),
        }
    }

    fn is_zero(&self) -> bool {
        self.num_vars == 0 && self.evaluations.is_empty()
    }
}
impl<R: Ring> Add for SparseMultilinearExtension<R> {
    type Output = SparseMultilinearExtension<R>;

    fn add(self, other: SparseMultilinearExtension<R>) -> Self {
        &self + &other
    }
}
impl<'a, 'b, R: Ring> Add<&'a SparseMultilinearExtension<R>> for &'b SparseMultilinearExtension<R> {
    type Output = SparseMultilinearExtension<R>;

    fn add(self, rhs: &'a SparseMultilinearExtension<R>) -> Self::Output {
        // handle zero case
        if self.is_zero() {
            return rhs.clone();
        }
        if rhs.is_zero() {
            return self.clone();
        }

        assert_eq!(
            rhs.num_vars, self.num_vars,
            "trying to add non-zero polynomial with different number of variables"
        );
        // simply merge the evaluations
        let mut evaluations =
            HashMap::with_hasher(core::hash::BuildHasherDefault::<DefaultHasher>::default());
        for (&i, &v) in self.evaluations.iter().chain(rhs.evaluations.iter()) {
            *evaluations.entry(i).or_insert(R::zero()) += v;
        }
        let evaluations: Vec<_> = evaluations
            .into_iter()
            .filter(|(_, v)| !v.is_zero())
            .collect();

        Self::Output {
            evaluations: tuples_to_treemap(&evaluations),
            num_vars: self.num_vars,
            zero: R::zero(),
        }
    }
}

impl<R: Ring> AddAssign for SparseMultilinearExtension<R> {
    fn add_assign(&mut self, other: Self) {
        *self = &*self + &other;
    }
}
impl<'a, R: Ring> AddAssign<&'a SparseMultilinearExtension<R>> for SparseMultilinearExtension<R> {
    fn add_assign(&mut self, rhs: &'a SparseMultilinearExtension<R>) {
        *self = &*self + rhs;
    }
}
impl<'a, R: Ring> AddAssign<(R, &'a SparseMultilinearExtension<R>)>
    for SparseMultilinearExtension<R>
{
    fn add_assign(&mut self, (r, other): (R, &SparseMultilinearExtension<R>)) {
        if !self.is_zero() && !other.is_zero() {
            assert_eq!(
                other.num_vars, self.num_vars,
                "trying to add non-zero polynomial with different number of variables"
            );
        }
        let ev: Vec<_> = cfg_iter!(other.evaluations)
            .map(|(i, v)| (*i, r * v))
            .collect();
        let other = Self {
            num_vars: other.num_vars,
            evaluations: tuples_to_treemap(&ev),
            zero: R::zero(),
        };
        *self += &other;
    }
}
impl<R: Ring> Neg for SparseMultilinearExtension<R> {
    type Output = SparseMultilinearExtension<R>;

    fn neg(self) -> Self {
        let ev: Vec<_> = cfg_iter!(self.evaluations)
            .map(|(i, v)| (*i, -*v))
            .collect();
        Self::Output {
            num_vars: self.num_vars,
            evaluations: tuples_to_treemap(&ev),
            zero: R::zero(),
        }
    }
}
impl<R: Ring> Sub for SparseMultilinearExtension<R> {
    type Output = SparseMultilinearExtension<R>;

    fn sub(self, other: SparseMultilinearExtension<R>) -> Self {
        &self - &other
    }
}
impl<'a, 'b, R: Ring> Sub<&'a SparseMultilinearExtension<R>> for &'b SparseMultilinearExtension<R> {
    type Output = SparseMultilinearExtension<R>;

    #[allow(clippy::suspicious_arithmetic_impl)]
    fn sub(self, rhs: &'a SparseMultilinearExtension<R>) -> Self::Output {
        self + &rhs.clone().neg()
    }
}
impl<R: Ring> SubAssign for SparseMultilinearExtension<R> {
    fn sub_assign(&mut self, other: SparseMultilinearExtension<R>) {
        *self = &*self - &other;
    }
}
impl<'a, R: Ring> SubAssign<&'a SparseMultilinearExtension<R>> for SparseMultilinearExtension<R> {
    fn sub_assign(&mut self, rhs: &'a SparseMultilinearExtension<R>) {
        *self = &*self - rhs;
    }
}
impl<R: Ring> Index<usize> for SparseMultilinearExtension<R> {
    type Output = R;

    /// Returns the evaluation of the polynomial at a point represented by
    /// index.
    ///
    /// Index represents a vector in {0,1}^`num_vars` in little endian form. For
    /// example, `0b1011` represents `P(1,1,0,1)`
    ///
    /// For Sparse multilinear polynomial, Lookup_evaluation takes log time to
    /// the size of polynomial.
    fn index(&self, index: usize) -> &Self::Output {
        if let Some(v) = self.evaluations.get(&index) {
            v
        } else {
            &self.zero
        }
    }
}

/// Utilities
fn tuples_to_treemap<R: Ring>(tuples: &[(usize, R)]) -> BTreeMap<usize, R> {
    BTreeMap::from_iter(tuples.iter().map(|(i, v)| (*i, *v)))
}

fn treemap_to_hashmap<R: Ring>(
    map: &BTreeMap<usize, R>,
) -> HashMap<usize, R, core::hash::BuildHasherDefault<DefaultHasher>> {
    HashMap::from_iter(map.iter().map(|(i, v)| (*i, *v)))
}
fn hashmap_to_treemap<R: Ring, S>(map: &HashMap<usize, R, S>) -> BTreeMap<usize, R> {
    BTreeMap::from_iter(map.iter().map(|(i, v)| (*i, *v)))
}

// precompute  f(x) = eq(g,x)
fn precompute_eq<R: Ring>(g: &[R]) -> Vec<R> {
    let dim = g.len();
    let mut dp = vec![R::zero(); 1 << dim];
    dp[0] = R::one() - g[0];
    dp[1] = g[0];
    for i in 1..dim {
        for b in 0..1 << i {
            let prev = dp[b];
            dp[b + (1 << i)] = prev * g[i];
            dp[b] = prev - dp[b + (1 << i)];
        }
    }
    dp
}

#[cfg(test)]
#[allow(non_snake_case)]
mod tests {
    use super::*;

    use ark_ff::Zero;
    use lattirust_ring::cyclotomic_ring::models::goldilocks::Fq;

    // Function to convert usize to a binary vector of Ring elements.
    fn usize_to_binary_vector<R: Ring>(n: usize, dimensions: usize) -> Vec<R> {
        let mut bits = Vec::with_capacity(dimensions);
        let mut current = n;

        for _ in 0..dimensions {
            if (current & 1) == 1 {
                bits.push(R::one());
            } else {
                bits.push(R::zero());
            }
            current >>= 1;
        }
        bits
    }

    // Wrapper function to generate a boolean hypercube.
    fn boolean_hypercube<R: Ring>(dimensions: usize) -> Vec<Vec<R>> {
        let max_val = 1 << dimensions; // 2^dimensions
        (0..max_val)
            .map(|i| usize_to_binary_vector::<R>(i, dimensions))
            .collect()
    }

    fn vec_cast<R: Ring>(v: &[usize]) -> Vec<R> {
        v.iter().map(|c| R::from(*c as u64)).collect()
    }

    fn matrix_cast<R: Ring>(m: &[Vec<usize>]) -> SparseMatrix<R> {
        m.iter()
            .map(|r| vec_cast::<R>(r))
            .collect::<Vec<Vec<R>>>()
            .as_slice()
            .into()
    }

    fn get_test_z<R: Ring>(input: usize) -> Vec<R> {
        vec_cast(&[
            input, // io
            1,
            input * input * input + input + 5, // x^3 + x + 5
            input * input,                     // x^2
            input * input * input,             // x^2 * x
            input * input * input + input,     // x^3 + x
        ])
    }

    #[test]
    fn test_matrix_to_mle() {
        type R = Fq;
        let A = matrix_cast::<R>(&[
            vec![2, 3, 4, 4],
            vec![4, 11, 14, 14],
            vec![2, 8, 17, 17],
            vec![420, 4, 2, 0],
        ]);

        let A_mle = SparseMultilinearExtension::from_matrix(&A);
        assert_eq!(A_mle.evaluations.len(), 15); // 15 non-zero elements
        assert_eq!(A_mle.num_vars, 4); // 4x4 matrix, thus 2bit x 2bit, thus 2^4=16 evals

        let A = matrix_cast::<R>(&[
            vec![2, 3, 4, 4, 1],
            vec![4, 11, 14, 14, 2],
            vec![2, 8, 17, 17, 3],
            vec![420, 4, 2, 0, 4],
            vec![420, 4, 2, 0, 5],
        ]);
        let A_mle = SparseMultilinearExtension::from_matrix(&A);
        assert_eq!(A_mle.evaluations.len(), 23); // 23 non-zero elements
        assert_eq!(A_mle.num_vars, 6); // 5x5 matrix, thus 3bit x 3bit, thus 2^6=64 evals
    }

    #[test]
    fn test_vec_to_mle() {
        type R = Fq;
        let z = get_test_z::<R>(3);
        let n_vars = 3;
        let z_mle = SparseMultilinearExtension::from_slice(n_vars, &z);

        // check that the z_mle evaluated over the boolean hypercube equals the vec z_i values
        let bhc = boolean_hypercube(z_mle.num_vars);

        for (i, z_i) in z.iter().enumerate() {
            let s_i = &bhc[i];
            assert_eq!(z_mle.evaluate(s_i), z_i.clone());
        }
        // for the rest of elements of the boolean hypercube, expect it to evaluate to zero
        for s_i in bhc.iter().take(1 << z_mle.num_vars).skip(z.len()) {
            assert_eq!(z_mle.fix_variables(s_i)[0], R::zero());
        }
    }
}
