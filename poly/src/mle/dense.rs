use ark_std::ops::{Add, AddAssign, Index, Neg, Sub, SubAssign};

use ark_ff::Zero;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::cfg_iter;

use super::{swap_bits, MultilinearExtension};
use lattirust_ring::Ring;
use rayon::iter::{IndexedParallelIterator, IntoParallelRefIterator, ParallelIterator};

#[derive(Debug, Clone, PartialEq, Eq, Hash, Default, CanonicalDeserialize, CanonicalSerialize)]
pub struct DenseMultilinearExtension<Rn: Ring> {
    /// The evaluation over {0,1}^`num_vars`
    pub evaluations: Vec<Rn>,
    /// Number of variables
    pub num_vars: usize,
}

impl<R: Ring> DenseMultilinearExtension<R> {
    pub fn from_evaluations_slice(num_vars: usize, evaluations: &[R]) -> Self {
        Self::from_evaluations_vec(num_vars, evaluations.to_vec())
    }

    pub fn evaluate(&self, point: &[R]) -> Option<R> {
        if point.len() == self.num_vars {
            Some(self.fix_variables(point)[0])
        } else {
            None
        }
    }

    pub fn from_evaluations_vec(num_vars: usize, evaluations: Vec<R>) -> Self {
        // assert that the number of variables matches the size of evaluations
        assert_eq!(
            evaluations.len(),
            1 << num_vars,
            "The size of evaluations should be 2^num_vars."
        );

        Self {
            num_vars,
            evaluations,
        }
    }

    pub fn relabel_in_place(&mut self, mut a: usize, mut b: usize, k: usize) {
        // enforce order of a and b
        if a > b {
            std::mem::swap(&mut a, &mut b);
        }
        if a == b || k == 0 {
            return;
        }
        assert!(b + k <= self.num_vars, "invalid relabel argument");
        assert!(a + k <= b, "overlapped swap window is not allowed");
        for i in 0..self.evaluations.len() {
            let j = swap_bits(i, a, b, k);
            if i < j {
                self.evaluations.swap(i, j);
            }
        }
    }
}

impl<R: Ring> MultilinearExtension<R> for DenseMultilinearExtension<R> {
    fn num_vars(&self) -> usize {
        self.num_vars
    }

    fn rand<Rn: rand::Rng>(num_vars: usize, rng: &mut Rn) -> Self {
        Self::from_evaluations_vec(num_vars, (0..1 << num_vars).map(|_| R::rand(rng)).collect())
    }

    fn relabel(&self, a: usize, b: usize, k: usize) -> Self {
        let mut copy = self.clone();
        copy.relabel_in_place(a, b, k);
        copy
    }

    fn fix_variables(&self, partial_point: &[R]) -> Self {
        assert!(
            partial_point.len() <= self.num_vars,
            "too many partial points"
        );

        let mut poly = self.evaluations.to_vec();
        let nv = self.num_vars;
        let dim = partial_point.len();

        for i in 1..dim + 1 {
            let r = partial_point[i - 1];
            for b in 0..1 << (nv - i) {
                let left = poly[b << 1];
                let right = poly[(b << 1) + 1];
                poly[b] = left + r * (right - left);
            }
        }
        Self::from_evaluations_slice(nv - dim, &poly[..1 << (nv - dim)])
    }

    fn to_evaluations(&self) -> Vec<R> {
        self.evaluations.to_vec()
    }
}
impl<R: Ring> Zero for DenseMultilinearExtension<R> {
    fn zero() -> Self {
        Self {
            num_vars: 0,
            evaluations: vec![R::zero()],
        }
    }

    fn is_zero(&self) -> bool {
        self.num_vars == 0 && self.evaluations[0].is_zero()
    }
}
impl<R: Ring> Add for DenseMultilinearExtension<R> {
    type Output = DenseMultilinearExtension<R>;

    fn add(self, other: DenseMultilinearExtension<R>) -> Self {
        &self + &other
    }
}
impl<'a, 'b, R: Ring> Add<&'a DenseMultilinearExtension<R>> for &'b DenseMultilinearExtension<R> {
    type Output = DenseMultilinearExtension<R>;

    fn add(self, rhs: &'a DenseMultilinearExtension<R>) -> Self::Output {
        if rhs.is_zero() {
            return self.clone();
        }
        if self.is_zero() {
            return rhs.clone();
        }
        assert_eq!(
            self.num_vars, rhs.num_vars,
            "trying to add two dense MLEs with different numbers of variables"
        );
        let result: Vec<R> = cfg_iter!(self.evaluations)
            .zip(cfg_iter!(rhs.evaluations))
            .map(|(a, b)| *a + *b)
            .collect();

        Self::Output::from_evaluations_vec(self.num_vars, result)
    }
}
impl<R: Ring> AddAssign for DenseMultilinearExtension<R> {
    fn add_assign(&mut self, rhs: DenseMultilinearExtension<R>) {
        *self = &*self + &rhs;
    }
}
impl<'a, R: Ring> AddAssign<&'a DenseMultilinearExtension<R>> for DenseMultilinearExtension<R> {
    fn add_assign(&mut self, rhs: &'a DenseMultilinearExtension<R>) {
        *self = &*self + rhs;
    }
}
impl<R: Ring> AddAssign<(R, &DenseMultilinearExtension<R>)> for DenseMultilinearExtension<R>
where
    R: Copy + std::ops::AddAssign,
{
    fn add_assign(&mut self, (r, other): (R, &DenseMultilinearExtension<R>)) {
        let other = Self {
            num_vars: other.num_vars,
            evaluations: cfg_iter!(other.evaluations).map(|x| r + x).collect(),
        };
        *self = &*self + &other;
    }
}
impl<R: Ring> Neg for DenseMultilinearExtension<R> {
    type Output = DenseMultilinearExtension<R>;

    fn neg(self) -> Self {
        Self::Output {
            num_vars: self.num_vars,
            evaluations: cfg_iter!(self.evaluations).map(|x| -*x).collect(),
        }
    }
}
impl<R: Ring> Sub for DenseMultilinearExtension<R> {
    type Output = DenseMultilinearExtension<R>;

    fn sub(self, other: DenseMultilinearExtension<R>) -> Self {
        &self - &other
    }
}
impl<'a, 'b, R: Ring> Sub<&'a DenseMultilinearExtension<R>> for &'b DenseMultilinearExtension<R> {
    type Output = DenseMultilinearExtension<R>;

    #[allow(clippy::suspicious_arithmetic_impl)]
    fn sub(self, rhs: &'a DenseMultilinearExtension<R>) -> Self::Output {
        self + &rhs.clone().neg()
    }
}
impl<R: Ring> SubAssign for DenseMultilinearExtension<R> {
    fn sub_assign(&mut self, other: DenseMultilinearExtension<R>) {
        *self = &*self - &other;
    }
}
impl<'a, R: Ring> SubAssign<&'a DenseMultilinearExtension<R>> for DenseMultilinearExtension<R> {
    fn sub_assign(&mut self, rhs: &'a DenseMultilinearExtension<R>) {
        *self = &*self - rhs;
    }
}
impl<R: Ring> Index<usize> for DenseMultilinearExtension<R> {
    type Output = R;

    fn index(&self, index: usize) -> &Self::Output {
        &self.evaluations[index]
    }
}
