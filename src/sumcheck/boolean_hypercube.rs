use ark_ff::Field;

pub struct BooleanHypercube<F: Field> {
    dim: u8,
    current: u128,
    _marker: std::marker::PhantomData<F>,
}

impl<F: Field> BooleanHypercube<F> {
    pub fn new(dim: u8) -> BooleanHypercube<F> {
        assert!(dim < 128, "only boolean hypercubes {{0,1}}^n with n < 128 are supported, got {dim}");
        BooleanHypercube {
            dim,
            current: 0,
            _marker: std::marker::PhantomData,
        }
    }
}

impl<F: Field> Iterator for BooleanHypercube<F> {
    type Item = Vec<F>;

    fn next(&mut self) -> Option<Self::Item> {
        let max = 2u128.pow(self.dim as u32);
        if self.current == max {
            None
        } else {
            let res = (0..self.dim).map(|i| F::from((self.current >> i) & 1)).collect();
            self.current += 1;
            Some(res)
        }
    }
}

#[cfg(test)]
mod tests {
    use num_traits::{One, Zero};
    use crate::lattice_arithmetic::ring::Fq;
    use super::BooleanHypercube;

    type F = Fq<655357>;

    #[test]
    pub fn test_boolean_hypercube() {
        let n: u8 = 16;
        let pows = (0..n).map(|i| F::from(2u128.pow(i as u32)));
        assert!(BooleanHypercube::<F>::new(n).enumerate().all(
            |(i, b)| b.iter().all(|b_i| b_i.is_zero() || b_i.is_one())
                && b.iter().zip(pows.to_owned()).map(|(b_i, t_i)| b_i * &t_i).sum::<F>() == F::from(i as u128)
        ));
    }
}