use lattirust_arithmetic::ring::Ring;

pub struct BooleanHypercube<F: Ring> {
    dim: u8,
    current: u128,
    _marker: std::marker::PhantomData<F>,
}

impl<F: Ring> BooleanHypercube<F> {
    pub fn new(dim: u8) -> BooleanHypercube<F> {
        assert!(
            dim < 128,
            "only boolean hypercubes {{0,1}}^n with n < 128 are supported, got {dim}"
        );
        BooleanHypercube {
            dim,
            current: 0,
            _marker: std::marker::PhantomData,
        }
    }
}

impl<F: Ring> Iterator for BooleanHypercube<F> {
    type Item = Vec<F>;

    fn next(&mut self) -> Option<Self::Item> {
        let max = 2u128.pow(self.dim as u32);
        if self.current == max {
            None
        } else {
            let res = (0..self.dim)
                .map(|i| F::try_from((self.current >> i) & 1).unwrap())
                .collect();
            self.current += 1;
            Some(res)
        }
    }
}

#[cfg(test)]
mod tests {
    use lattirust_arithmetic::ring::Zq1;
    use num_traits::identities::{One, Zero};
    use super::BooleanHypercube;

    type F = Zq1<655357>;

    #[test]
    pub fn test_boolean_hypercube() {
        let n: u8 = 16;
        let pows = (0..n).map(|i| F::try_from(2u128.pow(i as u32)).unwrap());
        assert!(BooleanHypercube::<F>::new(n).enumerate().all(|(i, b)| b
            .iter()
            .all(|b_i| b_i.is_zero() || b_i.is_one())
            && b.iter()
                .zip(pows.to_owned())
                .map(|(b_i, t_i)| b_i * &t_i)
                .sum::<F>()
                == F::try_from(i as u128).unwrap()));
    }
}
