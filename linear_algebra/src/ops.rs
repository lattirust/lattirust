use ark_std::{
    ops::{Add, BitXor, Div, Sub},
    Zero,
};
use num_integer::Integer;

pub trait Transpose {
    fn transpose(&self) -> Self;
}

impl<F: Zero + Clone> Transpose for Vec<Vec<F>> {
    fn transpose(&self) -> Self {
        let nrows = self.len();
        let ncols = self.iter().map(|d_i| d_i.len()).max().unwrap_or(0);

        let mut res: Vec<Vec<_>> = (0..ncols).map(|_| Vec::with_capacity(nrows)).collect();

        for row in self {
            // Copy existing values from original rows to new cols
            for (c, value) in row.iter().enumerate() {
                res[c].push(value.clone());
            }

            // Pad the shorter rows with zeroes
            for res_row in res.iter_mut().take(ncols).skip(row.len()) {
                res_row.push(F::zero());
            }
        }

        res
    }
}

pub fn rounded_div<T, D>(dividend: T, divisor: D) -> T
where
    T: Integer
        + BitXor<Output = T>
        + Add<D, Output = T>
        + Div<D, Output = T>
        + Sub<D, Output = T>
        + From<D>
        + Clone,
    D: Clone + Div<i128, Output = D>,
{
    if dividend.clone() ^ divisor.clone().into() >= T::zero() {
        (dividend + (divisor.clone() / 2)) / divisor
    } else {
        (dividend - (divisor.clone() / 2)) / divisor
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[rustfmt::skip]
    fn transpose_vec_of_vecs() {
        let v = vec![
            vec![1, 2, 3],
            vec![9],
            vec![7, 7, 7]
        ].transpose();

        #[rustfmt::skip]
        let r = vec![
            vec![1, 9, 7],
            vec![2, 0, 7],
            vec![3, 0, 7]
        ];

        assert_eq!(v, r);
    }
}
