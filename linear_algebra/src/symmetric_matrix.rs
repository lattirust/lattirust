#![allow(non_snake_case)]

use ark_std::io::{Read, Write};
use ark_std::ops::{Index, IndexMut};

use ark_ff::Zero;
use ark_serialize::{
    CanonicalDeserialize, CanonicalSerialize, Compress, SerializationError, Valid, Validate,
};
use ark_std::{rand, UniformRand};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use crate::Matrix;
use crate::Scalar;

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, Hash)]
pub struct SymmetricMatrix<F: Clone>(Vec<Vec<F>>);

impl<F: Clone> From<Vec<Vec<F>>> for SymmetricMatrix<F> {
    fn from(value: Vec<Vec<F>>) -> Self {
        assert!(value.iter().enumerate().all(|(i, v_i)| v_i.len() == i + 1), "cannot convert value: Vec<Vec<F>> to SymmetricMatrix<F>, row has wrong number of entries");
        Self(value)
    }
}

impl<F: Clone + Scalar> From<Matrix<F>> for SymmetricMatrix<F> {
    fn from(value: Matrix<F>) -> Self {
        assert_eq!(value.transpose(), value);
        Self(
            value
                .row_iter()
                .enumerate()
                .map(|(i, v_i)| {
                    v_i.iter()
                        .take(i + 1)
                        .into_iter()
                        .map(|v| v.clone())
                        .collect()
                })
                .collect(),
        )
    }
}

impl<F: Clone + Scalar> Into<Matrix<F>> for SymmetricMatrix<F> {
    fn into(self) -> Matrix<F> {
        Matrix::<F>::from_fn(self.size(), self.size(), |i, j| self.at(i, j).clone())
    }
}

impl<F: Zero + Clone> SymmetricMatrix<F> {
    pub fn zero(n: usize) -> SymmetricMatrix<F> {
        SymmetricMatrix::<F>((0..n).map(|i| vec![F::zero(); i + 1]).collect())
    }
}

impl<F: Scalar> PartialEq<Matrix<F>> for SymmetricMatrix<F> {
    fn eq(&self, other: &Matrix<F>) -> bool {
        self.0.iter().enumerate().all(|(i, self_i)| {
            self_i
                .into_iter()
                .enumerate()
                .all(|(j, self_ij)| other[(i, j)] == *self_ij)
        })
    }
}

impl<F: Clone> SymmetricMatrix<F> {
    #[inline]
    pub fn size(&self) -> usize {
        self.0.len()
    }

    #[inline]
    pub fn at(&self, i: usize, j: usize) -> &F {
        debug_assert!(i < self.0.len() && j < self.0.len());
        if j <= i {
            &self.0[i][j]
        } else {
            &self.0[j][i]
        }
    }
    #[inline]
    pub fn at_mut(&mut self, i: usize, j: usize) -> &mut F {
        debug_assert!(i < self.0.len() && j < self.0.len());
        if j <= i {
            &mut self.0[i][j]
        } else {
            &mut self.0[j][i]
        }
    }

    pub fn diag(&self) -> Vec<F> {
        (0..self.size()).map(|i| self.at(i, i).clone()).collect()
    }

    pub fn rows(&self) -> &Vec<Vec<F>> {
        &self.0
    }

    pub fn map<T, M>(&self, func: M) -> SymmetricMatrix<T>
    where
        T: Clone,
        M: Fn(&F) -> T,
    {
        SymmetricMatrix::<T>::from(
            self.rows()
                .into_iter()
                .map(|row| row.into_iter().map(&func).collect())
                .collect::<Vec<Vec<T>>>(),
        )
    }

    pub fn from_par_fn<Func>(size: usize, func: Func) -> Self
    where
        F: Send + Sync,
        Func: Send + Sync + Fn(usize, usize) -> F,
    {
        Self::from(
            (0..size)
                .into_par_iter()
                .map(|i| (0..i + 1).into_par_iter().map(|j| func(i, j)).collect())
                .collect::<Vec<Vec<F>>>(),
        )
    }
}
impl<F: Clone + Scalar> SymmetricMatrix<F> {
    pub fn from_blocks(
        top_left: SymmetricMatrix<F>,
        bottom_left: Matrix<F>,
        bottom_right: SymmetricMatrix<F>,
    ) -> Self {
        let n = top_left.size();
        assert_eq!(bottom_left.nrows(), n);
        assert_eq!(bottom_left.ncols(), n);
        assert_eq!(bottom_right.size(), n);

        let mut result = top_left.0;
        result.extend(
            bottom_left
                .row_iter()
                .zip(bottom_right.0.into_iter())
                .map(|(bl_i, br_i)| [&bl_i.0.into_owned().as_slice(), br_i.as_slice()].concat()),
        );
        Self(result)
    }
}

impl<F: Clone> Index<(usize, usize)> for SymmetricMatrix<F> {
    type Output = F;

    fn index(&self, index: (usize, usize)) -> &Self::Output {
        self.at(index.0, index.1)
    }
}

impl<F: Clone> IndexMut<(usize, usize)> for SymmetricMatrix<F> {
    fn index_mut(&mut self, index: (usize, usize)) -> &mut Self::Output {
        self.at_mut(index.0, index.1)
    }
}

impl<F: Clone + UniformRand> SymmetricMatrix<F> {
    pub fn rand<Rng: rand::Rng + ?Sized>(n: usize, rng: &mut Rng) -> SymmetricMatrix<F> {
        SymmetricMatrix::<F>(
            (0..n)
                .map(|i| (0..i + 1).map(|_| F::rand(rng)).collect())
                .collect(),
        )
    }
}

impl<F: Clone> CanonicalSerialize for SymmetricMatrix<F>
where
    Vec<Vec<F>>: CanonicalSerialize,
{
    fn serialize_with_mode<W: Write>(
        &self,
        writer: W,
        compress: Compress,
    ) -> Result<(), SerializationError> {
        self.0.serialize_with_mode(writer, compress)
    }

    fn serialized_size(&self, compress: Compress) -> usize {
        self.0.serialized_size(compress)
    }
}

impl<F: Clone> Valid for SymmetricMatrix<F>
where
    Vec<Vec<F>>: CanonicalDeserialize,
{
    fn check(&self) -> Result<(), SerializationError> {
        self.0.check()
    }
}

impl<F: Clone> CanonicalDeserialize for SymmetricMatrix<F>
where
    Vec<Vec<F>>: CanonicalDeserialize,
{
    fn deserialize_with_mode<R: Read>(
        reader: R,
        compress: Compress,
        validate: Validate,
    ) -> Result<Self, SerializationError> {
        Vec::<Vec<F>>::deserialize_with_mode(reader, compress, validate).map(Self)
    }
}

// impl<F: Clone> ToBytes for SymmetricMatrix<F>
// where
//     Vec<Vec<F>>: ToBytes,
// {
//     type ToBytesError = <Vec<Vec<F>> as ToBytes>::ToBytesError;
//
//     fn to_bytes(&self) -> Result<Vec<u8>, Self::ToBytesError> {
//         self.0.to_bytes()
//     }
// }
//
// impl<F: Clone> FromBytes for SymmetricMatrix<F>
// where
//     Vec<Vec<F>>: FromBytes,
// {
//     type FromBytesError = <Vec<Vec<F>> as FromBytes>::FromBytesError;
//
//     fn from_bytes(bytes: &[u8]) -> Result<Self, Self::FromBytesError> {
//         Vec::<Vec<F>>::from_bytes(bytes).map(Self)
//     }
// }
