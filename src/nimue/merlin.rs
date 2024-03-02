use std::ops::Deref;

use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, Compress};
use bincode;
use bincode::Options;
use delegate::delegate;
use nalgebra::Scalar;
use nimue::{DefaultHash, DuplexHash, IOPatternError, Merlin, UnitTranscript};
use nimue::hash::Unit;
use serde::{Deserialize, Serialize};

use crate::lattice_arithmetic::matrix::{Matrix, Vector};
use crate::lattice_arithmetic::poly_ring::PolyRing;
use crate::lattice_arithmetic::ring::{Fq, Ring};
use crate::lattice_arithmetic::traits::FromRandomBytes;
use crate::nimue::iopattern::LatticeIOPattern;

pub struct LatticeMerlin<'a, PR, H = DefaultHash, U = u8>
    where
        PR: PolyRing,
        H: DuplexHash<U>,
        U: Unit
{
    merlin: Merlin<'a, H, U>,
    _polyring: std::marker::PhantomData<PR>,
}


impl<'a, PR, H, U> LatticeMerlin<'a, PR, H, U>
    where
        PR: PolyRing,
        H: DuplexHash<U>,
        U: Unit
{
    pub fn new(io_pattern: &LatticeIOPattern<PR, H, U>, transcript: &'a [u8]) -> Self {
        Self { merlin: Merlin::<H, U>::new(&io_pattern.deref(), transcript), _polyring: std::marker::PhantomData }
    }

    delegate! {
        to self.merlin {
            pub fn public_units(&mut self, input: &[U]) -> Result<(), IOPatternError>;
            pub fn fill_challenge_units(&mut self, input: &mut [U]) -> Result<(), IOPatternError>;
            pub fn ratchet(&mut self) -> Result<(), IOPatternError>;
            pub fn preprocess(self) -> Result<&'static [U], IOPatternError>;
        }
    }
}

impl<'a, PR, H> LatticeMerlin<'a, PR, H, u8>
    where
        PR: PolyRing,
        H: DuplexHash<u8>
{
    pub fn fill_next_elems<T: Ring>(&mut self, output: &mut [T]) -> Result<(), IOPatternError> {
        let mut buf = vec![0u8; T::byte_size()];
        for o in output.iter_mut() {
            self.merlin.fill_next_units(&mut buf)?;
            *o = bincode::deserialize(&buf).expect("Invalid");
        }
        Ok(())
    }

    pub fn fill_next_challenges<T: Ring>(&mut self, output: &mut [T]) -> Result<(), IOPatternError> {
        let mut buf = vec![0u8; T::byte_size()];
        for o in output.iter_mut() {
            self.merlin.fill_challenge_units(&mut buf)?;
            *o = T::from_bytes(&buf).unwrap(); // TODO: or from_random_bytes? or something else?
            // *o = F::from_be_bytes_mod_order(&buf);
        }
        Ok(())
    }

    pub fn next<const N: usize>(&mut self) -> Result<[u8; N], IOPatternError> {
        let mut output = [0u8; N];
        self.merlin.fill_next_units(&mut output).map(|()| output)
    }

    pub fn next_elems<T: Ring, const N: usize>(&mut self) -> Result<[T; N], IOPatternError> {
        let mut output = [T::default(); N];
        self.fill_next_elems(&mut output).map(|()| output)
    }

    fn next_like_serializable<S>(&mut self, like: &S) -> Result<S, IOPatternError>
        where S: Serialize + for<'de> Deserialize<'de>
    {
        let s = bincode::serialized_size(&like).unwrap();
        let mut buf = vec![0u8; s as usize];
        self.merlin.fill_next_units(&mut buf)?;
        let res = bincode::deserialize(buf.as_slice());
        Ok(res.expect("Invalid"))
    }

    fn next_like_canonical_serializable<S: CanonicalSerialize + CanonicalDeserialize>(&mut self, like: &S) -> Result<S, IOPatternError> {
        let mut buf = vec![0u8; like.serialized_size(Compress::Yes)];
        self.merlin.fill_next_units(&mut buf)?;
        let res = S::deserialize_compressed(buf.as_slice());
        Ok(res.expect("Invalid"))
    }

    pub fn next_vec(&mut self, n: usize) -> Result<Vec<PR>, IOPatternError> {
        self.next_like_serializable(&vec![PR::default(); n])
    }

    pub fn next_symmetric_matrix(&mut self, n: usize) -> Result<Vec<Vec<PR>>, IOPatternError> {
        let mut res = Vec::<Vec<PR>>::with_capacity(n);
        for i in 0..n {
            res.push(vec![PR::default(); i + 1]);
        }
        self.next_like_serializable(&res)
    }

    pub fn next_vector(&mut self, n: usize) -> Result<Vector<PR>, IOPatternError> {
        self.next_like_serializable(&Vector::<PR>::zeros(n))
    }

    pub fn next_matrix(&mut self, m: usize, n: usize) -> Result<Matrix<PR>, IOPatternError> {
        self.next_like_serializable(&Matrix::<PR>::zeros(m, n))
    }

    pub fn next_vectors(&mut self, n: usize, num_vectors: usize) -> Result<Vec<Vector<PR>>, IOPatternError> {
        self.next_like_serializable(&vec![Vector::<PR>::zeros(n); num_vectors])
    }

    pub fn next_matrices(&mut self, m: usize, n: usize, num_matrices: usize) -> Result<Vec<Matrix<PR>>, IOPatternError> {
        self.next_like_serializable(&vec![Matrix::<PR>::zeros(m, n); num_matrices])
    }

    pub fn next_vec_baseringelem(&mut self, n: usize) -> Result<Vec<PR::BaseRing>, IOPatternError> {
        self.next_like_canonical_serializable(&vec![PR::BaseRing::default(); n])
    }


    pub fn next_vector_baseringelem(&mut self, n: usize) -> Result<Vector<PR::BaseRing>, IOPatternError> {
        self.next_vec_baseringelem(n).map(|v| Vector::<PR::BaseRing>::from_vec(v))
    }

    fn challenge<S, A: FromRandomBytes<S>>(&mut self) -> Result<S, IOPatternError> {
        let mut buf = vec![0u8; A::byte_size()];
        self.merlin.fill_challenge_units(&mut buf)?;
        A::try_from_random_bytes(&buf).ok_or(IOPatternError::from("error while generating ring element from random bytes"))
    }

    pub fn challenge_vec<T, A: FromRandomBytes<T>>(&mut self, size: usize) -> Result<Vec<T>, IOPatternError> {
        let mut vals = Vec::<T>::with_capacity(size);
        for _ in 0..size {
            match self.challenge::<T, A>() {
                Ok(v) => vals.push(v),
                Err(s) => return Err(IOPatternError::from(s))
            }
        }
        Ok(vals)
    }

    pub fn challenge_vector<T: Scalar, A: FromRandomBytes<T>>(&mut self, size: usize) -> Result<Vector<T>, IOPatternError> {
        self.challenge_vec::<T, A>(size).map(|v| Vector::<T>::from_vec(v))
    }

    pub fn challenge_vectors<T: Scalar, A: FromRandomBytes<T>>(&mut self, size: usize, num_vectors: usize) -> Result<Vec<Vector<T>>, IOPatternError> {
        let mut vals = Vec::<Vector<T>>::with_capacity(size);
        for _ in 0..num_vectors {
            match self.challenge_vector::<T, A>(size) {
                Ok(v) => vals.push(v),
                Err(s) => return Err(IOPatternError::from(s))
            }
        }
        Ok(vals)
    }

    pub fn challenge_matrix<T: Scalar, A: FromRandomBytes<T>>(&mut self, n_rows: usize, n_cols: usize) -> Result<Matrix<T>, IOPatternError> {
        self.challenge_vec::<T, A>(n_rows * n_cols).map(|v| Matrix::<T>::from_vec(n_rows, n_cols, v))
    }

    pub fn challenge_matrices<T: Scalar, A: FromRandomBytes<T>>(&mut self, n_rows: usize, n_cols: usize, n_matrices: usize) -> Result<Vec<Matrix<T>>, IOPatternError> {
        let mut vals = Vec::<Matrix<T>>::with_capacity(n_matrices);
        for _ in 0..n_matrices {
            match self.challenge_matrix::<T, A>(n_rows, n_cols) {
                Ok(val) => vals.push(val),
                Err(e) => return Err(IOPatternError::from(e))
            }
        }
        Ok(vals)
    }

    pub fn challenge_binary_matrix(&mut self, n_rows: usize, n_cols: usize) -> Result<Matrix<Fq<2>>, IOPatternError> {
        let mut vals = Vec::<Fq<2>>::with_capacity(n_cols * n_rows);
        let mut bytes = vec![0u8; (n_cols * n_rows).div_ceil(8)];
        self.merlin.fill_challenge_units(&mut bytes)?;
        for i in 0..n_cols * n_rows {
            vals.push(Fq::<2>::from(bytes[i / 8] >> (i % 8) & 1)); // TODO: check
        }

        Ok(Matrix::<Fq<2>>::from_vec(n_rows, n_cols, vals))
    }
}