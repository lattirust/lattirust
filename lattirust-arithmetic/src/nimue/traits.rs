use crate::challenge_set::binary::BinaryChallengeSet;
use ark_std::rand::{CryptoRng, RngCore};
use nimue::{Arthur, ByteChallenges, DuplexHash, IOPatternError, Merlin};

use crate::linear_algebra::Matrix;
use crate::linear_algebra::Scalar;
use crate::linear_algebra::Vector;
use crate::ring::Zq;
use crate::traits::FromRandomBytes;

pub trait ChallengeFromRandomBytes
where
    Self: ByteChallenges,
{
    fn challenge<T, A: FromRandomBytes<T>>(&mut self) -> Result<T, IOPatternError> {
        let mut bytes = vec![0u8; A::byte_size()];
        match self.fill_challenge_bytes(&mut bytes) {
            Ok(_) => A::try_from_random_bytes(bytes.as_slice()).ok_or(IOPatternError::from(
                "error while generating ring element from random bytes",
            )),
            Err(s) => Err(s),
        }
    }

    fn challenge_vec<T, A: FromRandomBytes<T>>(
        &mut self,
        size: usize,
    ) -> Result<Vec<T>, IOPatternError> {
        let mut vals = Vec::<T>::with_capacity(size);
        for _ in 0..size {
            match self.challenge::<T, A>() {
                Ok(v) => vals.push(v),
                Err(s) => return Err(s),
            }
        }
        Ok(vals)
    }

    fn challenge_vector<T: Scalar, A: FromRandomBytes<T>>(
        &mut self,
        size: usize,
    ) -> Result<Vector<T>, IOPatternError> {
        self.challenge_vec::<T, A>(size)
            .map(|v| Vector::<T>::from_vec(v))
    }

    fn challenge_vectors<T: Scalar, A: FromRandomBytes<T>>(
        &mut self,
        size: usize,
        n_vectors: usize,
    ) -> Result<Vec<Vector<T>>, IOPatternError> {
        let mut vals = Vec::<Vector<T>>::with_capacity(n_vectors);
        for _ in 0..n_vectors {
            match self.challenge_vector::<T, A>(size) {
                Ok(val) => vals.push(val),
                Err(e) => return Err(e),
            }
        }
        Ok(vals)
    }

    fn challenge_matrix<T: Scalar, A: FromRandomBytes<T>>(
        &mut self,
        n_rows: usize,
        n_cols: usize,
    ) -> Result<Matrix<T>, IOPatternError> {
        self.challenge_vec::<T, A>(n_rows * n_cols)
            .map(|v| Matrix::<T>::from_vec(n_rows, n_cols, v))
    }

    fn challenge_matrices<T: Scalar, A: FromRandomBytes<T>>(
        &mut self,
        n_rows: usize,
        n_cols: usize,
        n_matrices: usize,
    ) -> Result<Vec<Matrix<T>>, IOPatternError> {
        let mut vals = Vec::<Matrix<T>>::with_capacity(n_matrices);
        for _ in 0..n_matrices {
            match self.challenge_matrix::<T, A>(n_rows, n_cols) {
                Ok(val) => vals.push(val),
                Err(e) => return Err(e),
            }
        }
        Ok(vals)
    }

    fn challenge_binary_matrix(
        &mut self,
        n_rows: usize,
        n_cols: usize,
    ) -> Result<Matrix<Zq<2>>, IOPatternError> {
        // TODO: use more efficient sampling, this is currently using one byte where one bit would be enough
        self.challenge_matrix::<Zq<2>, BinaryChallengeSet<Zq<2>>>(n_rows, n_cols)
    }
}

impl<H: DuplexHash<u8>> ChallengeFromRandomBytes for Arthur<'_, H, u8> {}

impl<H: DuplexHash<u8>, R: RngCore + CryptoRng> ChallengeFromRandomBytes for Merlin<H, u8, R> {}
