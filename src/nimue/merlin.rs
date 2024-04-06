use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, Compress};
use bincode;
use bincode::Options;
use nalgebra::Scalar;
use nimue::{ByteReader, DuplexHash, IOPatternError, Merlin};
use num_traits::Zero;
use serde::{Deserialize, Serialize};

use crate::lattice_arithmetic::matrix::{Matrix, SymmetricMatrix, Vector};

pub trait SerMerlin<H>
    where H: DuplexHash<u8>,
          Self: ByteReader
{
    fn next_serializable_size<S>(&mut self, size: usize) -> Result<S, IOPatternError>
        where S: Serialize + for<'de> Deserialize<'de>
    {
        let mut buf = vec![0u8; size];
        self.fill_next_bytes(&mut buf)?;
        let res = bincode::deserialize(buf.as_slice());
        Ok(res.expect("Invalid"))
    }
    fn next_like_serializable<S>(&mut self, like: &S) -> Result<S, IOPatternError>
        where S: Serialize + for<'de> Deserialize<'de>
    {
        let size = bincode::serialized_size(&like).unwrap();
        self.next_serializable_size(size as usize)
    }

    fn next_canonical_serializable_size<S: CanonicalSerialize + CanonicalDeserialize>(&mut self, size: usize) -> Result<S, IOPatternError> {
        let mut buf = vec![0u8; size];
        self.fill_next_bytes(&mut buf)?;
        let res = S::deserialize_compressed(buf.as_slice());
        Ok(res.expect("Invalid"))
    }

    fn next_like_canonical_serializable<S: CanonicalSerialize + CanonicalDeserialize>(&mut self, like: &S) -> Result<S, IOPatternError> {
        self.next_canonical_serializable_size(like.serialized_size(Compress::Yes))
    }

    fn next_vector<F: Scalar + CanonicalSerialize + CanonicalDeserialize + Zero>(&mut self, n: usize) -> Result<Vector<F>, IOPatternError> {
        let size = F::zero().serialized_size(Compress::Yes);
        Ok(Vector::<F>::from_fn(n, |_, _| self.next_canonical_serializable_size(size).unwrap()))
    }

    fn next_symmetric_matrix<F: Zero + Clone>(&mut self, size: usize) -> Result<SymmetricMatrix<F>, IOPatternError>
        where SymmetricMatrix<F>: Serialize + for<'de> Deserialize<'de>
    {
        self.next_like_serializable(&SymmetricMatrix::<F>::zero(size))
    }

    fn next_matrix<F: Scalar + Zero + CanonicalSerialize + CanonicalDeserialize>(&mut self, m: usize, n: usize) -> Result<Matrix<F>, IOPatternError> {
        let size = F::zero().serialized_size(Compress::Yes);
        Ok(Matrix::<F>::from_fn(m, n, |_, _| self.next_canonical_serializable_size(size).unwrap()))
    }
}

impl<H> SerMerlin<H> for Merlin<'_, H, u8> where H: DuplexHash<u8> {}