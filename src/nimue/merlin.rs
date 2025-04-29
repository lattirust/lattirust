use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use nimue::{ByteReader, DuplexHash, IOPatternError, Merlin};
use num_traits::Zero;

use crate::linear_algebra::{Matrix, Scalar, SymmetricMatrix, Vector};
use crate::nimue::serialization::{FromBytes, ToBytes};

pub trait SerMerlin<H>
where
    H: DuplexHash<u8>,
    Self: ByteReader,
{
    fn err_to_io_pattern_error<E>(e: E) -> IOPatternError
    where
        E: std::fmt::Debug,
    {
        IOPatternError::from(format!("{:?}", e))
    }
    fn next_with_size<S: FromBytes>(&mut self, size: usize) -> Result<S, IOPatternError> {
        let mut buf = vec![0u8; size];
        self.fill_next_bytes(&mut buf)?;
        Ok(S::from_bytes(buf.as_slice()).map_err(Self::err_to_io_pattern_error)?)
    }

    fn next_like<S: FromBytes + ToBytes>(&mut self, like: &S) -> Result<S, IOPatternError> {
        let size = like
            .to_bytes()
            .map_err(Self::err_to_io_pattern_error)?
            .len();
        self.next_with_size(size)
    }

    fn next_canonical_serializable_size<S: CanonicalSerialize + CanonicalDeserialize>(
        &mut self,
        size: usize,
    ) -> Result<S, IOPatternError> {
        let mut buf = vec![0u8; size];
        self.fill_next_bytes(&mut buf)?;
        S::deserialize_compressed(buf.as_slice()).map_err(Self::err_to_io_pattern_error)
    }

    fn next_like_canonical_serializable<S: CanonicalSerialize + CanonicalDeserialize>(
        &mut self,
        like: &S,
    ) -> Result<S, IOPatternError> {
        self.next_canonical_serializable_size(like.compressed_size())
    }

    fn next_vec<F: Zero + Clone>(&mut self, n: usize) -> Result<Vec<F>, IOPatternError>
    where
        Vec<F>: ToBytes + FromBytes,
    {
        self.next_like(&vec![F::zero(); n])
    }
    fn next_vector<F: Scalar + Zero>(&mut self, n: usize) -> Result<Vector<F>, IOPatternError>
    where
        Vector<F>: ToBytes + FromBytes,
    {
        self.next_like(&Vector::<F>::zeros(n))
    }

    fn next_vector_canonical<F: Scalar + Zero>(
        &mut self,
        n: usize,
    ) -> Result<Vector<F>, IOPatternError>
    where
        Vector<F>: CanonicalSerialize + CanonicalDeserialize,
    {
        self.next_like_canonical_serializable(&Vector::<F>::zeros(n))
    }

    fn next_vectors<F: Scalar + Zero>(
        &mut self,
        n: usize,
        num_vectors: usize,
    ) -> Result<Vec<Vector<F>>, IOPatternError>
    where
        Vec<Vector<F>>: ToBytes + FromBytes,
    {
        self.next_like(&vec![Vector::<F>::zeros(n); num_vectors])
    }

    fn next_symmetric_matrix<F: Zero + Clone>(
        &mut self,
        size: usize,
    ) -> Result<SymmetricMatrix<F>, IOPatternError>
    where
        SymmetricMatrix<F>: ToBytes + FromBytes,
    {
        self.next_like(&SymmetricMatrix::<F>::zero(size))
    }

    fn next_matrix<F: Scalar + Zero>(
        &mut self,
        m: usize,
        n: usize,
    ) -> Result<Matrix<F>, IOPatternError>
    where
        Matrix<F>: CanonicalSerialize + CanonicalDeserialize,
    {
        self.next_like_canonical_serializable(&Matrix::<F>::zeros(m, n))
    }

    fn next_matrix_ser<F: Scalar + Zero>(
        &mut self,
        m: usize,
        n: usize,
    ) -> Result<Matrix<F>, IOPatternError>
    where
        Matrix<F>: FromBytes + ToBytes,
    {
        self.next_like(&Matrix::<F>::zeros(m, n))
    }
}

impl<H> SerMerlin<H> for Merlin<'_, H, u8> where H: DuplexHash<u8> {}
