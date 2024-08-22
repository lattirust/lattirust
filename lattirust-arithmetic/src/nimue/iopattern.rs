use ark_serialize::{CanonicalSerialize, Compress};
use nalgebra::Scalar;
use nimue::{ByteIOPattern, DuplexHash, IOPattern};
use num_traits::Zero;

use crate::linear_algebra::{Matrix, SymmetricMatrix, Vector};
use crate::nimue::serialization::ToBytes;
use crate::traits::FromRandomBytes;

pub trait SerIOPattern
where
    Self: ByteIOPattern + Sized,
{
    fn absorb_serializable_like<S: ToBytes + Sized>(self, like: &S, label: &'static str) -> Self {
        let s = like.to_bytes().unwrap().len();
        self.add_bytes(s, label)
    }

    fn absorb_canonical_serializable_like<S: CanonicalSerialize>(
        self,
        like: &S,
        label: &'static str,
    ) -> Self {
        let s = like.serialized_size(Compress::Yes);
        self.add_bytes(s, label)
    }

    fn absorb_vec<S>(self, size: usize, label: &'static str) -> Self
    where
        S: Zero + Clone,
        Vec<S>: ToBytes,
    {
        let s = vec![S::zero(); size].to_bytes().unwrap().len();
        self.add_bytes(s, label)
    }

    fn absorb_vector<S: Scalar + Clone + Zero>(self, size: usize, label: &'static str) -> Self
    where
        Vector<S>: ToBytes,
    {
        self.absorb_serializable_like(&Vector::<S>::zeros(size), label)
    }

    fn absorb_vectors<S: Scalar + Clone + Zero>(
        self,
        size: usize,
        num_vectors: usize,
        label: &'static str,
    ) -> Self
    where
        Vec<Vector<S>>: ToBytes,
    {
        self.absorb_serializable_like(&vec![Vector::<S>::zeros(size); num_vectors], label)
    }

    fn absorb_vector_canonical<S: Scalar + Clone + Zero>(
        self,
        size: usize,
        label: &'static str,
    ) -> Self
    where
        Vector<S>: CanonicalSerialize,
    {
        self.absorb_canonical_serializable_like(&Vector::<S>::zeros(size), label)
    }

    fn absorb_symmetric_matrix<S: Clone + Zero>(self, size: usize, label: &'static str) -> Self
    where
        SymmetricMatrix<S>: ToBytes,
    {
        self.absorb_serializable_like(&SymmetricMatrix::<S>::zero(size), label)
    }

    fn absorb_matrix<S: Scalar + Zero>(
        self,
        num_rows: usize,
        num_cols: usize,
        label: &'static str,
    ) -> Self
    where
        Matrix<S>: CanonicalSerialize,
    {
        self.absorb_canonical_serializable_like(&Matrix::<S>::zeros(num_rows, num_cols), label)
    }

    fn absorb_matrix_ser<S: Scalar + Zero>(
        self,
        num_rows: usize,
        num_cols: usize,
        label: &'static str,
    ) -> Self
    where
        Matrix<S>: ToBytes,
    {
        self.absorb_serializable_like(&Matrix::<S>::zeros(num_rows, num_cols), label)
    }
}

impl<H: DuplexHash<u8>> SerIOPattern for IOPattern<H> {}

pub trait SqueezeFromRandomBytes
where
    Self: ByteIOPattern + Sized,
{
    fn squeeze_elem<T, A: FromRandomBytes<T>>(self, label: &'static str) -> Self {
        self.challenge_bytes(A::byte_size(), label)
    }

    fn squeeze_vec<T, A: FromRandomBytes<T>>(self, size: usize, label: &'static str) -> Self {
        self.challenge_bytes(size * A::byte_size(), label)
    }

    fn squeeze_vector<T, A: FromRandomBytes<T>>(self, size: usize, label: &'static str) -> Self {
        self.challenge_bytes(size * A::byte_size(), label)
    }

    fn squeeze_vectors<T, A: FromRandomBytes<T>>(
        self,
        size: usize,
        n_vectors: usize,
        label: &'static str,
    ) -> Self {
        self.challenge_bytes(n_vectors * size * A::byte_size(), label)
    }

    fn squeeze_matrix<T, A: FromRandomBytes<T>>(
        self,
        num_rows: usize,
        num_cols: usize,
        label: &'static str,
    ) -> Self {
        self.challenge_bytes(num_rows * num_cols * A::byte_size(), label)
    }

    fn squeeze_matrices<T, A: FromRandomBytes<T>>(
        self,
        num_rows: usize,
        num_cols: usize,
        num_matrices: usize,
        label: &'static str,
    ) -> Self {
        self.challenge_bytes(num_matrices * num_rows * num_cols * A::byte_size(), label)
    }

    fn squeeze_binary_matrix(self, nrows: usize, ncols: usize, label: &'static str) -> Self {
        self.challenge_bytes(nrows * ncols, label)
    }
}

impl<H: DuplexHash<u8>> SqueezeFromRandomBytes for IOPattern<H> {}

// NOTE: In nimue, ratchet() is not exposed through a trait. This trait allows other traits to call ratchet as part of their implementation.
pub trait RatchetIOPattern {
    fn ratchet(self) -> Self;
}

impl<H: DuplexHash<u8>> RatchetIOPattern for IOPattern<H> {
    fn ratchet(self) -> Self {
        self.ratchet()
    }
}
