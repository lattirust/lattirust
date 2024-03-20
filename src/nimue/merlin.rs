use std::ops::Deref;

use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, Compress};
use bincode;
use bincode::Options;
use delegate::delegate;
use nimue::{ByteChallenges, BytePublic, ByteReader, ByteTranscript, DefaultHash, DuplexHash, IOPatternError, Merlin, UnitTranscript};
use serde::{Deserialize, Serialize};

use crate::lattice_arithmetic::matrix::{Matrix, Vector};
use crate::lattice_arithmetic::poly_ring::PolyRing;
use crate::lattice_arithmetic::ring::{Fq, Ring};
use crate::lattice_arithmetic::traits::FromRandomBytes;
use crate::nimue::iopattern::LatticeIOPattern;
use crate::nimue::traits::ChallengeFromRandomBytes;

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
}

impl<H> SerMerlin<H> for Merlin<'_, H, u8> where H: DuplexHash<u8> {}

impl<PR, H> SerMerlin<H> for LatticeMerlin<'_, PR, H> where PR: PolyRing, H: DuplexHash<u8> {}


impl<PR: PolyRing, H: DuplexHash<u8>> ChallengeFromRandomBytes for LatticeMerlin<'_, PR, H> {}

pub struct LatticeMerlin<'a, PR, H = DefaultHash>
    where
        PR: PolyRing,
        H: DuplexHash<u8>,
{
    merlin: Merlin<'a, H, u8>,
    _polyring: std::marker::PhantomData<PR>,
}

impl<'a, PR: PolyRing, H: DuplexHash<u8>> ByteReader for LatticeMerlin<'a, PR, H> {
    fn fill_next_bytes(&mut self, input: &mut [u8]) -> Result<(), IOPatternError> {
        self.merlin.fill_next_bytes(input)
    }
}

impl<'a, PR: PolyRing, H: DuplexHash<u8>> BytePublic for LatticeMerlin<'a, PR, H> {
    fn public_bytes(&mut self, input: &[u8]) -> Result<(), IOPatternError> {
        self.merlin.public_bytes(input)
    }
}

impl<'a, PR: PolyRing, H: DuplexHash<u8>> ByteChallenges for LatticeMerlin<'a, PR, H> {
    fn fill_challenge_bytes(&mut self, output: &mut [u8]) -> Result<(), IOPatternError> {
        self.merlin.fill_challenge_bytes(output)
    }
}

impl<'a, PR: PolyRing, H: DuplexHash<u8>> ByteTranscript for LatticeMerlin<'a, PR, H> {}

impl<'a, PR, H> LatticeMerlin<'a, PR, H>
    where
        PR: PolyRing,
        H: DuplexHash<u8>,
{
    pub fn new(io_pattern: &LatticeIOPattern<PR, H>, transcript: &'a [u8]) -> Self {
        Self { merlin: Merlin::<H, u8>::new(&io_pattern.deref(), transcript), _polyring: std::marker::PhantomData }
    }

    delegate! {
        to self.merlin {
            pub fn public_units(&mut self, input: &[u8]) -> Result<(), IOPatternError>;
            pub fn fill_challenge_units(&mut self, input: &mut [u8]) -> Result<(), IOPatternError>;
            pub fn ratchet(&mut self) -> Result<(), IOPatternError>;
            pub fn preprocess(self) -> Result<&'static [u8], IOPatternError>;
        }
    }
}

impl<'a, PR, H> LatticeMerlin<'a, PR, H>
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