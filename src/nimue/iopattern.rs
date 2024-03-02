use ark_serialize::{CanonicalSerialize, Compress};
use bincode::{DefaultOptions, Options};
use derive_more::Deref;
use nimue::{DefaultHash, DefaultRng, DuplexHash, IOPattern};
use nimue::hash::Unit;
use num_traits::Zero;
use serde::Serialize;

use crate::lattice_arithmetic::matrix::Vector;
use crate::lattice_arithmetic::poly_ring::PolyRing;
use crate::lattice_arithmetic::traits::FromRandomBytes;
use crate::nimue::arthur::LatticeArthur;
use crate::nimue::merlin::LatticeMerlin;

#[derive(Deref)]
pub struct LatticeIOPattern<R, H = DefaultHash, U = u8>
    where
        R: PolyRing,
        H: DuplexHash<U>,
        U: Unit,
{
    #[deref]
    io: IOPattern<H, U>,
    _base: std::marker::PhantomData<R>,
}

impl<R, H, U> LatticeIOPattern<R, H, U> where
    R: PolyRing,
    H: DuplexHash<U>,
    U: Unit,
{
    pub fn new(label: &'static str) -> Self {
        Self { io: IOPattern::new(label).into(), _base: std::marker::PhantomData::default() }
    }

    pub fn to_arthur(&self) -> LatticeArthur<R, H, DefaultRng, U> {
        LatticeArthur::new(self, DefaultRng::default())
    }

    pub fn to_merlin<'a>(&self, transcript: &'a [u8]) -> LatticeMerlin<'a, R, H, U> {
        LatticeMerlin::<R, H, U>::new(self, transcript)
    }

    pub fn add_bytes(self, count: usize, label: &'static str) -> Self {
        self.io.absorb(count, label).into()
    }

    pub fn challenge_bytes(self, count: usize, label: &'static str) -> Self {
        self.io.squeeze(count, label).into()
    }

    pub fn ratchet(self) -> Self {
        self.io.ratchet().into()
    }

    pub fn as_bytes(&self) -> &[u8] {
        self.io.as_bytes()
    }

    pub fn absorb_serializable_like<S: Serialize + Sized>(self, like: &S, label: &'static str) -> Self {
        let s = bincode::serialized_size(&like).unwrap();
        self.add_bytes(s as usize, label)
    }

    pub fn absorb_canonical_serializable_like<S: CanonicalSerialize>(self, like: &S, label: &'static str) -> Self {
        let s = like.serialized_size(Compress::Yes);
        self.add_bytes(s, label)
    }

    pub fn absorb_ringelem(self, label: &'static str) -> Self {
        self.absorb_serializable_like(&R::zero(), label)
    }

    pub fn absorb_vec(self, size: usize, label: &'static str) -> Self {
        self.absorb_serializable_like(&vec![R::zero(); size], label)
    }

    pub fn absorb_vector(self, size: usize, label: &'static str) -> Self {
        self.absorb_serializable_like(&Vector::<R>::zeros(size), label)
    }

    pub fn absorb_vectors(self, size: usize, num_vectors: usize, label: &'static str) -> Self {
        self.absorb_serializable_like(&vec![&Vector::<R>::zeros(size); num_vectors], label)
    }

    pub fn absorb_lower_triangular_matrix(self, n: usize, label: &'static str) -> Self {
        let mat = (0..n).map(
            |i| (0..i + 1).map(
                |_| R::zero()
            ).collect::<Vec<R>>()
        ).collect::<Vec<Vec<R>>>();
        self.absorb_serializable_like(&mat, label)
    }

    pub fn absorb_baseringelem(self, label: &'static str) -> Self {
        self.absorb_canonical_serializable_like(&R::BaseRing::zero(), label)
    }

    pub fn absorb_vec_baseringelem(self, size: usize, label: &'static str) -> Self {
        self.absorb_canonical_serializable_like(&vec![R::BaseRing::zero(); size], label)
    }

    pub fn absorb_vector_baseringelem(self, size: usize, label: &'static str) -> Self {
        // TODO: we use Vec here instead of Vector, because Vector does not implement CanonicalSerialize, which is not ideal
        self.absorb_vec_baseringelem(size, label)
    }

    pub fn squeeze_elem<T, A: FromRandomBytes<T>>(self, label: &'static str) -> Self {
        self.challenge_bytes(A::byte_size(), label)
    }

    pub fn squeeze_vec<T, A: FromRandomBytes<T>>(self, size: usize, label: &'static str) -> Self {
        self.challenge_bytes(size * A::byte_size(), label)
    }

    pub fn squeeze_vector<T, A: FromRandomBytes<T>>(self, size: usize, label: &'static str) -> Self {
        self.challenge_bytes(size * A::byte_size(), label)
    }

    pub fn squeeze_vectors<T, A: FromRandomBytes<T>>(self, size: usize, n_vectors: usize, label: &'static str) -> Self {
        self.challenge_bytes(n_vectors * size * A::byte_size(), label)
    }

    pub fn squeeze_matrix<T, A: FromRandomBytes<T>>(self, num_rows: usize, num_cols: usize, label: &'static str) -> Self {
        self.challenge_bytes(num_rows * num_cols * A::byte_size(), label)
    }

    pub fn squeeze_matrices<T, A: FromRandomBytes<T>>(self, num_rows: usize, num_cols: usize, num_matrices: usize, label: &'static str) -> Self {
        self.challenge_bytes(num_matrices * num_rows * num_cols * A::byte_size(), label)
    }

    pub fn squeeze_binary_matrix(self, num_rows: usize, num_cols: usize, label: &'static str) -> Self {
        self.challenge_bytes((num_rows * num_cols).div_ceil(8), label)
    }
}


impl<R, H, U> From<IOPattern<H, U>> for LatticeIOPattern<R, H, U> where
    R: PolyRing,
    H: DuplexHash<U>,
    U: Unit,
{
    fn from(value: IOPattern<H, U>) -> Self {
        Self { io: value, _base: std::marker::PhantomData::default() }
    }
}