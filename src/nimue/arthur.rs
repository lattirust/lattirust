use std::ops::Deref;

use bincode;
use crypto_bigint::rand_core::{CryptoRng, RngCore};
use delegate::delegate;
use nimue::{Arthur, DefaultHash, DefaultRng, DuplexHash, IOPatternError, UnitTranscript};
use nimue::hash::Unit;
use crate::labrador::common_reference_string::CommonReferenceString;

use crate::lattice_arithmetic::matrix::{Matrix, Vector};
use crate::lattice_arithmetic::poly_ring::PolyRing;
use crate::lattice_arithmetic::ring::{Ring, Zq};
use crate::lattice_arithmetic::traits::FromRandomBytes;
use crate::nimue::iopattern::LatticeIOPattern;
use crate::relations::labrador::principal_relation::PrincipalRelation;

pub struct LatticeArthur<PR, H = DefaultHash, R = DefaultRng, U = u8>
    where
        PR: PolyRing,
        H: DuplexHash<U>,
        R: RngCore + CryptoRng,
        U: Unit
{
    arthur: Arthur<H, U, R>,
    _polyring: std::marker::PhantomData<PR>,
}


impl<PR, H, R, U> LatticeArthur<PR, H, R, U>
    where
        PR: PolyRing,
        H: DuplexHash<U>,
        R: RngCore + CryptoRng,
        U: Unit
{
    pub fn new(io_pattern: &LatticeIOPattern<PR, H, U>, csrng: R) -> Self {
        Self { arthur: Arthur::<H, U, R>::new(io_pattern.deref(), csrng), _polyring: std::marker::PhantomData }
    }

    delegate! {
        to self.arthur {
        pub fn add_units(&mut self, input: &[U]) -> Result<(), IOPatternError>;
        pub fn public_units(&mut self, input: &[U]) -> Result<(), IOPatternError>;
        pub fn fill_challenge_units(&mut self, output: &mut [U]) -> Result<(), IOPatternError>;
        pub fn ratchet(&mut self) -> Result<(), IOPatternError>;
        pub fn rng(&mut self) -> &mut (impl CryptoRng + RngCore);
        pub fn transcript(&self) -> &[u8];
        }
    }
}

impl<PR, H, R> LatticeArthur<PR, H, R, u8>
    where
        PR: PolyRing,
        H: DuplexHash<u8>,
        R: RngCore + CryptoRng,
{
    pub fn absorb<S: serde::Serialize>(&mut self, msg: &S) -> Result<(), IOPatternError> {
        match bincode::serialize(&msg) {
            Ok(bytes) => self.add_units(bytes.as_slice()),
            Err(e) => Err(IOPatternError::from(e.to_string()))
        }
    }

    pub fn absorb_crs(&mut self, crs: &CommonReferenceString<PR>) -> Result<(), IOPatternError> {
        self.absorb(&crs)
    }

    pub fn absorb_instance(&mut self, instance: &PrincipalRelation<PR>) -> Result<(), IOPatternError> {
        self.absorb(&instance)
    }

    pub fn absorb_vec<A: Ring>(&mut self, r: &Vec<A>) -> Result<(), IOPatternError> {
        self.absorb(&r)
    }

    pub fn absorb_vector<A: Ring>(&mut self, r: &Vector<A>) -> Result<(), IOPatternError> {
        self.absorb(&r)
    }

    pub fn absorb_vectors<A: Ring>(&mut self, r: &Vec<Vector<A>>) -> Result<(), IOPatternError> {
        self.absorb(&r)
    }

    pub fn absorb_lower_triangular_matrix<A: Ring>(&mut self, r: &Vec<Vec<A>>) -> Result<(), IOPatternError> {
        self.absorb(&r)
    }

    pub fn squeeze_elem<T, A: FromRandomBytes<T>>(&mut self) -> Result<T, IOPatternError> {
        let mut bytes = vec![0u8; A::byte_size()];
        match self.fill_challenge_units(&mut bytes) {
            Ok(_) => A::from_random_bytes(bytes.as_slice()).ok_or(IOPatternError::from("error while generating ring element from random bytes")),
            Err(s) => Err(s)
        }
    }

    pub fn squeeze_vec<T, A: FromRandomBytes<T>>(&mut self, size: usize) -> Result<Vec<T>, IOPatternError> {
        let mut vals = Vec::<T>::with_capacity(size);
        for _ in 0..size {
            match self.squeeze_elem::<T, A>() {
                Ok(v) => vals.push(v),
                Err(s) => return Err(s)
            }
        }
        Ok(vals)
    }

    pub fn squeeze_vector<T: Ring, A: FromRandomBytes<T>>(&mut self, size: usize) -> Result<Vector<T>, IOPatternError> {
        self.squeeze_vec::<T, A>(size).map(|v| Vector::<T>::from_vec(v))
    }

    pub fn squeeze_vectors<T: Ring, A: FromRandomBytes<T>>(&mut self, size: usize, n_vectors: usize) -> Result<Vec<Vector<T>>, IOPatternError> {
        let mut vals = Vec::<Vector<T>>::with_capacity(n_vectors);
        for _ in 0..n_vectors {
            match self.squeeze_vector::<T, A>(size) {
                Ok(val) => vals.push(val),
                Err(e) => return Err(e)
            }
        }
        Ok(vals)
    }

    pub fn squeeze_matrix<T: Ring, A: FromRandomBytes<T>>(&mut self, n_rows: usize, n_cols: usize) -> Result<Matrix<T>, IOPatternError> {
        self.squeeze_vec::<T, A>(n_rows * n_cols).map(|v| Matrix::<T>::from_vec(n_rows, n_cols, v))
    }

    pub fn squeeze_matrices<T: Ring, A: FromRandomBytes<T>>(&mut self, n_rows: usize, n_cols: usize, n_matrices: usize) -> Result<Vec<Matrix<T>>, IOPatternError> {
        let mut vals = Vec::<Matrix<T>>::with_capacity(n_matrices);
        for _ in 0..n_matrices {
            match self.squeeze_matrix::<T, A>(n_rows, n_cols) {
                Ok(val) => vals.push(val),
                Err(e) => return Err(e)
            }
        }
        Ok(vals)
    }

    pub fn squeeze_binary_matrix(&mut self, n_rows: usize, n_cols: usize) -> Result<Matrix<Zq<2>>, IOPatternError> {
        let mut vals = Vec::<Zq<2>>::with_capacity(n_cols * n_rows);
        let mut bytes = vec![0u8; (n_cols * n_rows).div_ceil(8)];
        self.fill_challenge_units(&mut bytes)?;

        for i in 0..n_cols * n_rows {
            vals.push(Zq::<2>::from(bytes[i / 8] >> (i % 8) & 1)); // TODO: check
        }

        Ok(Matrix::<Zq<2>>::from_vec(n_rows, n_cols, vals))
    }
}
