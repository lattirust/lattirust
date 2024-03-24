use std::ops::Deref;

use crypto_bigint::rand_core::{CryptoRng, RngCore};
use delegate::delegate;
use nimue::{Arthur, ByteWriter, DefaultHash, DefaultRng, DuplexHash, IOPatternError, Unit, UnitTranscript};

use crate::labrador::common_reference_string::CommonReferenceString;
use crate::lattice_arithmetic::matrix::{Matrix, Vector};
use crate::lattice_arithmetic::poly_ring::PolyRing;
use crate::lattice_arithmetic::ring::Fq;
use crate::nimue::arthur::SerArthur;
use crate::nimue::iopattern::LatticeIOPattern;
use crate::nimue::traits::ChallengeFromRandomBytes;
use crate::relations::labrador::principal_relation::PrincipalRelation;

impl<PR, H, R> SerArthur<H, R> for LatticeArthur<PR, H, R>
    where
        PR: PolyRing,
        H: DuplexHash<u8>,
        R: RngCore + CryptoRng
{}

impl<PR, H, R> ChallengeFromRandomBytes for LatticeArthur<PR, H, R>
    where
        PR: PolyRing,
        H: DuplexHash<u8>,
        R: RngCore + CryptoRng
{}

pub struct LatticeArthur<PR, H = DefaultHash, R = DefaultRng>
    where
        PR: PolyRing,
        H: DuplexHash<u8>,
        R: RngCore + CryptoRng,
        u8: Unit
{
    arthur: Arthur<H, u8, R>,
    _polyring: std::marker::PhantomData<PR>,
}


impl<PR, H, R> LatticeArthur<PR, H, R>
    where
        PR: PolyRing,
        H: DuplexHash<u8>,
        R: RngCore + CryptoRng,
{
    pub fn new(io_pattern: &LatticeIOPattern<PR, H>, csrng: R) -> Self {
        Self { arthur: Arthur::<H, u8, R>::new(io_pattern.deref(), csrng), _polyring: std::marker::PhantomData }
    }

    delegate! {
        to self.arthur {
            pub fn add_units(&mut self, input: &[u8]) -> Result<(), IOPatternError>;
            pub fn public_units(&mut self, input: &[u8]) -> Result<(), IOPatternError>;
            pub fn fill_challenge_units(&mut self, output: &mut [u8]) -> Result<(), IOPatternError>;
            pub fn ratchet(&mut self) -> Result<(), IOPatternError>;
            pub fn rng(&mut self) -> &mut (impl CryptoRng + RngCore);
            pub fn transcript(&self) -> &[u8];
        }
    }
}

impl<PR, H, R> LatticeArthur<PR, H, R>
    where
        PR: PolyRing,
        H: DuplexHash<u8>,
        R: RngCore + CryptoRng,
{
    pub fn absorb_crs(&mut self, crs: &CommonReferenceString<PR>) -> Result<(), IOPatternError> {
        self.absorb_serializable(&crs)
    }

    pub fn absorb_instance(&mut self, instance: &PrincipalRelation<PR>) -> Result<(), IOPatternError> {
        self.absorb_serializable(&instance)
    }

    pub fn absorb_vec(&mut self, r: &Vec<PR>) -> Result<(), IOPatternError> {
        self.absorb_serializable(&r)
    }

    pub fn absorb_vector(&mut self, r: &Vector<PR>) -> Result<(), IOPatternError> {
        self.absorb_serializable(&r)
    }

    pub fn absorb_vectors(&mut self, r: &Vec<Vector<PR>>) -> Result<(), IOPatternError> {
        self.absorb_serializable(&r)
    }

    pub fn absorb_lower_triangular_matrix(&mut self, r: &Vec<Vec<PR>>) -> Result<(), IOPatternError> {
        self.absorb_serializable(&r)
    }

    pub fn absorb_vec_baseringelem(&mut self, r: &Vec<PR::BaseRing>) -> Result<(), IOPatternError> {
        self.absorb_canonical_serializable(r)
    }

    pub fn absorb_vector_baseringlem(&mut self, r: &Vector<PR::BaseRing>) -> Result<(), IOPatternError> {
        self.absorb_vec_baseringelem(&r.as_slice().to_vec())
    }

    pub fn squeeze_binary_matrix(&mut self, n_rows: usize, n_cols: usize) -> Result<Matrix<Fq<2>>, IOPatternError> {
        let mut vals = Vec::<Fq<2>>::with_capacity(n_cols * n_rows);
        let mut bytes = vec![0u8; (n_cols * n_rows).div_ceil(8)];
        self.fill_challenge_units(&mut bytes)?;

        for i in 0..n_cols * n_rows {
            vals.push(Fq::<2>::from(bytes[i / 8] >> (i % 8) & 1)); // TODO: check
        }

        Ok(Matrix::<Fq<2>>::from_vec(n_rows, n_cols, vals))
    }
}

impl<PR, H, R> ByteWriter for LatticeArthur<PR, H, R>
    where
        PR: PolyRing,
        H: DuplexHash<u8>,
        R: RngCore + CryptoRng {
    fn add_bytes(&mut self, input: &[u8]) -> Result<(), IOPatternError> {
        self.arthur.add_bytes(input)
    }
}

impl<PR, H, R> UnitTranscript<u8> for LatticeArthur<PR, H, R>
    where
        PR: PolyRing,
        H: DuplexHash<u8>,
        R: RngCore + CryptoRng {
    fn public_units(&mut self, input: &[u8]) -> Result<(), IOPatternError> {
        self.public_units(input)
    }

    fn fill_challenge_units(&mut self, output: &mut [u8]) -> Result<(), IOPatternError> {
        self.fill_challenge_units(output)
    }
}