
use std::ops::Deref;

use bincode;
use delegate::delegate;
use nimue::{DefaultHash, DuplexHash, InvalidTag, Merlin};
use nimue::hash::Unit;
use serde::{Deserialize, Serialize};

use crate::lattice_arithmetic::matrix::{Matrix, Vector};
use crate::lattice_arithmetic::poly_ring::PolyRing;
use crate::lattice_arithmetic::ring::Ring;
use crate::lattice_arithmetic::traits::FromRandomBytes;
use crate::nimue::iopattern::LatticeIOPattern;

pub struct LatticeMerlin<'a, H = DefaultHash, U = u8>
    where
        H: DuplexHash<U>,
        U: Unit
{
    merlin: Merlin<'a, H, U>,
}


impl<'a, H, U> LatticeMerlin<'a, H, U>
    where
        H: DuplexHash<U>, U: Unit
{
    pub fn new<PR: PolyRing>(io_pattern: &LatticeIOPattern<PR, H, U>, transcript: &'a [u8]) -> Self {
        Self { merlin: Merlin::<H, U>::new(&io_pattern.deref(), transcript) }
    }

    delegate! {
        to self.merlin {
            pub fn public_input(&mut self, input: &[U]) -> Result<(), InvalidTag>;
            pub fn fill_challenges(&mut self, input: &mut [U]) -> Result<(), InvalidTag>;
            pub fn ratchet(&mut self) -> Result<(), InvalidTag>;
            pub fn preprocess(self) -> Result<&'static [U], InvalidTag>;
        }
    }
}

impl<'a, H> LatticeMerlin<'a, H, u8>
    where H: DuplexHash<u8>
{
    pub fn fill_next_elems<A: Ring>(&mut self, output: &mut [A]) -> Result<(), InvalidTag> {
        let mut buf = vec![0u8; A::byte_size()];
        for o in output.iter_mut() {
            self.merlin.fill_next(&mut buf)?;
            *o = bincode::deserialize(&buf).expect("Invalid");
        }
        Ok(())
    }

    pub fn fill_next_challenges<A: Ring>(&mut self, output: &mut [A]) -> Result<(), InvalidTag> {
        let mut buf = vec![0u8; A::byte_size()];
        for o in output.iter_mut() {
            self.merlin.fill_challenges(&mut buf)?;
            *o = A::from_bytes(&buf).unwrap(); // TODO: or from_random_bytes? or something else?
            // *o = F::from_be_bytes_mod_order(&buf);
        }
        Ok(())
    }

    pub fn next<const N: usize>(&mut self) -> Result<[u8; N], InvalidTag> {
        let mut output = [0u8; N];
        self.merlin.fill_next(&mut output).map(|()| output)
    }

    pub fn next_elems<A: Ring, const N: usize>(&mut self) -> Result<[A; N], InvalidTag> {
        let mut output = [A::default(); N];
        self.fill_next_elems(&mut output).map(|()| output)
    }


    fn next_like<S>(&mut self, like: &S) -> Result<S, InvalidTag>
        where S: Serialize + for<'de> Deserialize<'de>
    {
        let s = bincode::serialized_size(&like).unwrap();
        let mut buf = vec![0u8; s as usize];
        self.merlin.fill_next(&mut buf)?;
        let res = bincode::deserialize(buf.as_slice());
        Ok(res.expect("Invalid"))
    }

    pub fn next_vec<A: Ring>(&mut self, n: usize) -> Result<Vec<A>, InvalidTag> {
        self.next_like(&vec![A::default(); n])
    }

    pub fn next_symmetric_matrix<A: Ring>(&mut self, n: usize) -> Result<Vec<Vec<A>>, InvalidTag> {
        let mut res = Vec::<Vec<A>>::with_capacity(n);
        for i in 0..n {
            res.push(self.next_vec(i)?);
        }
        Ok(res)
    }

    pub fn next_vector<A: Ring>(&mut self, n: usize) -> Result<Vector<A>, InvalidTag> {
        self.next_like(&Vector::<A>::zeros(n))
    }

    pub fn next_matrix<A: Ring>(&mut self, m: usize, n: usize) -> Result<Matrix<A>, InvalidTag> {
        self.next_like(&Matrix::<A>::zeros(m, n))
    }

    pub fn next_vectors<A: Ring>(&mut self, n: usize, num_vectors: usize) -> Result<Vec<Vector<A>>, InvalidTag> {
        self.next_like(&vec![Vector::<A>::zeros(n); num_vectors])
    }

    pub fn next_matrices<A: Ring>(&mut self, m: usize, n: usize, num_matrices: usize) -> Result<Vec<Matrix<A>>, InvalidTag> {
        self.next_like(&vec![Matrix::<A>::zeros(m, n); num_matrices])
    }

    // pub fn challenge<const N: usize>(&mut self) -> Result<[u8; N], InvalidTag> {
    //     let mut output = [0u8; N];
    //     self.fill_challenges(&mut output).map(|()| output)
    // }

    fn challenge<S: Serialize + Deserialize<'a>, C: FromRandomBytes<S>>(&mut self) -> Result<S, InvalidTag> {
        let mut buf = vec![0u8; C::byte_size()];
        self.merlin.fill_challenges(&mut buf)?;
        C::from_random_bytes(&buf).ok_or(InvalidTag::from("error while generating ring element from random bytes"))
    }

    pub fn challenge_vec<A: Ring, T: FromRandomBytes<A>>(&mut self, size: usize) -> Result<Vec<A>, InvalidTag> {
        let mut vals = Vec::<A>::with_capacity(size);
        for _ in 0..size {
            match self.challenge::<A, T>() {
                Ok(v) => vals.push(v),
                Err(s) => return Err(InvalidTag::from(s))
            }
        }
        Ok(vals)
    }

    pub fn challenge_vector<A: Ring, T: FromRandomBytes<A>>(&mut self, size: usize) -> Result<Vector<A>, InvalidTag> {
        self.challenge_vec::<A, T>(size).map(|v| Vector::<A>::from_vec(v))
    }

    pub fn challenge_vectors<A: Ring, T: FromRandomBytes<A>>(&mut self, size: usize, num_vectors: usize) -> Result<Vec<Vector<A>>, InvalidTag> {
        let mut vals = Vec::<Vector<A>>::with_capacity(size);
        for _ in 0..num_vectors {
            match self.challenge_vector::<A, T>(size) {
                Ok(v) => vals.push(v),
                Err(s) => return Err(InvalidTag::from(s))
            }
        }
        Ok(vals)
    }

    pub fn challenge_matrix<A: Ring, T: FromRandomBytes<A>>(&mut self, n_rows: usize, n_cols: usize) -> Result<Matrix<A>, InvalidTag> {
        self.challenge_vec::<A, T>(n_rows * n_cols).map(|v| Matrix::<A>::from_vec(n_rows, n_cols, v))
    }

    pub fn challenge_matrices<A: Ring, T: FromRandomBytes<A>>(&mut self, n_rows: usize, n_cols: usize, n_matrices: usize) -> Result<Vec<Matrix<A>>, InvalidTag> {
        let mut vals = Vec::<Matrix<A>>::with_capacity(n_matrices);
        for _ in 0..n_matrices {
            match self.challenge_matrix::<A, T>(n_rows, n_cols) {
                Ok(val) => vals.push(val),
                Err(e) => return Err(InvalidTag::from(e))
            }
        }
        Ok(vals)
    }
}