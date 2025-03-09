use std::array::from_fn;
use std::io::{Read, Write};

use ark_serialize::{
    CanonicalDeserialize, CanonicalSerialize, Compress, SerializationError, Valid, Validate,
};
use nalgebra::allocator::Allocator;
use nalgebra::{ArrayStorage, Const, DefaultAllocator, Dim, Dyn, RawStorage, Scalar, VecStorage};
use serde::{Deserialize, Deserializer, Serialize, Serializer};

use crate::linear_algebra::generic_matrix::GenericMatrix;

/// Serialize dynamically-sized/fixed-sized matrix/vector/row-vector
impl<T: Scalar, R: Dim, C: Dim, S: RawStorage<T, R, C>> Serialize for GenericMatrix<T, R, C, S>
where
    nalgebra::Matrix<T, R, C, S>: Serialize,
{
    fn serialize<Ser>(&self, serializer: Ser) -> Result<Ser::Ok, Ser::Error>
    where
        Ser: Serializer,
    {
        self.0.serialize(serializer)
    }
}

/// Deserialize dynamically-sized/fixed-sized matrix/vector/row-vector
impl<'de, T: Scalar, R: Dim, C: Dim, S: RawStorage<T, R, C>> Deserialize<'de>
    for GenericMatrix<T, R, C, S>
where
    nalgebra::Matrix<T, R, C, S>: Deserialize<'de>,
{
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        nalgebra::Matrix::<T, R, C, S>::deserialize(deserializer).map(|x| x.into())
    }
}

/// CanonicalSerialize non-fixed-sized matrix/vector/row-vector
impl<T: Scalar, R: Dim, C: Dim> CanonicalSerialize for GenericMatrix<T, R, C, VecStorage<T, R, C>>
where
    T: CanonicalSerialize,
    VecStorage<T, R, C>: RawStorage<T, R, C>,
    DefaultAllocator: Allocator<Dyn, Dyn>,
{
    fn serialize_with_mode<W: Write>(
        &self,
        mut writer: W,
        compress: Compress,
    ) -> Result<(), SerializationError> {
        let nrows = self.nrows() as u64;
        let ncols = self.ncols() as u64;

        nrows.serialize_with_mode(&mut writer, compress)?;
        ncols.serialize_with_mode(&mut writer, compress)?;
        for item in self.iter() {
            item.serialize_with_mode(&mut writer, compress)?;
        }
        Ok(())
    }

    fn serialized_size(&self, compress: Compress) -> usize {
        8 + 8
            + self
                .iter()
                .map(|x| x.serialized_size(compress))
                .sum::<usize>()
    }
}

/// CanonicalSerialize fixed-sized matrix/vector/row-vector
impl<T: Scalar, const R: usize, const C: usize> CanonicalSerialize
    for GenericMatrix<T, Const<R>, Const<C>, ArrayStorage<T, R, C>>
where
    T: CanonicalSerialize + Send,
    ArrayStorage<T, R, C>: RawStorage<T, Const<R>, Const<C>>,
{
    fn serialize_with_mode<W: Write>(
        &self,
        mut writer: W,
        compress: Compress,
    ) -> Result<(), SerializationError> {
        for item in self.iter() {
            item.serialize_with_mode(&mut writer, compress)?;
        }
        Ok(())
    }

    fn serialized_size(&self, compress: Compress) -> usize {
        self.iter()
            .map(|x| x.serialized_size(compress))
            .sum::<usize>()
    }
}

/// Valid-ate non-fixed-sized matrix/vector/row-vector
impl<T: Scalar, R: Dim, C: Dim> Valid for GenericMatrix<T, R, C, VecStorage<T, R, C>>
where
    T: Valid,
    VecStorage<T, R, C>: RawStorage<T, R, C>,
{
    fn check(&self) -> Result<(), SerializationError> {
        T::batch_check(self.0.as_slice().iter())
    }

    fn batch_check<'a>(
        batch: impl Iterator<Item = &'a Self> + Send,
    ) -> Result<(), SerializationError>
    where
        Self: 'a,
    {
        T::batch_check(batch.flat_map(|x| x.0.as_slice().iter()))
    }
}

/// Valid-ate fixed-sized matrix/vector/row-vector
impl<T: Scalar, const R: usize, const C: usize> Valid
    for GenericMatrix<T, Const<R>, Const<C>, ArrayStorage<T, R, C>>
where
    T: Valid,
    ArrayStorage<T, R, C>: RawStorage<T, Const<R>, Const<C>>,
{
    fn check(&self) -> Result<(), SerializationError> {
        T::batch_check(self.0.as_slice().iter())
    }

    fn batch_check<'a>(
        batch: impl Iterator<Item = &'a Self> + Send,
    ) -> Result<(), SerializationError>
    where
        Self: 'a,
    {
        T::batch_check(batch.flat_map(|x| x.0.as_slice().iter()))
    }
}

/// CanonicalDeserialize non-fixed-sized matrix/vector/row-vector
impl<T: Scalar, R: Dim, C: Dim> CanonicalDeserialize for GenericMatrix<T, R, C, VecStorage<T, R, C>>
where
    T: CanonicalDeserialize + Send,
    VecStorage<T, R, C>: RawStorage<T, R, C>,
    DefaultAllocator: Allocator<Dyn, Dyn>,
{
    fn deserialize_with_mode<Re: Read>(
        mut reader: Re,
        compress: Compress,
        validate: Validate,
    ) -> Result<Self, SerializationError> {
        let nrows = u64::deserialize_with_mode(&mut reader, compress, validate)? as usize;
        let ncols = u64::deserialize_with_mode(&mut reader, compress, validate)? as usize;

        let mut data = Vec::<T>::with_capacity(nrows * ncols);
        for _ in 0..nrows * ncols {
            data.push(T::deserialize_with_mode(&mut reader, compress, validate)?);
        }

        let vec_storage = VecStorage::new(R::from_usize(nrows), C::from_usize(ncols), data);
        // Ok(Self(Self::Inner::from_vec_storage(vec_storage)))
        Ok(Self(Self::Inner::from_data(vec_storage)))
    }
}

/// CanonicalDeserialize fixed-sized matrix/vector/row-vector
impl<T: Scalar, const R: usize, const C: usize> CanonicalDeserialize
    for GenericMatrix<T, Const<R>, Const<C>, ArrayStorage<T, R, C>>
where
    T: CanonicalDeserialize + Send,
    ArrayStorage<T, R, C>: RawStorage<T, Const<R>, Const<C>>,
{
    fn deserialize_with_mode<Re: Read>(
        mut reader: Re,
        compress: Compress,
        validate: Validate,
    ) -> Result<Self, SerializationError> {
        let arr: [[T; R]; C] = from_fn(|_| {
            from_fn(|_| T::deserialize_with_mode(&mut reader, compress, validate).unwrap())
        });
        Ok(Self(Self::Inner::from_data(ArrayStorage(arr))))
    }
}

#[cfg(test)]
mod test {
    use std::fmt::Debug;

    use crate::linear_algebra::{Matrix, SMatrix, SRowVector, Vector};

    use super::*;

    const M: usize = 101;
    const N: usize = 42;

    fn test_canonical_serialization_deserialization<
        T: Scalar,
        R: Dim,
        C: Dim,
        S: nalgebra::RawStorage<T, R, C>,
    >(
        mat: GenericMatrix<T, R, C, S>,
    ) where
        GenericMatrix<T, R, C, S>: CanonicalSerialize + CanonicalDeserialize + PartialEq + Debug,
    {
        for mode in [Compress::No, Compress::Yes] {
            let mut bytes = vec![];
            mat.serialize_with_mode(&mut bytes, mode).unwrap();

            let mat2 = GenericMatrix::<T, R, C, S>::deserialize_with_mode(
                bytes.as_slice(),
                mode,
                Validate::Yes,
            )
            .unwrap();
            assert_eq!(mat, mat2);
        }
    }

    #[test]
    fn test_canonical_serialization_deserialization_dyn_dyn() {
        let rng = &mut ark_std::test_rng();
        let mat = Matrix::<u64>::rand(M, N, rng);
        test_canonical_serialization_deserialization(mat);
    }

    #[test]
    fn test_canonical_serialization_deserialization_const_dyn() {
        let rng = &mut ark_std::test_rng();
        let mat = SRowVector::<u64, M>::rand(rng);
        test_canonical_serialization_deserialization(mat);
    }

    #[test]
    fn test_canonical_serialization_deserialization_dyn_const() {
        let rng = &mut ark_std::test_rng();
        let mat = Vector::<u64>::rand(N, rng);
        test_canonical_serialization_deserialization(mat);
    }

    #[test]
    fn test_canonical_serialization_deserialization_const_const() {
        let rng = &mut ark_std::test_rng();
        let mat = SMatrix::<u64, M, N>::rand(rng);
        test_canonical_serialization_deserialization(mat);
    }
}
