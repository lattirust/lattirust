use std::io::{Read, Write};

use ark_serialize::{Compress, SerializationError, Validate};
use nalgebra::{Const, DefaultAllocator, Dim, Dyn, RawStorage, Scalar, VecStorage};
use nalgebra::allocator::Allocator;
use serde::{Deserialize, Deserializer, Serialize, Serializer};

use crate::linear_algebra::generic_matrix::GenericMatrix;
use crate::linear_algebra::matrix::Matrix;
use crate::linear_algebra::Vector;

impl<T: Scalar, R: Dim, C: Dim, S> Serialize for GenericMatrix<T, R, C, S>
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

impl<'de, T: Scalar, R: Dim, C: Dim, S> Deserialize<'de> for GenericMatrix<T, R, C, S>
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

impl<T: Scalar, R: Dim, C: Dim, S: RawStorage<T, R, C>> ark_serialize::CanonicalSerialize
    for GenericMatrix<T, R, C, S>
where
    T: ark_serialize::CanonicalSerialize,
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

impl<T: Scalar, R: Dim, C: Dim> ark_serialize::Valid for GenericMatrix<T, R, C, VecStorage<T, R, C>>
where
    T: ark_serialize::CanonicalDeserialize,
{
    fn check(&self) -> Result<(), SerializationError> {
        T::batch_check(self.0.data.as_slice().into_iter())
    }

    fn batch_check<'a>(
        batch: impl Iterator<Item = &'a Self> + Send,
    ) -> Result<(), SerializationError>
    where
        Self: 'a,
    {
        T::batch_check(batch.flat_map(|x| x.0.data.as_slice().into_iter()))
    }
}

impl<T: Scalar> ark_serialize::CanonicalDeserialize for Matrix<T>
where
    T: ark_serialize::CanonicalDeserialize + Send,
    DefaultAllocator: Allocator<T, Dyn, Dyn>,
{
    fn deserialize_with_mode<Re: Read>(
        mut reader: Re,
        compress: Compress,
        validate: Validate,
    ) -> Result<Self, SerializationError> {
        let nrows = u64::deserialize_with_mode(&mut reader, compress, validate)? as usize;
        let ncols = u64::deserialize_with_mode(&mut reader, compress, validate)? as usize;

        let data = Vec::<T>::deserialize_with_mode(&mut reader, compress, validate)?;

        if nrows * ncols != data.len() {
            return Err(SerializationError::InvalidData);
        }
        let vec_storage = VecStorage::new(Dyn::from_usize(nrows), Dyn::from_usize(ncols), data);
        Ok(Self {
            0: Self::Inner::from_vec_storage(vec_storage),
        })
    }
}

impl<T: Scalar> ark_serialize::CanonicalDeserialize for Vector<T>
where
    T: ark_serialize::CanonicalDeserialize + Send,
    DefaultAllocator: Allocator<T, Dyn, Const<1>>,
{
    fn deserialize_with_mode<Re: Read>(
        mut reader: Re,
        compress: Compress,
        validate: Validate,
    ) -> Result<Self, SerializationError> {
        let nrows = u64::deserialize_with_mode(&mut reader, compress, validate)? as usize;
        let ncols = u64::deserialize_with_mode(&mut reader, compress, validate)? as usize;

        let data = Vec::<T>::deserialize_with_mode(&mut reader, compress, validate)?;

        if ncols != 1 || nrows * ncols != data.len() {
            return Err(SerializationError::InvalidData);
        }
        let vec_storage = VecStorage::new(Dyn::from_usize(nrows), Const::<1>, data);
        Ok(Self {
            0: Self::Inner::from_vec_storage(vec_storage),
        })
    }
}
