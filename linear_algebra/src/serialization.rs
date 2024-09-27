use ark_serialize::{
    CanonicalDeserialize, CanonicalSerialize, Compress, SerializationError, Valid, Validate,
};
use ark_std::io::{Read, Write};
use nalgebra::{
    allocator::Allocator, Const, DefaultAllocator, Dim, Dyn, IsContiguous, RawStorage, Scalar,
    VecStorage,
};
use serde::{Deserialize, Deserializer, Serialize, Serializer};

use crate::{generic_matrix::GenericMatrix, matrix::Matrix, vector::GenericVector};

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

impl<T: Scalar, R: Dim, C: Dim, S: RawStorage<T, R, C>> CanonicalSerialize
    for GenericMatrix<T, R, C, S>
where
    T: CanonicalSerialize,
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

impl<T: Scalar, R: Dim, C: Dim, S: RawStorage<T, R, C> + IsContiguous + Sync> Valid
    for GenericMatrix<T, R, C, S>
where
    T: CanonicalDeserialize,
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

impl<T: Scalar> CanonicalDeserialize for Matrix<T>
where
    T: CanonicalDeserialize + Send,
    DefaultAllocator: Allocator<T, Dyn, Dyn>,
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

        let vec_storage = VecStorage::new(Dyn::from_usize(nrows), Dyn::from_usize(ncols), data);
        Ok(Self(Self::Inner::from_vec_storage(vec_storage)))
    }
}

impl<T: Scalar, R: Dim, S: RawStorage<T, R, Const<1>> + IsContiguous + Sync> CanonicalDeserialize
    for GenericVector<T, R, S>
where
    T: CanonicalDeserialize + Send,
    DefaultAllocator: Allocator<T, Dyn, Const<1>>,
    Self: TryFrom<Vec<T>>,
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
        // let vec_storage = VecStorage::new(Dyn::from_usize(nrows), Const::<1>, data);
        Self::try_from(data).map_err(|_| SerializationError::InvalidData)
    }
}

#[cfg(test)]
mod test {
    use super::*;

    const M: usize = 101;
    const N: usize = 42;
    #[test]
    fn test_canonical_serialization_deserialization() {
        let rng = &mut ark_std::test_rng();
        let mat = Matrix::<u64>::rand(M, N, rng);

        for mode in [Compress::No, Compress::Yes] {
            let mut bytes = vec![];
            mat.serialize_with_mode(&mut bytes, mode).unwrap();

            let mat2 = Matrix::<u64>::deserialize_with_mode(bytes.as_slice(), mode, Validate::Yes)
                .unwrap();
            assert_eq!(mat, mat2);
        }
    }
}
