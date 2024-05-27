//! This adds a few utility functions for serializing and deserializing
//! [arkworks](http://arkworks.rs/) types that implement [CanonicalSerialize] and [CanonicalDeserialize].
//! Adapted from [o1-labs/proof-systems](https://raw.githubusercontent.com/o1-labs/proof-systems/31c76ceae3122f0ce09cded8260960ed5cbbe3d8/utils/src/serialization.rs).

use std::fmt::Debug;

use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, SerializationError};
use serde;
use serde_with::Bytes;

pub mod ser {
    use serde_with::{DeserializeAs, SerializeAs};

    use super::*;

    pub fn serialize<S>(val: &impl CanonicalSerialize, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        let mut bytes = vec![];
        val.serialize_uncompressed(&mut bytes)
            .map_err(serde::ser::Error::custom)?;

        Bytes::serialize_as(&bytes, serializer)
    }

    pub fn deserialize<'de, T, D>(deserializer: D) -> Result<T, D::Error>
    where
        T: CanonicalDeserialize,
        D: serde::Deserializer<'de>,
    {
        let bytes: Vec<u8> = Bytes::deserialize_as(deserializer)?;
        T::deserialize_uncompressed(&mut &bytes[..]).map_err(serde::de::Error::custom)
    }
}

pub struct SerdeAs;

impl<T> serde_with::SerializeAs<T> for SerdeAs
where
    T: CanonicalSerialize,
{
    fn serialize_as<S>(val: &T, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        let mut bytes = vec![];
        val.serialize_uncompressed(&mut bytes)
            .map_err(serde::ser::Error::custom)?;

        Bytes::serialize_as(&bytes, serializer)
    }
}

impl<'de, T> serde_with::DeserializeAs<'de, T> for SerdeAs
where
    T: CanonicalDeserialize,
{
    fn deserialize_as<D>(deserializer: D) -> Result<T, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        let bytes: Vec<u8> = Bytes::deserialize_as(deserializer)?;
        T::deserialize_uncompressed(&mut &bytes[..]).map_err(serde::de::Error::custom)
    }
}

pub trait ToBytes {
    type ToBytesError: Debug;
    fn to_bytes(&self) -> Result<Vec<u8>, Self::ToBytesError>;
}

pub trait FromBytes: Sized {
    type FromBytesError: Debug;
    fn from_bytes(bytes: &[u8]) -> Result<Self, Self::FromBytesError>;
}

// impl<T: Serialize> ToBytes for T {
//     type ToBytesError = bincode::Error;
//     fn to_bytes(&self) -> Result<Vec<u8>, Self::ToBytesError> {
//         bincode::serialize(self)
//     }
// }
//
// impl<T: for<'de> Deserialize<'de>> FromBytes for T {
//     type FromBytesError = bincode::Error;
//     fn from_bytes(bytes: &[u8]) -> Result<Self, bincode::Error> {
//         bincode::deserialize(bytes)
//     }
// }

impl<T: CanonicalSerialize> ToBytes for T {
    type ToBytesError = SerializationError;
    fn to_bytes(&self) -> Result<Vec<u8>, Self::ToBytesError> {
        let mut bytes = vec![];
        self.serialize_compressed(&mut bytes)?;
        Ok(bytes)
    }
}

impl<T: CanonicalDeserialize> FromBytes for T {
    type FromBytesError = SerializationError;
    fn from_bytes(bytes: &[u8]) -> Result<Self, Self::FromBytesError> {
        T::deserialize_compressed(bytes)
    }
}

// impl ToBytes for u64 {
//     type ToBytesError = ();
//     fn to_bytes(&self) -> Result<Vec<u8>, Self::ToBytesError> {
//         Ok(self.to_be_bytes().to_vec())
//     }
// }

// impl FromBytes for u64 {
//     type FromBytesError = ();
//     fn from_bytes(bytes: &[u8]) -> Result<Self, Self::FromBytesError> {
//         if bytes.len() != 8 {
//             return Err(());
//         }
//         let mut buf = [0u8; 8];
//         buf.copy_from_slice(&bytes[..8]);
//         Ok(u64::from_be_bytes(buf))
//     }
// }
//
// impl<T: ToBytes> ToBytes for Vec<T> {
//     type ToBytesError = T::ToBytesError;
//     fn to_bytes(&self) -> Result<Vec<u8>, Self::ToBytesError> {
//         let mut bytes = vec![];
//         bytes.extend_from_slice(&(self.len() as u64).to_bytes()?);
//         for elem in self {
//             bytes.extend_from_slice(&elem.to_bytes()?);
//         }
//         Ok(bytes)
//     }
// }
//
// impl<T: FromBytes> FromBytes for Vec<T> {
//     type FromBytesError = T::FromBytesError;
//     fn from_bytes(bytes: &[u8]) -> Result<Self, Self::FromBytesError> {
//         let mut bytes = bytes;
//         let len = u64::from_bytes(bytes)?;
//         let mut vec = Vec::with_capacity(len as usize);
//         for _ in 0..len {
//             let (elem_bytes, rest) = T::from_bytes(bytes)?;
//             vec.push(elem_bytes);
//             bytes = rest;
//         }
//         if !bytes.is_empty() {
//             return Err(());
//         }
//         Ok(vec)
//     }
// }
