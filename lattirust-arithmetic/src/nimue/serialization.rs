//! This adds a few utility functions for serializing and deserializing
//! [arkworks](http://arkworks.rs/) types that implement [CanonicalSerialize] and [CanonicalDeserialize].
//! Adapted from [o1-labs/proof-systems](https://raw.githubusercontent.com/o1-labs/proof-systems/31c76ceae3122f0ce09cded8260960ed5cbbe3d8/utils/src/serialization.rs).

use std::fmt::Debug;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use serde::{self, Deserialize, Serialize};
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

impl<T: Serialize> ToBytes for T {
    type ToBytesError = bincode::Error;
    fn to_bytes(&self) -> Result<Vec<u8>, Self::ToBytesError> {
        bincode::serialize(self)
    }
}

impl<T: for<'de> Deserialize<'de>> FromBytes for T {
    type FromBytesError = bincode::Error;
    fn from_bytes(bytes: &[u8]) -> Result<Self, bincode::Error> {
        bincode::deserialize(bytes)
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[derive(Serialize, Deserialize, Debug, PartialEq)]
    struct Test {
        a: u32,
        b: u32,
    }

    #[test]
    fn test_serde() {
        let test = Test { a: 1, b: 2 };
        let bytes = test.to_bytes().unwrap();
        let test2 = Test::from_bytes(&bytes).unwrap();
        assert_eq!(test, test2);
    }
}
