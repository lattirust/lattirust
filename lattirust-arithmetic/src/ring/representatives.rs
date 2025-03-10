use std::fmt::Debug;

use num_traits::Signed;

pub trait WithSignedRepresentative: Sized + Clone {
    type SignedRepresentative: Signed
        + From<Self>
        + Into<Self>
        + TryFrom<i128, Error: Debug>
        + Clone
        + Debug
        + Send
        + Sync;

    fn as_signed_representative(&self) -> Self::SignedRepresentative {
        (*self).clone().into()
    }
}
