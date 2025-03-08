#![feature(associated_type_defaults)]

use std::error::Error;

pub mod ajtai_cm;
pub mod principal_relation;
pub mod r1cs;
pub mod reduction;

pub trait Relation {
    type Size;
    type Index;
    type Instance;
    type Witness;

    /// Returns true iff the index `i` and instance `x` (and witness `w`, if not `None`) are well-defined.
    /// For example, for R1CS, this function should check that the dimensions of the matrices A, B, and C are the same, are consistent with the public parameters, and that the witness has the correct length.
    fn is_well_defined(i: &Self::Index, x: &Self::Instance, w: Option<&Self::Witness>) -> bool;

    /// Returns true iff the index `i` and instance `x` (and witness `w`, if not `None`) are well-defined.
    /// For example, for R1CS, this function should check that the dimensions of the matrices A, B, and C are the same, are consistent with the public parameters, and that the witness has the correct length.
    fn is_well_defined_err(
        i: &Self::Index,
        x: &Self::Instance,
        w: Option<&Self::Witness>,
    ) -> anyhow::Result<()>;

    /// Return true iff the index `i` and instance `x` and witness `w` satisfy the relation.
    /// For example, for R1CS, this function should check that Aw * Bw = Cw.
    fn is_satisfied(i: &Self::Index, x: &Self::Instance, w: &Self::Witness) -> bool;

    /// Return true iff the index `i` and instance `x` and witness `w` satisfy the relation.
    /// For example, for R1CS, this function should check that Aw * Bw = Cw.
    fn is_satisfied_err(
        i: &Self::Index,
        x: &Self::Instance,
        w: &Self::Witness,
    ) -> anyhow::Result<()>;

    /// Generate a (possibly random) instance-witness that are in the relation for a given size `size`.
    /// This is used in particular for testing that `Reduction` implementations are complete (where we require an instance-witness pair in the relation as input).
    fn generate_satisfied_instance(
        size: &Self::Size,
    ) -> (Self::Index, Self::Instance, Self::Witness);

    /// Generate a (possibly random) instance-witness that are not in the relation for a given size `size`.
    /// This is used in particular for testing that `Reduction` implementations are sound (where we require an instance-witness pair not in the relation as input).
    /// The instance may not even be in the language of the relation, in which case the output witness can be any well-defined witness.
    fn generate_unsatisfied_instance(
        size: &Self::Size,
    ) -> (Self::Index, Self::Instance, Self::Witness);
}

#[macro_export]
macro_rules! test_generate_satisfied_instance {
    ($T:tt, $size:expr) => {
        #[test]
        fn test_generate_satisfied_instance() {
            let (index, instance, witness) = <$T>::generate_satisfied_instance(&$size);
            assert!(<$T>::is_well_defined(&index, &instance, Some(&witness)));
            assert!(<$T>::is_satisfied(&index, &instance, &witness));
        }
    };
}

#[macro_export]
macro_rules! test_generate_unsatisfied_instance {
    ($T:tt, $size:expr) => {
        #[test]
        fn test_generate_unsatisfied_instance() {
            let (index, instance, witness) = <$T>::generate_unsatisfied_instance(&$size);
            assert!(<$T>::is_well_defined(&index, &instance, Some(&witness)));
            assert!(!<$T>::is_satisfied(&index, &instance, &witness));
        }
    };
}
