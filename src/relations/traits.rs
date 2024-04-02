pub trait Relation {
    type PublicParameters;
    type Instance;
    type Witness;

    fn is_well_defined(pp: &Self::PublicParameters, x: &Self::Instance, w: Option<&Self::Witness>) -> bool;

    fn is_satisfied(pp: &Self::PublicParameters, x: &Self::Instance, w: &Self::Witness) -> bool;
}