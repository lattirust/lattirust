pub trait Relation {
    type Instance;
    type Witness;
    type Crs;

    fn is_satisfied(crs: &Self::Crs, x: &Self::Instance, w: &Self::Witness) -> bool;
}