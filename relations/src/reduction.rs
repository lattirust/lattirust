use nimue::{Arthur, Merlin, ProofResult};
use crate::traits::Relation;

pub trait Reduction<RelationIn: Relation, RelationOut: Relation> {
    fn prove(crs: RelationIn::PublicParameters, instance: RelationIn::Instance, witness: RelationIn::Witness, merlin: &mut Merlin) -> ProofResult<(RelationOut::Instance, RelationOut::Witness)>;
    fn verify(crs: RelationIn::PublicParameters, instance: RelationIn::Instance, arthur: &mut Arthur) -> ProofResult<RelationOut::Instance>;
}