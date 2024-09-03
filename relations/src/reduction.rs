use std::error::Error;

use nimue::{Arthur, IOPattern, Merlin, ProofResult};

use crate::Relation;

pub trait Reduction<RelationIn: Relation, RelationOut: Relation, PublicParameters>
where
    RelationOut::Index: PartialEq,
    RelationOut::Instance: PartialEq,
{
    type IndexIn = RelationIn::Index;
    type InstanceIn = RelationIn::Instance;
    type WitnessIn = RelationIn::Witness;

    type IndexOut = RelationOut::Index;
    type InstanceOut = RelationOut::Instance;
    type WitnessOut = RelationOut::Witness;

    fn iopattern(
        pp: &PublicParameters,
        index_in: &RelationIn::Index,
        instance_in: &RelationIn::Instance,
    ) -> IOPattern;

    fn prove(
        pp: &PublicParameters,
        index_in: &RelationIn::Index,
        instance: &RelationIn::Instance,
        witness: &RelationIn::Witness,
        merlin: &mut Merlin,
    ) -> ProofResult<(
        RelationOut::Index,
        RelationOut::Instance,
        RelationOut::Witness,
    )>;

    fn verify(
        pp: &PublicParameters,
        index_in: &RelationIn::Index,
        instance_in: &RelationIn::Instance,
        arthur: &mut Arthur,
    ) -> ProofResult<(RelationOut::Index, RelationOut::Instance)>;

    fn test_completeness(
        pp: &PublicParameters,
        size: &RelationIn::Size,
    ) -> Result<(), Box<dyn Error>> {
        let (index_in, instance_in, witness_in) = RelationIn::generate_satisfied_instance(&size);
        debug_assert_eq!(
            RelationIn::is_satisfied(&index_in, &instance_in, &witness_in),
            true,
            "generated instance is not satisfied, aborting test"
        );
        let io = Self::iopattern(&pp, &index_in, &instance_in);

        let mut merlin = io.to_merlin();

        let prover_result = Self::prove(&pp, &index_in, &instance_in, &witness_in, &mut merlin);

        let (pp_out_prover, instance_out_prover, witness_out) = match prover_result {
            Ok(result) => result,
            Err(e) => {
                return Err(From::from(format!(
                    "reduction is not complete; prover failed, returned {}",
                    e
                )))
            }
        };

        if !RelationOut::is_satisfied(&pp_out_prover, &instance_out_prover, &witness_out) {
            return Err("reduction is not complete; the output witness is not a valid witness for the output instance and public parameters".into());
        }

        let proof = merlin.transcript();
        let mut arthur = io.to_arthur(proof);
        let verifier_result = Self::verify(&pp, &index_in, &instance_in, &mut arthur);

        let (pp_out_verifier, instance_out_verifier) = verifier_result
            .map_err(|e| format!("reduction is not complete; verifier failed, returned {e}"))?;

        if pp_out_prover != pp_out_verifier {
            return Err(From::from("reduction is not complete; the prover and verifier output different public parameters"));
        }

        if instance_out_prover != instance_out_verifier {
            return Err(From::from(
                "reduction is not complete; the prover and verifier output different instances",
            ));
        }
        Ok(())
    }

    fn test_soundness(
        pp: &PublicParameters,
        size: &RelationIn::Size,
    ) -> Result<(), Box<dyn Error>> {
        let (index_in, instance_in, witness_in) = RelationIn::generate_unsatisfied_instance(&size);
        debug_assert_eq!(
            !RelationIn::is_satisfied(&index_in, &instance_in, &witness_in),
            false,
            "generated instance is satisfied, aborting test"
        );
        let io = Self::iopattern(&pp, &index_in, &instance_in);

        let mut merlin = io.to_merlin();

        let prover_result = Self::prove(&pp, &index_in, &instance_in, &witness_in, &mut merlin);
        let (pp_out_prover, instance_out_prover, witness_out) = match prover_result {
            Ok(result) => result,
            Err(e) => {
                return Err(From::from(format!(
                    "reduction is not complete; prover failed, returned {}",
                    e
                )))
            }
        };

        if RelationOut::is_satisfied(&pp_out_prover, &instance_out_prover, &witness_out) {
            return Err(From::from("reduction is not sound; the output witness is a valid witness for the output instance and public parameters"));
        }

        let proof = merlin.transcript();
        let mut arthur = io.to_arthur(proof);
        let verifier_result = Self::verify(&pp, &index_in, &instance_in, &mut arthur);

        if verifier_result.is_ok() {
            return Err(From::from(
                "reduction is not sound; the verifier accepted the proof",
            ));
        }

        Ok(())
    }
}

#[macro_export]
macro_rules! test_completeness {
    ($Reduction:tt, $pp:expr, $size:expr) => {
        #[test]
        fn test_completeness() {
            $Reduction::test_completeness(&$pp, &$size).unwrap()
        }
    };
}

#[macro_export]
macro_rules! test_soundness {
    ($Reduction:tt, $pp:expr, $size:expr) => {
        #[test]
        fn test_soundness() {
            $Reduction::test_soundness(&$pp, &$size).unwrap()
        }
    };
}
