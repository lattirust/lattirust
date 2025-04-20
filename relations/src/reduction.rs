use anyhow::bail;
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

    fn test_completeness(pp: &PublicParameters, size: &RelationIn::Size) -> anyhow::Result<()>
    where
        <RelationOut as Relation>::Instance: std::fmt::Display,
    {
        let (index_in, instance_in, witness_in) = RelationIn::generate_satisfied_instance(size);
        match RelationIn::is_satisfied_err(&index_in, &instance_in, &witness_in) {
            Ok(_) => {}
            Err(e) => {
                bail!("generated instance is not satisfied: {}", e)
            }
        }
        let io = Self::iopattern(pp, &index_in, &instance_in);

        let mut merlin = io.to_merlin();

        let prover_result = Self::prove(pp, &index_in, &instance_in, &witness_in, &mut merlin);

        let (pp_out_prover, instance_out_prover, witness_out) = match prover_result {
            Ok(result) => result,
            Err(e) => {
                bail!("reduction is not complete; prover failed, returned {}", e)
            }
        };

        let sat = RelationOut::is_satisfied_err(&pp_out_prover, &instance_out_prover, &witness_out);
        if sat.is_err() {
            bail!("reduction is not complete; the output witness is not a valid witness for the output instance and public parameters: {}", sat.err().unwrap());
        }

        let proof = merlin.transcript();
        let mut arthur = io.to_arthur(proof);
        let verifier_result = Self::verify(pp, &index_in, &instance_in, &mut arthur);

        let (pp_out_verifier, instance_out_verifier) = match verifier_result {
            Ok(result) => result,
            Err(e) => {
                bail!("reduction is not complete; verifier failed, returned {}", e)
            }
        };

        if pp_out_prover != pp_out_verifier {
            bail!("reduction is not complete; the prover and verifier output different public parameters");
        }

        if instance_out_prover != instance_out_verifier {
            bail!("reduction is not complete; the prover and verifier output different instances");
        }
        Ok(())
    }

    fn test_soundness(pp: &PublicParameters, size: &RelationIn::Size) -> anyhow::Result<()> {
        let (index_in, instance_in, witness_in) = RelationIn::generate_unsatisfied_instance(size);
        debug_assert!(
            !RelationIn::is_satisfied(&index_in, &instance_in, &witness_in),
            "generated instance is satisfied, aborting test"
        );
        let io = Self::iopattern(pp, &index_in, &instance_in);

        let mut merlin = io.to_merlin();

        let prover_result = Self::prove(pp, &index_in, &instance_in, &witness_in, &mut merlin);
        let (_, _, witness_out) = match prover_result {
            Ok(result) => result,
            Err(e) => {
                bail!("Unable to provide meaningful soundness test, the honest prover failed when run on an unsatisfied statement and returned {}", e)
            }
        };

        let proof = merlin.transcript();
        let mut arthur = io.to_arthur(proof);
        let verifier_result = Self::verify(pp, &index_in, &instance_in, &mut arthur);
        match verifier_result {
            Ok((index_out_verifier, instance_out_verifier)) => {
                // verifier accepted the proof, check that the output witness is not valid for the verifier's output instance
                let sat = RelationOut::is_satisfied_err(
                    &index_out_verifier,
                    &instance_out_verifier,
                    &witness_out,
                );
                match sat {
                    Ok(()) => bail!("reduction is not sound; the output witness is a valid witness for the output instance and public parameters"),
                    Err(_) => Ok(()) // verifier accepted the proof, but the new statement is not in the relation
                }
            }
            Err(_) => Ok(()), // verifier rejected the proof
        }
    }
}

#[macro_export]
macro_rules! test_completeness {
    ($Reduction:ty, $pp:expr, $size:expr) => {
        #[test]
        fn test_completeness() {
            <$Reduction>::test_completeness(&$pp, &$size).unwrap()
        }
    };
}

#[macro_export]
macro_rules! test_completeness_with_init {
    ($Reduction:ty, $pp:expr, $size:expr, $init:expr) => {
        #[test]
        fn test_completeness() {
            $init();
            <$Reduction>::test_completeness(&$pp, &$size).unwrap()
        }
    };
}

#[macro_export]
macro_rules! test_soundness {
    ($Reduction:ty, $pp:expr, $size:expr) => {
        #[test]
        fn test_soundness() {
            <$Reduction>::test_soundness(&$pp, &$size).unwrap()
        }
    };
}

#[macro_export]
macro_rules! test_soundness_with_init {
    ($Reduction:ty, $pp:expr, $size:expr, $init:expr) => {
        #[test]
        fn test_soundness() {
            $init();
            <$Reduction>::test_soundness(&$pp, &$size).unwrap()
        }
    };
}
