use criterion::{BenchmarkId, black_box, Criterion, criterion_group, criterion_main};
use criterion::BatchSize::PerIteration;
use humansize::DECIMAL;
use log::info;
use nimue::IOPattern;

use lattirust_arithmetic::ring::Z2_64;
use lova::prover::Prover;
use lova::util::{
    BaseRelation, Instance, LovaIOPattern, OptimizationMode, PublicParameters,
    rand_matrix_with_bounded_column_norms,
};
use lova::verifier::Verifier;
use relations::traits::Relation;

type F = Z2_64;

const WITNESS_SIZES: [usize; 3] = [1 << 16, 1 << 18, 1 << 20];

pub fn criterion_benchmark(c: &mut Criterion) {
    env_logger::builder().is_test(true).try_init().unwrap();

    let mut group = c.benchmark_group("lova");
    for witness_size in WITNESS_SIZES {
        for opt_mode in [
            OptimizationMode::OptimizeForSpeed,
            OptimizationMode::OptimizeForSize,
        ] {
            info!(
                "Theoretical proof size for N={}, mode={opt_mode}: {}",
                witness_size,
                humansize::format_size(
                    PublicParameters::<F>::proof_size_bytes_with_mode(witness_size, opt_mode),
                    DECIMAL
                )
            );
        }
        let pp = PublicParameters::<F>::new(witness_size, OptimizationMode::OptimizeForSpeed);

        let witness_1 = rand_matrix_with_bounded_column_norms(
            pp.witness_len(),
            pp.security_parameter,
            pp.norm_bound as i128,
        );
        let instance_1 = Instance::new(&pp, &witness_1);
        debug_assert!(BaseRelation::is_satisfied(&pp, &instance_1, &witness_1));

        let witness_2 = rand_matrix_with_bounded_column_norms(
            pp.witness_len(),
            pp.security_parameter,
            pp.norm_bound as i128,
        );
        let instance_2 = Instance::new(&pp, &witness_2);
        debug_assert!(BaseRelation::is_satisfied(&pp, &instance_2, &witness_2));

        let proof = &mut vec![];
        group.bench_with_input(
            BenchmarkId::new("prover", witness_size),
            &(pp.clone(), witness_1, witness_2),
            |b, (pp, witness_1, witness_2)| {
                b.iter_batched(
                    || {
                        (
                            IOPattern::new("lova").fold(&pp).to_arthur(),
                            witness_1.clone(),
                            witness_2.clone(),
                        )
                    },
                    |(mut arthur, witness_1, witness_2)| {
                        // Prove folding
                        let new_witness =
                            Prover::fold(&mut arthur, &pp, witness_1, witness_2).unwrap();
                        black_box(new_witness);
                        let folding_proof = arthur.transcript();

                        // Save proof globally for verifier
                        if proof.len() == 0 {
                            proof.extend_from_slice(folding_proof);
                        }
                    },
                    PerIteration,
                )
            },
        );

        info!(
            "Theoretical proof size for N={}: {}",
            witness_size,
            humansize::format_size(pp.proof_size_bytes(), DECIMAL)
        );
        info!(
            "Actual proof size for N={}:      {}",
            witness_size,
            humansize::format_size(proof.len(), DECIMAL)
        );

        group.bench_with_input(
            BenchmarkId::new("verifier", witness_size),
            &(pp, instance_1, instance_2),
            |b, (pp, instance_1, instance_2)| {
                b.iter_batched(
                    || {
                        (
                            IOPattern::new("lova").fold(&pp).to_merlin(proof),
                            instance_1.clone(),
                            instance_2.clone(),
                        )
                    },
                    |(mut merlin, instance_1, instance_2)| {
                        // Verify folding
                        let new_instance =
                            Verifier::fold(&mut merlin, &pp, instance_1, instance_2).unwrap();
                        black_box(new_instance);
                    },
                    PerIteration,
                )
            },
        );
    }
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
