use criterion::BatchSize::PerIteration;
use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use humansize::DECIMAL;
use log::info;
use nimue::IOPattern;

use lattirust_arithmetic::ring::Z2_64;
use lova::prover::Prover;
use lova::util::OptimizationMode::{OptimizeForSpeed, OptimizeForSpeedWithCompletenessError};
use lova::util::{
    rand_matrix_with_bounded_column_norms, BaseRelation, Instance, LovaIOPattern, OptimizationMode,
    PublicParameters,
};
use lova::verifier::Verifier;
use relations::traits::Relation;

type F = Z2_64;
const SECURITY_PARAMETER: usize = 128;
const LOG_FIAT_SHAMIR: usize = 64;

const WITNESS_SIZES: [usize; 3] = [1 << 17, 1 << 18, 1 << 19];

fn pretty_print(param: f64) -> String {
    format!("{param} = 2^{}", param.log2())
}

pub fn criterion_benchmark(c: &mut Criterion) {
    env_logger::builder().is_test(true).try_init().unwrap();

    let mut group = c.benchmark_group("lova");
    for witness_size in WITNESS_SIZES {
        for mode in [OptimizeForSpeedWithCompletenessError, OptimizeForSpeed] {
            let pp =
                PublicParameters::<F>::new(witness_size, mode, SECURITY_PARAMETER, LOG_FIAT_SHAMIR);

            info!(
                "Theoretical proof size for N={}, mode={}: {}",
                pretty_print(witness_size as f64),
                mode,
                humansize::format_size(pp.proof_size_bytes(), DECIMAL)
            );

            let witness_1 = rand_matrix_with_bounded_column_norms(
                pp.witness_len(),
                pp.inner_security_parameter,
                pp.norm_bound as i128,
            );
            let instance_1 = Instance::new(&pp, &witness_1);
            debug_assert!(BaseRelation::is_satisfied(&pp, &instance_1, &witness_1));

            let witness_2 = rand_matrix_with_bounded_column_norms(
                pp.witness_len(),
                pp.inner_security_parameter,
                pp.norm_bound as i128,
            );
            let instance_2 = Instance::new(&pp, &witness_2);
            debug_assert!(BaseRelation::is_satisfied(&pp, &instance_2, &witness_2));

            let proof = &mut vec![];
            group.bench_with_input(
                BenchmarkId::new(format!("prover_{}", mode), witness_size),
                &(pp.clone(), witness_1, witness_2),
                |b, (pp, witness_1, witness_2)| {
                    b.iter_batched(
                        || {
                            (
                                IOPattern::new("lova").fold(&pp).to_merlin(),
                                witness_1.clone(),
                                witness_2.clone(),
                            )
                        },
                        |(mut merlin, witness_1, witness_2)| {
                            // Prove folding
                            let new_witness =
                                Prover::fold(&mut merlin, &pp, witness_1, witness_2).unwrap();
                            black_box(new_witness);
                            let folding_proof = merlin.transcript();

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
                "Actual proof size for N={}, mode={}:      {}",
                pretty_print(witness_size as f64),
                mode,
                humansize::format_size(proof.len(), DECIMAL)
            );

            group.bench_with_input(
                BenchmarkId::new(format!("verifier_{}", mode), witness_size),
                &(pp, instance_1, instance_2),
                |b, (pp, instance_1, instance_2)| {
                    b.iter_batched(
                        || {
                            (
                                IOPattern::new("lova").fold(&pp).to_arthur(proof),
                                instance_1.clone(),
                                instance_2.clone(),
                            )
                        },
                        |(mut arthur, instance_1, instance_2)| {
                            // Verify folding
                            let new_instance =
                                Verifier::fold(&mut arthur, &pp, instance_1, instance_2).unwrap();
                            black_box(new_instance);
                        },
                        PerIteration,
                    )
                },
            );
        }
    }
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
