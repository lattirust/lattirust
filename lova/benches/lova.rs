use ark_std::test_rng;
use criterion::BatchSize::PerIteration;
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use nimue::IOPattern;

use lattirust_arithmetic::ring::Z2_64;
use lova::prover::Prover;
use lova::util::{
    rand_matrix_with_bounded_column_norms, BaseRelation, Instance, LovaIOPattern, OptimizationMode,
    PublicParameters,
};
use lova::verifier::Verifier;
use relations::traits::Relation;

type F = Z2_64;

const WITNESS_SIZES: [usize; 5] = [1 << 4, 1 << 8, 1 << 12, 1 << 16, 1 << 20];

pub fn criterion_benchmark(c: &mut Criterion) {
    env_logger::builder().is_test(true).try_init().unwrap();
    let rng = &mut test_rng();

    let mut group = c.benchmark_group("lova");
    for witness_size in WITNESS_SIZES {
        let pp = PublicParameters::<F>::new(witness_size, OptimizationMode::OptimizeForSpeed);

        let witness_1 = rand_matrix_with_bounded_column_norms(
            pp.witness_len(),
            pp.security_parameter,
            rng,
            pp.norm_bound as i128,
        );
        let instance_1 = Instance::new(&pp, &witness_1);
        debug_assert!(BaseRelation::is_satisfied(&pp, &instance_1, &witness_1));

        let witness_2 = rand_matrix_with_bounded_column_norms(
            pp.witness_len(),
            pp.security_parameter,
            rng,
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
                    },
                    PerIteration,
                )
            },
        );
    }

    // let new_witness = Prover::fold(&mut arthur, &pp, witness_1.clone(), witness_2.clone()).unwrap();
    // let folding_proof = arthur.transcript();
    //
    // c.bench_function("verifier", |b| {
    //     b.iter(|| {
    //         let instance_1 = instance_1.clone();
    //         let instance_2 = instance_2.clone();
    //
    //
    //     })
    // });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
