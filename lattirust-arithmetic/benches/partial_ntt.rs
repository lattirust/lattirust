//! Compare the reference transform, compatible wrappers, and a reused plan.
use std::{hint::black_box, time::Instant};

use ark_std::{
    rand::{rngs::StdRng, SeedableRng},
    UniformRand,
};
use lattirust_arithmetic::{
    partial_ntt::{PartialNTT, PartialNttPlan},
    ring::{CyclotomicPolyRingSplittedNTT, Zq1},
};

#[path = "support/partial_ntt_reference.rs"]
mod reference;

struct Reference;
impl<const Q: u64, const N: usize, const D: usize, const Z: usize, const F: usize>
    reference::PartialNTT<Q, N, D, Z, F> for Reference
{
}

/// Median of nine batches, with approximately five milliseconds per batch.
fn measure(mut operation: impl FnMut()) -> f64 {
    let start = Instant::now();
    operation();
    let elapsed = start.elapsed().as_secs_f64().max(1e-9);
    let iterations = (0.005 / elapsed).clamp(1.0, 16_384.0) as usize;
    let mut samples = Vec::with_capacity(9);
    for _ in 0..9 {
        let start = Instant::now();
        for _ in 0..iterations {
            operation();
        }
        samples.push(start.elapsed().as_secs_f64() * 1e6 / iterations as f64);
    }
    samples.sort_by(f64::total_cmp);
    samples[4]
}

fn benchmark<const Q: u64, const N: usize, const D: usize, const Z: usize, const F: usize>() {
    type Split<const Q: u64, const N: usize, const D: usize, const Z: usize, const F: usize> =
        CyclotomicPolyRingSplittedNTT<Q, N, D, Z, F>;
    let root = Split::<Q, N, D, Z, F>::primitive_zth_root();
    let plan = PartialNttPlan::<Q, N, D, Z, F>::new(root).unwrap();
    let mut rng = StdRng::seed_from_u64(123);
    let input: [Zq1<Q>; N] = core::array::from_fn(|_| Zq1::<Q>::rand(&mut rng));
    let mut transformed = input;
    plan.ntt(&mut transformed);
    let mut original = input;
    <Reference as reference::PartialNTT<Q, N, D, Z, F>>::ntt(&mut original, root);
    assert_eq!(original, transformed);
    plan.intt(&mut original);
    assert_eq!(original, input);

    let setup = measure(|| {
        black_box(PartialNttPlan::<Q, N, D, Z, F>::new(black_box(root)).unwrap());
    });
    let reference_forward = measure(|| {
        let mut values = black_box(input);
        <Reference as reference::PartialNTT<Q, N, D, Z, F>>::ntt(&mut values, black_box(root));
        black_box(values);
    });
    let wrapper_forward = measure(|| {
        let mut values = black_box(input);
        <Split<Q, N, D, Z, F> as PartialNTT<Q, N, D, Z, F>>::ntt(&mut values, black_box(root));
        black_box(values);
    });
    let planned_forward = measure(|| {
        let mut values = black_box(input);
        black_box(&plan).ntt(&mut values);
        black_box(values);
    });
    let reference_inverse = measure(|| {
        let mut values = black_box(transformed);
        <Reference as reference::PartialNTT<Q, N, D, Z, F>>::intt(&mut values, black_box(root));
        black_box(values);
    });
    let wrapper_inverse = measure(|| {
        let mut values = black_box(transformed);
        <Split<Q, N, D, Z, F> as PartialNTT<Q, N, D, Z, F>>::intt(&mut values, black_box(root));
        black_box(values);
    });
    let planned_inverse = measure(|| {
        let mut values = black_box(transformed);
        black_box(&plan).intt(&mut values);
        black_box(values);
    });
    println!(
        "q={Q}, N={N}, D={D}, factors={F}, plan_bytes={}",
        std::mem::size_of_val(&plan)
    );
    println!("plan_setup_us={setup:.3}");
    println!("operation reference_us wrapper_us reused_plan_us reference_to_plan_speedup");
    println!(
        "forward {reference_forward:.3} {wrapper_forward:.3} {planned_forward:.3} {:.2}",
        reference_forward / planned_forward
    );
    println!(
        "inverse {reference_inverse:.3} {wrapper_inverse:.3} {planned_inverse:.3} {:.2}",
        reference_inverse / planned_inverse
    );
}

fn main() {
    benchmark::<17, 16, 8, 4, 2>();
    benchmark::<3329, 256, 2, 256, 128>();
    benchmark::<12289, 256, 1, 512, 256>();
    benchmark::<274877906837, 256, 128, 4, 2>();
}
