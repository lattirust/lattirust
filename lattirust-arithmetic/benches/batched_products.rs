//! Compare complete batched operations with the pre-generalization backend.
//! Each sample runs in a fresh process under /usr/bin/time to measure peak RSS.
use std::{
    collections::hash_map::DefaultHasher,
    hash::{Hash, Hasher},
    hint::black_box,
    process::Command,
    time::Instant,
};

use ark_std::{
    rand::{rngs::StdRng, SeedableRng},
    UniformRand,
};
use lattirust_arithmetic::{
    linear_algebra::{inner_products, Matrix, SymmetricMatrix, Vector},
    ring::{Pow2CyclotomicPolyRing, Ring, Zq1},
};
use rayon::prelude::*;

#[allow(dead_code)]
#[path = "support/batched_linear_algebra_reference.rs"]
mod reference_linear_algebra;
#[allow(dead_code)]
#[path = "support/batched_products_reference.rs"]
mod reference_products;

type Polynomial<const Q: u64, const N: usize> = Pow2CyclotomicPolyRing<Zq1<Q>, N>;

#[derive(Clone, Copy)]
enum Operation {
    Dot,
    Gram,
    Cross,
    MatrixVector,
}

fn flatten<R: Ring>(matrix: SymmetricMatrix<R>) -> Vec<R> {
    let rows: Vec<Vec<R>> = matrix.into();
    rows.into_iter().flatten().collect()
}

fn worker<const Q: u64, const N: usize>(
    baseline: bool,
    ordinary: bool,
    operation: Operation,
    rows: usize,
    width: usize,
    small: bool,
) {
    let mut rng = StdRng::seed_from_u64(0x165e1);
    let mut sample = || {
        Polynomial::<Q, N>::from_fn(|_| {
            if small {
                match ark_std::rand::Rng::gen_range(&mut rng, 0..3) {
                    0 => -Zq1::<Q>::ONE,
                    1 => Zq1::<Q>::ZERO,
                    _ => Zq1::<Q>::ONE,
                }
            } else {
                Zq1::<Q>::rand(&mut rng)
            }
        })
    };
    let left: Vec<Vector<Polynomial<Q, N>>> = (0..rows)
        .map(|_| Vector::from_vec((0..width).map(|_| sample()).collect()))
        .collect();
    let right_count = match operation {
        Operation::Gram => 0,
        Operation::MatrixVector | Operation::Dot => 1,
        Operation::Cross => rows,
    };
    let right: Vec<Vector<Polynomial<Q, N>>> = (0..right_count)
        .map(|_| Vector::from_vec((0..width).map(|_| sample()).collect()))
        .collect();
    let matrix = matches!(operation, Operation::MatrixVector)
        .then(|| Matrix::from_fn(rows, width, |i, j| left[i][j]));
    // Matrix operations need only the matrix, not a second copy of its rows.
    let left = if matrix.is_some() { Vec::new() } else { left };
    let compute = || -> Vec<Polynomial<Q, N>> {
        if ordinary {
            return match operation {
                Operation::Dot => vec![left[0].dot(&right[0])],
                Operation::Gram => inner_products::lower_triang_indices(rows)
                    .into_par_iter()
                    .map(|(i, j)| left[i].dot(&left[j]))
                    .collect(),
                Operation::Cross => inner_products::lower_triang_indices(rows)
                    .into_par_iter()
                    .map(|(i, j)| {
                        if i == j {
                            left[i].dot(&right[j])
                        } else {
                            left[i].dot(&right[j]) + left[j].dot(&right[i])
                        }
                    })
                    .collect(),
                Operation::MatrixVector => (matrix.as_ref().unwrap() * &right[0])
                    .iter()
                    .copied()
                    .collect(),
            };
        }
        match (baseline, operation) {
            (true, Operation::Dot) => vec![reference_linear_algebra::dot_product(
                left[0].as_slice(),
                right[0].as_slice(),
            )],
            (false, Operation::Dot) => vec![inner_products::dot_product(
                left[0].as_slice(),
                right[0].as_slice(),
            )],
            (true, Operation::Gram) => flatten(reference_linear_algebra::inner_products(&left)),
            (false, Operation::Gram) => flatten(inner_products::inner_products(&left)),
            (true, Operation::Cross) => flatten(
                reference_linear_algebra::symmetric_cross_products(&left, &right),
            ),
            (false, Operation::Cross) => {
                flatten(inner_products::symmetric_cross_products(&left, &right))
            }
            (true, Operation::MatrixVector) => {
                reference_linear_algebra::matrix_vector_product(matrix.as_ref().unwrap(), &right[0])
                    .iter()
                    .copied()
                    .collect()
            }
            (false, Operation::MatrixVector) => {
                inner_products::matrix_vector_product(matrix.as_ref().unwrap(), &right[0])
                    .iter()
                    .copied()
                    .collect()
            }
        }
    };
    let warmup = black_box(compute());
    let mut hasher = DefaultHasher::new();
    warmup.hash(&mut hasher);
    let checksum = hasher.finish();
    drop(warmup);
    let mut times = Vec::new();
    for _ in 0..3 {
        let start = Instant::now();
        let output = black_box(compute());
        times.push(start.elapsed().as_secs_f64() * 1000.0);
        let mut hasher = DefaultHasher::new();
        output.hash(&mut hasher);
        assert_eq!(checksum, hasher.finish());
        drop(output);
    }
    times.sort_by(f64::total_cmp);
    println!("elapsed_ms={:.6} checksum={checksum}", times[1]);
}

fn run_worker(backend: &str, case: &str) {
    let baseline = match backend {
        "previous" => true,
        "generic" | "ordinary" => false,
        _ => panic!("unknown backend"),
    };
    let ordinary = backend == "ordinary";
    match case {
        "q38-small-gram8x1024" => {
            worker::<274877906837, 256>(baseline, ordinary, Operation::Gram, 8, 1024, true)
        }
        "q38-full-gram8x1024" => {
            worker::<274877906837, 256>(baseline, ordinary, Operation::Gram, 8, 1024, false)
        }
        "q38-full-cross8x1024" => {
            worker::<274877906837, 256>(baseline, ordinary, Operation::Cross, 8, 1024, false)
        }
        "q38-full-dot8192" => {
            worker::<274877906837, 256>(baseline, ordinary, Operation::Dot, 1, 8192, false)
        }
        "q38-full-matvec32x512" => {
            worker::<274877906837, 256>(baseline, ordinary, Operation::MatrixVector, 32, 512, false)
        }
        "q46-full-gram32x128-n256" => {
            worker::<67280421310721, 256>(baseline, ordinary, Operation::Gram, 32, 128, false)
        }
        "q46-full-gram8x128-n256" => {
            worker::<67280421310721, 256>(baseline, ordinary, Operation::Gram, 8, 128, false)
        }
        "q65537-full-gram8x128-n256" => {
            worker::<65537, 256>(baseline, ordinary, Operation::Gram, 8, 128, false)
        }
        "q65537-full-gram4x128-n64" => {
            worker::<65537, 64>(baseline, ordinary, Operation::Gram, 4, 128, false)
        }
        _ => panic!("unknown case"),
    }
}

fn sample(executable: &std::path::Path, backend: &str, case: &str) -> (f64, f64, u64) {
    let mut command = Command::new("/usr/bin/time");
    if cfg!(target_os = "macos") {
        command.arg("-l");
    } else {
        command.arg("-v");
    }
    let output = command
        .arg(executable)
        .args(["--worker", backend, case])
        .env("RAYON_NUM_THREADS", "8")
        .output()
        .expect("launch benchmark worker");
    assert!(
        output.status.success(),
        "worker failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let stdout = String::from_utf8(output.stdout).unwrap();
    let mut fields = stdout.split_whitespace();
    let elapsed = fields
        .next()
        .unwrap()
        .strip_prefix("elapsed_ms=")
        .unwrap()
        .parse()
        .unwrap();
    let checksum = fields
        .next()
        .unwrap()
        .strip_prefix("checksum=")
        .unwrap()
        .parse()
        .unwrap();
    let stderr = String::from_utf8(output.stderr).unwrap();
    let memory_line = stderr
        .lines()
        .find(|line| {
            line.to_ascii_lowercase()
                .contains("maximum resident set size")
        })
        .expect("peak RSS from /usr/bin/time");
    let bytes: f64 = if cfg!(target_os = "macos") {
        memory_line
            .split_whitespace()
            .next()
            .unwrap()
            .parse()
            .unwrap()
    } else {
        memory_line
            .rsplit(':')
            .next()
            .unwrap()
            .trim()
            .parse::<f64>()
            .unwrap()
            * 1024.0
    };
    (elapsed, bytes / (1024.0 * 1024.0), checksum)
}

fn main() {
    let arguments: Vec<_> = std::env::args().collect();
    if arguments.get(1).map(String::as_str) == Some("--worker") {
        run_worker(&arguments[2], &arguments[3]);
        return;
    }
    let executable = std::env::current_exe().unwrap();
    println!("# 8 Rayon threads; 3 fresh processes/backend/case; 3 timed operations/process after warmup");
    println!("case,backend,median_ms,min_ms,max_ms,median_peak_rss_mib,min_peak_rss_mib,max_peak_rss_mib,checksum");
    for case in [
        "q38-small-gram8x1024",
        "q38-full-gram8x1024",
        "q38-full-cross8x1024",
        "q38-full-dot8192",
        "q38-full-matvec32x512",
        "q65537-full-gram4x128-n64",
        "q65537-full-gram8x128-n256",
        "q46-full-gram8x128-n256",
        "q46-full-gram32x128-n256",
    ] {
        if arguments.get(1).map(String::as_str) == Some("--case")
            && arguments.get(2).map(String::as_str) != Some(case)
        {
            continue;
        }
        let mut old = Vec::new();
        let mut new = Vec::new();
        for repeat in 0..3 {
            for backend in if repeat % 2 == 0 {
                ["previous", "generic"]
            } else {
                ["generic", "previous"]
            } {
                let value = sample(&executable, backend, case);
                if backend == "previous" {
                    old.push(value);
                } else {
                    new.push(value);
                }
            }
        }
        assert!(
            old.iter().chain(&new).all(|x| x.2 == old[0].2),
            "different results for {case}"
        );
        for (name, values) in [("previous", old), ("generic", new)] {
            let mut times: Vec<_> = values.iter().map(|x| x.0).collect();
            let mut memory: Vec<_> = values.iter().map(|x| x.1).collect();
            times.sort_by(f64::total_cmp);
            memory.sort_by(f64::total_cmp);
            println!(
                "{case},{name},{:.6},{:.6},{:.6},{:.3},{:.3},{:.3},{}",
                times[1], times[0], times[2], memory[1], memory[0], memory[2], values[0].2
            );
        }
    }
}
