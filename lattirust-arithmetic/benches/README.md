# Split NTT plans

Run from the workspace root:

```sh
cargo +nightly-2026-09-17 bench --locked -p lattirust-arithmetic --bench partial_ntt
```

The benchmark compares three paths with the same coefficient backend and compiler. `Reference` is the dense transform and per-call Vandermonde inversion from revision `83341df`, retained in `support/partial_ntt_reference.rs`. `Root API` uses the existing root-taking methods and constructs a plan on each call. `Reused plan` receives a previously constructed `PartialNttPlan`.

Each value is the median of nine batches. The benchmark calibrates each batch to approximately 5 ms, with at most 16,384 iterations. Each transform includes the same input-array copy. Plan construction is excluded from the reused-plan measurements and measured separately. Inputs, roots, plans, and outputs pass through `std::hint::black_box`. These are single-thread transform measurements, not proof-system measurements.

## Measurements

Measured on 2026-11-21 on an Apple M2 Max, arm64 macOS, using Rust `nightly-2026-09-17` and the Cargo bench profile. Times are in microseconds. Small timings are sensitive to scheduling and clock frequency.

| Q | N | D | Factors | Operation | Reference | Root API | Reused plan | Reference / plan |
| ---: | ---: | ---: | ---: | :--- | ---: | ---: | ---: | ---: |
| 17 | 16 | 8 | 2 | Forward | 0.178 | 0.564 | 0.034 | 5.21× |
| 17 | 16 | 8 | 2 | Inverse | 0.237 | 0.478 | 0.031 | 7.67× |
| 3329 | 256 | 2 | 128 | Forward | 105.654 | 4.755 | 1.605 | 65.84× |
| 3329 | 256 | 2 | 128 | Inverse | 5093.708 | 4.713 | 1.577 | 3230.96× |
| 12289 | 256 | 1 | 256 | Forward | 397.352 | 8.434 | 2.466 | 161.13× |
| 12289 | 256 | 1 | 256 | Inverse | 39713.625 | 8.595 | 2.450 | 16206.83× |
| 274877906837 | 256 | 128 | 2 | Forward | 0.846 | 1.215 | 0.528 | 1.60× |
| 274877906837 | 256 | 128 | 2 | Inverse | 1.118 | 1.205 | 0.543 | 2.06× |

| Q | N | Factors | Plan construction (µs) | Plan storage (bytes) |
| ---: | ---: | ---: | ---: | ---: |
| 17 | 16 | 2 | 0.805 | 104 |
| 3329 | 256 | 128 | 3.178 | 6152 |
| 12289 | 256 | 256 | 5.908 | 12296 |
| 274877906837 | 256 | 2 | 0.665 | 104 |

The speedups include an algorithm change: the power-of-two path replaces dense evaluation and matrix inversion with twisted radix-two transforms. Both transform directions take O(N log(Factors)) operations, apart from the linear scaling and permutation steps. Plan construction and storage take O(Factors). Reuse also avoids plan construction; rebuilding a plan can cost more than the original transform when there are few factors.

The plan contains only public transform constants. It has no polynomial or witness storage, heap allocations, or global cache. Existing root-taking signatures and factor order are preserved. Non-power-of-two split orders continue to use the dense implementation. Callers must use the same primitive root for a transform, its inverse, and multiplication in the transform domain.

Validation: all 480 arithmetic unit tests and the plan documentation example passed. The added tests compare factor remainders with direct polynomial evaluation, compare multiplication with coefficient-domain products, cover one, two, 128, and 256 factors, reject invalid parameters and roots, and check the existing APIs at split order three.
