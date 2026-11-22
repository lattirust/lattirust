# Generic batched products

## Interface and scope

The existing `dot_product`, `inner_products`, `inner_products2`, `symmetric_cross_products`, and `matrix_vector_product` signatures are unchanged. `Ring` has one additional associated type, `ProductBackend`, with a default backend that declines acceleration. Existing ring implementations can omit this type. Dispatch uses associated types; it does not inspect `TypeId` or downcast inputs.

The public module `linear_algebra::ntt_products` supplies `BatchedProductBackend`, `NoBatchedProducts`, `IntegerCoefficientBackend`, and `NegacyclicProductBackend`. The coefficient adapter `IntegerCoefficient` specifies a single-word modulus and conversion to and from canonical residues. `Zq<C, 1>` supplies this adapter. Coefficient-domain `Pow2CyclotomicPolyRing` forwards batches to its coefficient backend. Other scalar or polynomial representations retain ordinary ring arithmetic unless they supply a backend.

The CRT kernel supports power-of-two degrees from 2 through 2^24 and single-word target moduli. It uses a fixed pool of four auxiliary primes. Unsupported degrees, small batches, coefficient adapters without a single-word modulus, and unsafe integer reconstruction bounds return `None`; existing callers then use ordinary ring products. A direct caller of `products` or `negacyclic_products` must handle `None` in the same way.

## Exact arithmetic and temporary storage

For each requested pair `(i, j)`, the kernel computes `sum_k left[i][k] * right[j][k]` in the original negacyclic ring. If `right` is absent, both operands use `left`. It checks pair dimensions before accepting the batch.

Let `w` be the largest requested row width, `N` the polynomial degree, and `A` and `B` the maximum absolute centered input coefficients. Each signed output coefficient has absolute value at most `w * N * A * B`. The selected auxiliary-prime product `M` must satisfy `2 * w * N * A * B < M`. Every multiplication in this bound is checked for overflow. A batch that exceeds the available reconstruction range is declined.

Transforms use compile-time auxiliary root tables and shared forward/inverse butterfly kernels. A checked `FqConfig::from_canonical_u64` constructor converts already reduced residues with Montgomery multiplication and a constant `R^2 mod q`, avoiding runtime wide division. The original general `from_u64` constructor is unchanged; the comparison baseline retains that conversion path. The batch accumulates pointwise products before taking one inverse transform per output and auxiliary prime. Raw products fit in `u128`: residues are below 2^31, and even `usize::MAX` terms on a 64-bit target have a sum below 2^126.

Each call creates and releases its own transform tiles. The 512-column schedule is unchanged from the prior backend. Scratch storage scales with the number of input rows, polynomial degree, output count, and one tile; no full-row pointer tables are constructed. The backend stores no input coefficients or transforms between calls. It has no witness cache.

## Comparison baseline

The benchmark retains the backend and linear-algebra wrappers present immediately before this change, including the fixed degree-256, modulus-274877906837 dispatch, pointer tables, and 512-column tiles. These were local, uncommitted sources. Their SHA-256 hashes before adaptation were:

```text
ntt_products.rs   2cc9c014f44dfbe8aa0d3b3ba70b5e7b3e34a5c0a02becb7ee8e2f879902ded0
inner_products.rs 76a8a055ff1a1e597cbc96840eb66473da2d7c6f203cad7151c21226b3255bb0
```

The retained copies are in `benches/support/batched_products_reference.rs` and `benches/support/batched_linear_algebra_reference.rs`. Adaptations change imports and module names. The original crate-private coefficient-array copy is reproduced through the public allocation-free conversion into `SVector`. Both backends use the same coefficient arithmetic and shared NTT butterfly kernels in the final binary.

Automatic dispatch is conservative: it advertises the coefficient CRT backend only for target moduli larger than `u32::MAX` and polynomial degrees of at least 128. Benchmarks showed that auxiliary-prime conversion could cost more than ordinary multiplication for smaller moduli or degrees. The explicit `negacyclic_products` function still supports those cases; it uses only the dimension, batch-size, and reconstruction checks above. For automatic dispatch, a batch with a target modulus larger than 2^40 also requires at least eight requested products per supplied input row, counting both input families. Measurements showed that lower reuse could make this path slower than ordinary multiplication. This reuse check runs before input norm scans. The explicit kernel bypasses this additional policy. These cutoffs are conservative performance policy, not requirements for exact reconstruction, and are not a claim of optimal dispatch on all hardware.

For modulus 65537, both previous and final automatic paths use ordinary ring multiplication. The degree-256, modulus-67280421310721 cases check both low-reuse fallback at eight rows and newly available acceleration at 32 rows. The former fixed-ring backend declines both and uses ordinary products.

## Reproduction and measurement

Run from the workspace root:

```sh
cargo +nightly-2026-09-17 bench --locked -p lattirust-arithmetic --bench batched_products
```

Measurements use an Apple M2 Max on arm64 macOS, Rust `nightly-2026-09-17`, the Cargo bench profile, and eight Rayon threads. The driver alternates backend order across three fresh processes per backend and case. Each worker creates deterministic inputs, warms up once, and times three complete operations. Reported time is the median of the three worker medians; the raw CSV also reports their range. Input generation and output hashing are outside the timed interval. Coefficient conversion, transforms, accumulation, inverse transforms, CRT reconstruction, and wrapper allocations are inside it.

Peak RSS comes from `/usr/bin/time -l` for each worker. Reported RSS is the median of three process high-water marks and includes input generation, input storage, thread stacks, warmup, output storage, and temporary arithmetic allocations. It is not a measurement of scratch storage alone. The benchmark also supports Linux `/usr/bin/time -v`. Results from both backends must have equal hashes in every case.

`small` inputs have coefficients sampled uniformly from {-1, 0, 1}. `full` inputs have coefficients sampled uniformly from the target residue field. `Gram r x w` requests all lower-triangular row inner products. `Cross r x w` requests the symmetric cross products, with one term on the diagonal and two off the diagonal. `Matvec r x w` applies an `r`-row matrix to one vector. All cases have degree 256 and modulus 274877906837 unless stated otherwise. `q46` denotes modulus 67280421310721.

## Results

| Case | Previous (ms) | Generic (ms) | Previous RSS (MiB) | Generic RSS (MiB) |
| :--- | ---: | ---: | ---: | ---: |
| Gram 8×1024, small | 343.20 | 586.68 | 26.62 | 26.05 |
| Gram 8×1024, full | 1903.86 | 1325.21 | 27.50 | 27.61 |
| Cross 8×1024, full | 2544.68 | 2434.55 | 48.64 | 49.33 |
| Dot 8192, full | 2601.65 | 2312.42 | 39.53 | 40.03 |
| Matvec 32×512, full | 1928.84 | 2212.42 | 154.84 | 154.38 |
| q=65537, N=64, Gram 4×128 | 2.33 | 2.33 | 3.25 | 3.28 |
| q=65537, N=256, Gram 8×128 | 89.17 | 90.05 | 5.47 | 5.36 |
| q46, Gram 8×128 | 86.59 | 87.99 | 5.48 | 5.52 |
| q46, Gram 32×128 | 1183.77 | 820.40 | 13.58 | 22.55 |

The high-reuse q46 case is 1.44× faster and uses 8.97 MiB more peak RSS. Low-reuse q46 and small-modulus cases retain ordinary multiplication. Results on the previously accelerated ring are mixed: small-coefficient Gram products are 71% slower and matrix-vector products are 15% slower in this run; full-width Gram and dot products are faster. Peak RSS is approximately unchanged on that ring. Several timing ranges overlap, and the three-process samples do not establish a stable speedup across workloads. The small-coefficient and matrix-vector regressions remain performance limitations; the backend should not be described as a uniform optimization.

The [raw CSV](batched-products-results.csv) records all timing and memory ranges. Both q46 rows were remeasured with the final early reuse guard; the other rows use unchanged arithmetic and dispatch paths. The whole batch harness can be rerun, or one case can be selected with `--case CASE` when invoking the built executable.

These measurements cover arithmetic operations. They do not establish an end-to-end prover or verifier speedup. The CSV records process-to-process timing variation, so small differences should not be treated as stable improvements.

## Validation

All 486 arithmetic unit tests and both documentation tests pass. The added tests compare generic batches with independent coefficient-domain schoolbook products across several degrees and moduli. They cover signed and full-width coefficients, ragged but pair-compatible rows, repeated pairs, one through four auxiliary primes, invalid dimensions, zero and small batches, unsupported degrees, integer-bound overflow, and a finite bound larger than the four-prime reconstruction range. Existing linear-algebra entry points are checked against independent products. The checked Montgomery conversion is compared with the original division-based conversion at canonical boundaries, including a modulus near 2^64. The four-prime test has dense coefficients whose signed product exceeds half the three-prime product, so its fourth prime is necessary. The LattiRust workspace also passes `cargo check --locked --release --offline --workspace`.
