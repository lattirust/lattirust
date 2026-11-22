# LattiRust PCS

Research implementations of the [Hachi](https://eprint.iacr.org/2026/156) and [Greyhound](https://eprint.iacr.org/2024/1293) commitment and opening reductions. This crate uses the current LattiRust arithmetic and principal-relation types. It has no dependency on a particular closing argument.

```sh
cargo test -p lattirust-pcs
```

## Composition through `Reduction`

The PCS stages implement the existing `relations::Relation` and `relations::reduction::Reduction` traits. There is no separate closer trait.

| Reduction | Input relation | Output relation |
| :--- | :--- | :--- |
| `HachiToGreyhound<N>` | Original field polynomial commitment and evaluation | Successor Greyhound table commitment and evaluation |
| `GreyhoundToShort<R>` | Greyhound table commitment and evaluations | Short affine equations with an integer squared-norm bound |
| `ShortToPrincipal<R>` | Short affine relation | Existing `PrincipalRelation<R>` with zero padding |
| `labrador::reduction::LabradorRound<S>` | Principal relation | Folded principal relation |
| `labrador::zk::subring::IntoDegree256<B,N>` | Short affine relation, degree `N` dividing 256 | Degree-256 principal relation, same modulus |
| `labrador::zk::reduction::HidingCompression<R>` | Linked principal relation | Hiding-folded linked principal relation |
| `labrador::zk::reduction::LnpLite<R>` | Linked principal relation | Masked terminal principal relation |

Public keys are shared through `Arc`; reduction handoffs do not copy their matrices. No witness is stored in a key or retained between calls. The prover returns the output index, instance and private witness. The verifier returns the output index and instance. Framework tests compare the public outputs and check that the output witness satisfies the output relation. The final principal relation still needs an argument; returning a reduced relation is not verifier acceptance.

`HachiOpening<N>` uses `HachiIndex<N>` for verifier-selected keys and public dimensions, `HachiClaim<N>` for the commitment, evaluation and context, and base-field coefficients as the witness. `GreyhoundOpening<R>` uses a table witness and canonical unsigned radix digits for its commitment. Its generic implementation requires a prime coefficient field. The lower-level functions still accept an explicit decomposition function. `ShortAffine<R>` checks the existing affine rows and the integer norm bound.

Hachi and Greyhound retain their existing SHAKE challenge derivation. Their `Merlin`/`Arthur` implementations transport fixed-shape field messages. Decoder dimensions come from public parameters, not proof-supplied length headers. `ShortToPrincipal` sends no messages. A LaBRADOR round uses its native `Merlin`/`Arthur` transcript through `Reduction` and retains the existing PCS continuation binding. This is not a new common challenge derivation across all stages.

The companion `labrador::pcs::LabradorBackend` provides the complete argument as ordinary convenience methods. In a coordinated checkout with `lattirust/` and `labrador/` as siblings:

```rust,ignore
let proof = argument.prove_opening(&index, &claim, &coefficients)?;
argument.verify_opening(&index, &claim, &proof)?;
```

```sh
cargo test --manifest-path ../labrador/Cargo.toml --features pcs
```

The Cargo feature selects whether to build the PCS argument. It does not select a different composition interface. The core PCS crate has no dependency on LaBRADOR. The existing LaBRADOR entry points retain their challenge distributions; the PCS argument selects coefficient-wise ternary projections, fixed-budget unbiased prime and ring challenges, and the integer-filtered degree-64 fold sampler.

The lower-level `hachi_initial::prepare_hachi_reference` and `reduce_hachi` expose the Hachi handoff. `greyhound::commit`, `prepare_opening` and `reduce` expose the single-table protocol. `greyhound::grouped::reduce_grouped` remains a verifier reduction without a grouped prover helper.

## Protocol scope

Both field-to-ring maps retain the imported general trace construction over the full point. They do not use the fixed-basis combination of extension-valued partial evaluations described in the source patch's discussion of Hachi Sections 3.2 and 4.5. The quartic profile packs `D/4` base coefficients per word. Hachi Figure 7's endpoints are batched with a line restriction. SHAKE256 transcripts use domain separation and bind the inputs from which Equation (20) is derived. These encodings are not the reference implementations' wire formats.

The LaBRADOR argument requires verifier-fixed CRS matrices and an explicit nonempty round schedule. It checks canonical message encodings, coefficient-wise JL aggregation, the final equations, and an integer norm bound. Its terminal message discloses the final folded witness. **This closing argument is not zero knowledge.** The companion `zk::argument::Parameters` now composes hiding compression, LNP-Lite, and recursive LaBRADOR for a principal relation. Its `zk::subring::IntoDegree256` helpers also prove short affine relations at the same modulus. See the companion `src/zk/argument.rs` and `src/zk/subring.rs` for the interfaces and parameter limits. Hachi/Greyhound prefixes still send unmasked sumchecks, so using a ZK argument for the final affine relation alone does not make the full PCS opening zero knowledge.

The integration test uses deliberately small correctness parameters, including a real LaBRADOR fold. Its matrices, ranks, modulus, challenge weights, and norm bounds do not establish a security level. `prepare_opening_reference` is not an optimized large-instance prover. Large profiles in `tools/` are analytical circuit inventories, not measured runtimes or calibrated parameters. The Python scripts are optional analysis tools; they are not part of proving or verification.

## Provenance

Imported from the user-supplied `lattirust-hachi-greyhound.patch`, exported from the uncommitted Cloud branch `hachi-greyhound`. The patch applies to LattiRust commit `2aa8f46`. Integration updates arkworks 0.5 to the current 0.6 types, uses the existing reduction framework, and reuses the existing LaBRADOR implementation.

The framework refactor preserves the pre-refactor LaBRADOR closing bytes in the integrated arkworks 0.6 fixture. The new Hachi/Greyhound transport encoding omits redundant vector length headers and is not compatible with the former structured proof container. The original patch and its protocol equations remain the source of the imported construction; this refactor adds software interfaces and tests, not a new security theorem.
