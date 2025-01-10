# Lattirust - A Rust Library for Lattice-based Cryptography

Lattirust is a Rust library for safe and efficient lattice-based cryptography, in particular for zero-knowledge proofs. 

Lattirust aims to be like [arkworks](https://github.com/arkworks-rs) for lattice-based constructions, and like [lattigo](https://github.com/tune-insight/lattigo) for zero-knowledge and succinct proofs. 

This repository is primarily maintained by [Christian Knabenhans](https://cknabs.github.io); see [CONTRIBUTORS.md](/CONTRIBUTORS.md) for a list of contributors. 
If you encounter any correctness, security, or performance issues, please [open an issue](/issues/new/choose); pull requests are also always welcome!

## Why Lattices?
Lattices allow for a wide range of quantum-secure cryptography. Due to the additional structure present in lattice assumptions compared to other post-quantum assumptions (e.g., collision-resistant hash functions), lattice-based cryptographic protocols are often more efficient. 

## Why Rust?
- Safety: Rust enforces memory safety, but also allows us to enforce more fine-grained security properties through Rust’s strong and expressive type system. 
- Performance: Rust allows for close-to-the-metal performance; in particular it is not garbage-collected (in constrast to other languages like Go), which allows for reliable, predictable, and reproducible benchmarks, which is important for academic and more exploratory contexts. 
- ZK ecosystem: the succinct and zero-knowledge proofs community (from both academia and industry) has written many artifacts and tools in Rust; we aim to interface with many of those tools, e.g., by being compatible with various frontends, or by using higher-order safety crates like [nimue](https://github.com/arkworks/nimue) for secure Fiat-Shamir transformations.
- Formal verification: Rust is also the language of choice for modern formal verification; medium-term, we hope to formally verify (parts of) lattirust, e.g. by using [hax](https://github.com/hacspec/hax). 

## Structure
At the moment, this repo holds the following Rust crates: 

### Core
- [lattirust_arithmetic](lattirust_arithmetic/): 
Implementations of power-of-two cyclotomic rings $\mathbb{Z}_q[X]/(X^{2^k}+1)$, number-theoretic transforms, matrices and vectors, various norms and challenge sets. 
Also contains wrappers and helper traits around the excellent [nimue](https://github.com/arkworks-rs/nimue) library (which enforces secure instantiations of the Fiat-Shamir transformation), making it easier to work with in lattice protocols.
- [lattice-estimator](lattice-estimator/):
Wrappers around the [lattice-estimator](https://github.com/malb/lattice-estimator) (for SIS) and [pq-crystals/security-estimates](https://github.com/pq-crystals/security-estimates) for MSIS).

### Protocols
- [labrador](labrador/):
An implementation of the [LaBRADOR protocol](https://eprint.iacr.org/2022/1341) (concretely small proof sizes, linear verifier). 
- [lova](lova/):
An implementation of the Lova lattice folding scheme.

## Security
Lattirust is provided for research and prototyping purposes, and has not been audited nor is fit for real-world deployment. Always consult your trusted cryptographer before using lattirust in security-critical protocols. 

## Roadmap
- Fall 2024:
  - Open-source release of lattirust v1 with labrador and lova
- End-of-year 2024: 
  - Proof-of-encryption / Proof-of-decryption for LWE/RLWE encryption (e.g., for FHE)
  - GPU / Icicle backend
  - in-tree Rust estimator for SIS/MSIS/vSIS/k-R-ISIS/BASIS/Power-BASIS
  - Greyhound proof system implemented
  - LatticeFold folding scheme implemented
  - hax specification of core arithmetic available (but not verified)
- Spring 2025:
  - hax specification for core arithmetic verified
  - hax specification for LWE/RLWE encryption/decryption verified
  - hax specification for LWE/RLWE key generation available
- Summer 2025:
  - hax specification for LWE/RLWE key generation verified

## GPU Acceleration with ICICLE CUDA Backend

- The current implementation uses a specific ICICLE revision specified in the cargo.lock file of `lattirust-arithmetic`. To switch to the latest ICICLE version, update the lines 31-32 in the cargo.toml file to the following:
```
icicle-runtime = { git = "https://github.com/ingonyama-zk/icicle.git", branch = "main" }
icicle-core = { git = "https://github.com/ingonyama-zk/icicle.git", branch = "main" }
icicle-babybear = { git = "https://github.com/ingonyama-zk/icicle.git", branch = "main" }
```

- The package requires the environment variable $ICICLE_BACKEND_INSTALL_DIR to be set to the path of the ICICLE backend installation directory.

- To enable the CUDA backend and GPU accelerated functions to use the GPU implementation, you should specify the `--features GPU` flag when building the package.

- The commands to run the tests for NTT and Vector Operations with the GPU backend are as follows:
```
$ cargo test --features GPU test_ntt_intt
$ cargo test --features GPU test_inner_products_icicle
```

- (Linux environment) To resolve the SAGE dependency, we use a conda environment. In case of issues with the compilers for ICICLE, you can override the default compilers of the environment by setting the following environment variables. Note that in more recent revisions of the ICICLE library, default compilers are set to clang, so this is not applicable in that case.
```
$ export CC=/usr/bin/gcc && export CXX=/usr/bin/g++
```

