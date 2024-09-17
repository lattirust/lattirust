# Lattirust - A Rust Library for Lattice-based Cryptography

This repo is a (proof-of-concept, for now) Rust library for lattice-based cryptography, in particular for zero-knowledge proofs. 

Lattirust aims to be like [arkworks](https://github.com/arkworks-rs), but for lattice-based constructions, and like [lattigo](https://github.com/tune-insight/lattigo), but for zero-knowledge proofs. 

## Why Lattices?

Lattices allow for a wide range of quantum-secure cryptography. Due to the additional structure present in lattice assumptions compared to other post-quantum assumptions (e.g., collision-resistant hash functions), lattice-based cryptographic protocols are often more efficient. 

## Why Rust?
- Safety: Rust enforces memory safety, but also allows us to enforce more fine-grained security properties through Rust’s strong and expressive type system. 
- Performance: Rust allows for close-to-the-metal performance; in particular it is not garbage-collected (in contrast to other languages like Go), which allows for reliable, predictable, and reproducible benchmarks, which is important for academic and more exploratory contexts. 
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
