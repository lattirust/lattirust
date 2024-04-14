# Lattirust - A Rust Library for Lattice-based Cryptography

This repo is a (proof-of-concept, for now) Rust library for lattice-based cryptography, in particular for zero-knowledge proofs. 

Lattirust aims to be like [arkworks](https://github.com/arkworks-rs), but for lattice-based constructions, and like [lattigo](https://github.com/tune-insight/lattigo), but for zero-knowledge proofs. 

> [!CAUTION]
> This repository contains a proof-of-concept implementation of a lattice-based folding scheme, with a corresponding paper to be submitted soon; please don't circulate this code until this paper is made public on ePrint.

## Why Lattices?

Lattices allow for a wide range of quantum-secure cryptography. Due to the additional structure present in lattice assumptions compared to other post-quantum assumptions (e.g., collision-resistant hash functions), lattice-based cryptographic protocols are often more efficient. 

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
