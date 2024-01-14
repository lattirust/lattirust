# TODO
- [ ] R1CS prover/verifier
    - [ ] Check whether arkworks' normal form decomposition can be reused
- [ ] Have Matrix/Vector be structs instead of newtypes for better modularity
- [ ] Add some type-checking to nimue serialization; at the moment it's way too easy to accidentally deserialize to the wrong struct
- [ ] Investigate failing tests for operator norm
- [ ] Add parallelism to prover/verifier
    - [ ] ideally natively to matrix/vector as well