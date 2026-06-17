#!/bin/bash
cd "$(dirname "$0")"

# arkworks-rs/algebra BN254 pairing benchmarks
(cd algebra/curves/bn254 && cargo bench --bench bn254_pairing_only)

# arkworks-rs/algebra BN254 G1 EC ops (add / mixed add / double)
(cd algebra/curves/bn254 && cargo bench --bench bn254_g1_ec_only)

# zksync-crypto BN254 (bn256) pairing benchmarks
(cd zksync-crypto/crates/pairing && cargo bench --bench pairing_benches bn256 -- "pairing")

# zksync-crypto BN254 G1 EC ops (add / mixed add / double)
(cd zksync-crypto/crates/pairing && cargo bench --bench pairing_benches -- bn256::ec::g1::bench_g1_add_assign bn256::ec::g1::bench_g1_double)

# gnark-crypto BN254 pairing benchmarks
(cd gnark-crypto/ecc/bn254 && go test -bench="BenchmarkPairing|BenchmarkMillerLoop|BenchmarkFinalExponentiation" -benchmem .)

# gnark-crypto BN254 G1 EC ops (add / mixed add / double)
(cd gnark-crypto/ecc/bn254 && go test -bench="BenchmarkG1JacAdd|BenchmarkG1JacAddMixed|BenchmarkG1JacDouble" -benchmem .)
