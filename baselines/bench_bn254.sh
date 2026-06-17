#!/usr/bin/env bash
# BN254 external CPU baselines (pairing + G1 EC ops).

set -euo pipefail
cd "$(dirname "$0")"
source ./bench_lib.sh
bench_init bn254

run_timed "arkworks bn254 pairing" \
  bash -c 'cd algebra/curves/bn254 && cargo bench --bench bn254_pairing_only'

run_timed "arkworks bn254 g1 ec" \
  bash -c 'cd algebra/curves/bn254 && cargo bench --bench bn254_g1_ec_only'

run_timed "zksync bn254 pairing" \
  bash -c 'cd zksync-crypto/crates/pairing && cargo bench --bench pairing_benches bn256 -- "pairing"'

run_timed "zksync bn254 g1 ec" \
  bash -c 'cd zksync-crypto/crates/pairing && cargo bench --bench pairing_benches -- bn256::ec::g1::bench_g1_add_assign bn256::ec::g1::bench_g1_double'

run_timed "gnark bn254 pairing" \
  bash -c 'cd gnark-crypto/ecc/bn254 && go test -bench="BenchmarkPairing|BenchmarkMillerLoop|BenchmarkFinalExponentiation" -benchmem .'

run_timed "gnark bn254 g1 ec" \
  bash -c 'cd gnark-crypto/ecc/bn254 && go test -bench="BenchmarkG1JacAdd|BenchmarkG1JacAddMixed|BenchmarkG1JacDouble" -benchmem .'

bench_finish bn254
