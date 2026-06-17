#!/usr/bin/env bash
# BLS12-381 external CPU baselines (pairing libraries).
# Vendored trees under baselines/; use clone.sh only for a fresh upstream checkout.

set -euo pipefail
cd "$(dirname "$0")"
source ./bench_lib.sh
bench_init bls12

run_timed "arkworks bls12_381 pairing" \
  bash -c 'cd algebra/curves/bls12_381 && cargo bench --bench bls12_381_pairing_only'

run_timed "zkcrypto bls12_381 pairing" \
  bash -c 'cd bls12_381 && cargo bench --bench groups -- "pairing"'

run_timed "zksync bls12_381 pairing" \
  bash -c 'cd zksync-crypto/crates/pairing && cargo bench --bench pairing_benches bls12_381 -- "pairing"'

run_timed "gnark bls12_381 pairing" \
  bash -c 'cd gnark-crypto/ecc/bls12-381 && go test -bench="BenchmarkPairing|BenchmarkMillerLoop|BenchmarkFinalExponentiation" -benchmem .'

bench_finish bls12
