#!/usr/bin/env bash
# BLS12-381 external CPU baselines (pairing + ed25519 field/EC comparison).
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

# curve25519-dalek / ed25519-dalek with AVX512 IFMA backend (nightly + unstable_avx512).
# Requires Intel CPU with avx512ifma and avx512vl at runtime (e.g. Ice Lake, Sapphire Rapids).
# Do NOT add global -C target-feature=+avx512ifma,+avx512vl — it breaks build scripts on non-IFMA hosts.
# Backend is selected at runtime via CPUID; -C target-cpu=native is safe on any x86_64 host.
export RUSTFLAGS='--cfg curve25519_dalek_backend="unstable_avx512" -C target-cpu=native'

run_timed "dalek ec_comparison serial" \
  bash -c 'cd curve25519-dalek && cargo +nightly bench -p curve25519-dalek \
    --features "alloc,rand_core,bench_internals" \
    --bench ec_comparison'

run_timed "dalek ec_comparison_ifma" \
  bash -c 'cd curve25519-dalek && cargo +nightly bench -p curve25519-dalek \
    --features "alloc,rand_core,bench_internals" \
    --bench ec_comparison_ifma'

run_timed "dalek full benchmarks" \
  bash -c 'cd curve25519-dalek && cargo +nightly bench -p curve25519-dalek \
    --features "alloc,rand_core,digest,precomputed-tables" \
    --bench dalek_benchmarks'

run_timed "ed25519 benchmarks" \
  bash -c 'cd curve25519-dalek && cargo +nightly bench -p ed25519-dalek \
    --features batch \
    --bench ed25519_benchmarks'

unset RUSTFLAGS

bench_finish bls12
