#!/bin/bash
cd "$(dirname "$0")"

# Run arkworks-rs/algebra BLS12-381 pairing benchmarks only
(cd algebra/curves/bls12_381 && cargo bench --bench bls12_381_pairing_only)

# Run zkcrypto/bls12_381 pairing benchmarks only (from groups.rs)
(cd bls12_381 && cargo bench --bench groups -- "pairing")

# Run zksync-crypto pairing BLS12-381 pairing benchmarks only
(cd zksync-crypto/crates/pairing && cargo bench --bench pairing_benches bls12_381 -- "pairing")

# Run gnark-crypto BLS12-381 pairing benchmarks (Go)
(cd gnark-crypto/ecc/bls12-381 && go test -bench="BenchmarkPairing|BenchmarkMillerLoop|BenchmarkFinalExponentiation" -benchmem .)

# curve25519-dalek / ed25519-dalek with AVX512 IFMA backend (nightly + unstable_avx512).
# Requires Intel CPU with avx512ifma and avx512vl at runtime (e.g. Ice Lake, Sapphire Rapids).
# Do NOT add global -C target-feature=+avx512ifma,+avx512vl — it breaks build scripts on non-IFMA hosts.
# Backend is selected at runtime via CPUID; -C target-cpu=native is safe on any x86_64 host.
export RUSTFLAGS='--cfg curve25519_dalek_backend="unstable_avx512" -C target-cpu=native'

# ModMul / G1_Add / G1_Double aligned with v4/bench_ec_secp256k1.cpp (serial field/point ops)
(cd curve25519-dalek && cargo +nightly bench -p curve25519-dalek \
    --features 'alloc,rand_core,bench_internals' \
    --bench ec_comparison)

# G1_Add / G1_MixedAdd / G1_Double on the AVX-512 IFMA vector backend (SIGILL if CPU lacks IFMA)
(cd curve25519-dalek && cargo +nightly bench -p curve25519-dalek \
    --features 'alloc,rand_core,bench_internals' \
    --bench ec_comparison_ifma)

# Full upstream dalek benchmarks (scalar mul, MSM, compression, etc.)
(cd curve25519-dalek && cargo +nightly bench -p curve25519-dalek \
    --features 'alloc,rand_core,digest,precomputed-tables' \
    --bench dalek_benchmarks)

# Ed25519 sign / verify / batch-verify (uses curve25519-dalek backend above)
(cd curve25519-dalek && cargo +nightly bench -p ed25519-dalek \
    --features batch \
    --bench ed25519_benchmarks)

unset RUSTFLAGS
