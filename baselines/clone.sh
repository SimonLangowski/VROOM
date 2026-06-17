#!/usr/bin/env bash
# Clone optional external baselines at paper-pinned commits, then apply overlays/.
# Commits match the bibliography entries in ARTIFACT.md §8.
set -euo pipefail
cd "$(dirname "$0")"

clone_pinned() {
  local dir="$1"
  local url="$2"
  local commit="$3"
  echo "== $dir @ $commit =="
  rm -rf "$dir"
  git clone "$url" "$dir"
  git -C "$dir" checkout "$commit"
  echo "    checked out $(git -C "$dir" rev-parse --short HEAD)"
}

clone_pinned algebra \
  https://github.com/arkworks-rs/algebra.git \
  598a5fba

clone_pinned bls12_381 \
  https://github.com/zkcrypto/bls12_381.git \
  6bb9695

clone_pinned zksync-crypto \
  https://github.com/matter-labs/zksync-crypto.git \
  e770ffd

clone_pinned gnark-crypto \
  https://github.com/Consensys/gnark-crypto.git \
  d1dece6

# Local overlays: copies of files we add/modify in upstream repos (tracked under overlays/)
cp overlays/algebra/curves/bn254/Cargo.toml \
   algebra/curves/bn254/Cargo.toml
cp overlays/algebra/curves/bn254/benches/bn254_pairing_only.rs \
   algebra/curves/bn254/benches/bn254_pairing_only.rs
cp overlays/algebra/curves/bn254/benches/bn254_g1_ec_only.rs \
   algebra/curves/bn254/benches/bn254_g1_ec_only.rs

cp overlays/algebra/curves/bls12_381/Cargo.toml \
   algebra/curves/bls12_381/Cargo.toml
cp overlays/algebra/curves/bls12_381/benches/bls12_381_pairing_only.rs \
   algebra/curves/bls12_381/benches/bls12_381_pairing_only.rs

cp overlays/zksync-crypto/crates/pairing/benches/bn256/ec.rs \
   zksync-crypto/crates/pairing/benches/bn256/ec.rs

echo "Done. Overlays applied."
