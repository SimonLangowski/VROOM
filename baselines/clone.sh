#!/bin/bash
cd "$(dirname "$0")"

git clone https://github.com/arkworks-rs/algebra.git

git clone https://github.com/zkcrypto/bls12_381.git

git clone https://github.com/matter-labs/zksync-crypto.git

git clone https://github.com/Consensys/gnark-crypto.git

# git clone https://github.com/dalek-cryptography/curve25519-dalek.git

git clone https://github.com/SimonLangowski/ed25519-dalek.git

# Local overlays: copies of files we add/modify in upstream repos (tracked under overlays/)
cp overlays/algebra/curves/bn254/Cargo.toml \
   algebra/curves/bn254/Cargo.toml
cp overlays/algebra/curves/bn254/benches/bn254_pairing_only.rs \
   algebra/curves/bn254/benches/bn254_pairing_only.rs
cp overlays/algebra/curves/bn254/benches/bn254_g1_ec_only.rs \
   algebra/curves/bn254/benches/bn254_g1_ec_only.rs

cp overlays/zksync-crypto/crates/pairing/benches/bn256/ec.rs \
   zksync-crypto/crates/pairing/benches/bn256/ec.rs
