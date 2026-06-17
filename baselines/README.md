# External CPU baselines

Third-party pairing and EC libraries used for **optional** comparison against VROOM. The **primary** in-tree baseline for paper CPU numbers is **BLST** (`blst/`, compared via `scripts/reproduce_cpu_bench.sh`). These trees are supplementary context (arkworks, gnark, zksync, zkcrypto, curve25519-dalek).

Vendored copies live under this directory so reviewers do not need network access. To refresh from upstream, run `clone.sh` (then re-apply `overlays/` for BN254 bench hooks).

## Layout

| Path | Upstream | What we benchmark |
|------|----------|-------------------|
| `algebra/` | [arkworks-rs/algebra](https://github.com/arkworks-rs/algebra) | BLS12-381 / BN254 pairing; BN254 G1 EC |
| `bls12_381/` | [zkcrypto/bls12_381](https://github.com/zkcrypto/bls12_381) | BLS12-381 pairing |
| `zksync-crypto/` | [matter-labs/zksync-crypto](https://github.com/matter-labs/zksync-crypto) | BLS12-381 / BN256 pairing + BN256 G1 EC |
| `gnark-crypto/` | [Consensys/gnark-crypto](https://github.com/Consensys/gnark-crypto) | BLS12-381 / BN254 pairing + G1 EC (Go) |
| `curve25519-dalek/` | [SimonLangowski/ed25519-dalek](https://github.com/SimonLangowski/ed25519-dalek) fork | ed25519 field mul + G1 ops (serial and AVX-512 IFMA vector) |
| `overlays/` | — | Patched bench files applied by `clone.sh` (BN254 arkworks + zksync G1 benches) |

## Dependencies

| Tool | Required for |
|------|----------------|
| **Rust / cargo** | arkworks, zkcrypto, zksync, curve25519-dalek |
| **rustup + nightly** | `curve25519-dalek` IFMA benches (`cargo +nightly`) |
| **Go** | gnark-crypto |
| **AVX-512 IFMA** (runtime) | `ec_comparison_ifma` (SIGILL without IFMA) |

Install example (Ubuntu):

```bash
sudo apt install golang-go
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
rustup toolchain install nightly
```

## Running

**One command** (records wall time + peak RSS per phase, full bench log):

```bash
chmod +x baselines/reproduce_baselines.sh baselines/bench.sh baselines/bench_bn254.sh
./baselines/reproduce_baselines.sh          # BLS12-381 + BN254
./baselines/reproduce_baselines.sh bls12    # BLS12-381 only
./baselines/reproduce_baselines.sh bn254    # BN254 only
```

Outputs under `results/`:

| File | Contents |
|------|----------|
| `baselines_bls12_resources.txt` | Per-phase elapsed / max RSS + totals |
| `baselines_bls12_bench.log` | Raw `cargo bench` / `go test -bench` stdout |
| `baselines_bn254_resources.txt` | Same for BN254 |
| `baselines_bn254_bench.log` | Raw BN254 bench output |

Or run phases manually:

```bash
cd baselines
./bench.sh          # BLS12-381 pairing libs + dalek/ed25519
./bench_bn254.sh    # BN254 pairing + G1 EC
```

**First run** compiles Rust and Go dependencies (long wall time, several GiB disk under `baselines/*/target`). Subsequent runs are faster.

## What is measured

- **Pairing stacks** (arkworks, zkcrypto, zksync, gnark): Miller loop, final exponentiation, full pairing — via each project's own Criterion / Go benchmarks. Output is human-readable ns/op in the bench log, not JSON aligned to VROOM names.
- **curve25519-dalek** (`ec_comparison`, `ec_comparison_ifma`): `ModMul`, `G1_Add`, `G1_MixedAdd`, `G1_Double` on **ed25519** (pseudomersenne field), aligned with legacy `bench_ec_secp256k1` naming — **not** BLS12-381 pairing. Useful as a second CPU vector baseline for field/EC micro-ops.

There is **no** automated RNS↔arkworks ratio table (unlike `parse_bench_json.py` for BLST). Compare numbers manually from the bench logs.

## Resource measurement

Same pattern as `scripts/reproduce_cpu_bench.sh`: `/usr/bin/time` peak RSS and wall clock per phase, plus script totals suitable for ARTIFACT.md. Baselines are **optional** for the functional badge; record totals only if you reproduce them for the paper.

Expect rough order of magnitude on a clean build (IFMA host, network for first `cargo` fetch if vendored trees are missing):

| Script | Typical first-run wall time | Peak RSS (compile) |
|--------|----------------------------|--------------------|
| `bench.sh` (BLS12) | tens of minutes | several GiB |
| `bench_bn254.sh` | ~10–20 min | several GiB |

Exact numbers depend on CPU, disk, and whether `target/` already exists.

## Fresh clone from upstream

```bash
cd baselines
./clone.sh    # clones algebra, bls12_381, zksync-crypto, gnark-crypto, ed25519-dalek fork
              # applies overlays/ for BN254 bench files
```

`curve25519-dalek` is vendored in-repo (not cloned by `clone.sh`).
