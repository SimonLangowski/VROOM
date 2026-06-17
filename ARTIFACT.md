# Artifact guide

Step-by-step instructions for reviewers reproducing this artifact. The CPU path is the primary deliverable; GPU benchmarks are optional (see `scripts/README.md`).

## 1. Supported environment

| Requirement | Details |
|-------------|---------|
| OS | Linux x86_64 (**Amazon Linux 2023** for paper numbers; Ubuntu also works) |
| Cloud (paper) | AWS **c7i** family — e.g. `c7i.metal-24xl` (benchmarks), `c7i-flex.large` (development) |
| CPU | Intel **Xeon Platinum 8488C** (Sapphire Rapids) or equivalent with **AVX-512 IFMA** |
| Compiler | **clang 21.1.0** (LLVM 21 prebuilt); GCC / older clang work for correctness but are much slower on IFMA code |
| RAM | ≥ 8 GiB |
| Packages | C++ toolchain, GMP, **Google Benchmark built with clang 21** (see §2) |

Optional (GPU — **separate machine** from CPU artifact; see **[ARTIFACT_GPU.md](ARTIFACT_GPU.md)**):

- NVIDIA GPU (paper: RTX 3090), CUDA toolkit, CGBN — `scripts/setup_gpu_deps.sh`

### Verify your host

Paper CPU numbers assume an IFMA-capable **c7i** instance. On the reference host this looks like:

```
NAME="Amazon Linux"
VERSION="2023"
PRETTY_NAME="Amazon Linux 2023.7.20250331"
...
6.1.131-143.221.amzn2023.x86_64
Model name: Intel(R) Xeon(R) Platinum 8488C
Flags: ... avx512f avx512dq ... avx512ifma ... avx512vl ...
```

Quick check:

```bash
cat /etc/os-release && uname -r
lscpu | grep -E 'Model name|Flags' | head -5
grep -q avx512ifma /proc/cpuinfo && echo "IFMA ok" || echo "Use fallback build (see below)"
export PATH="$HOME/LLVM-21.1.0-Linux-X64/bin:$PATH" 2>/dev/null; $CXX --version 2>/dev/null || true
```

Without `avx512ifma` in CPU flags, use `Makefile.fallback` / `FALLBACK=1` (correctness only, not paper timings).

### Compiler and AVX-512 IFMA performance

The production build uses `-mavx512ifma` throughout the hot path (`cpu/vector/avx512ifma.hpp`, Montgomery reduction, pairing benches). **Newer clang releases generate substantially better code for these intrinsics** than older compilers or GCC on the same hardware. We observed a large gap between distro **GCC** / older **clang** and **clang 21.x** on Sapphire Rapids (8488C): not only absolute speed, but also how consistently the backend schedules `_mm256_madd52lo_epu64` / `_mm256_madd52hi_epu64` and related IFMA sequences.

For **correctness**, any recent C++ compiler with IFMA support is fine. For **performance numbers comparable to the paper**, use **clang 21.1.0** (or newer clang from the LLVM 21+ line) via `scripts/setup_toolchain.sh` — not the Amazon Linux default clang 15 / GCC.

## 2. Install dependencies

Do **not** rely on `google-benchmark-devel` with **clang 21** on the reference host — the distro package is often built with a different compiler/ABI. Use the automated toolchain script (pinned versions below) or the manual steps.

### Pinned toolchain (reference host)

| Component | Source | Pin |
|-----------|--------|-----|
| **clang / clang++** | [LLVM 21.1.0](https://github.com/llvm/llvm-project/releases/tag/llvmorg-21.1.0) Linux x86_64 prebuilt | llvm-project @ `3623fe661ae35c6c80ac221f14d85be76aa870f1` |
| **Google Benchmark** | [google/benchmark](https://github.com/google/benchmark) | `04ccbd86038796c319ea19987457e651a24f6b44` (`v1.9.4-89-g04ccbd8`) |

Install paths on the reference host: `~/LLVM-21.1.0-Linux-X64/`, `~/benchmark/build2/`.

### One command (recommended)

```bash
chmod +x scripts/setup_toolchain.sh
./scripts/setup_toolchain.sh
eval "$(./scripts/setup_toolchain.sh --print)"
```

The script installs OS packages (`dnf` on Amazon Linux 2023, `apt` on Ubuntu), fetches the LLVM tarball (with sha256 check), clones Google Benchmark at the pinned commit, builds `libbenchmark.a` with `clang++`, and prints `export` lines for your shell.

### Session environment

After setup, every `make` and `reproduce_cpu_bench.sh` invocation needs:

```bash
eval "$(./scripts/setup_toolchain.sh --print)"
# or add the printed lines to ~/.bashrc
```

Expected `clang++ --version` on the reference host:

```
clang version 21.1.0 (https://github.com/llvm/llvm-project 3623fe661ae35c6c80ac221f14d85be76aa870f1)
Target: x86_64-unknown-linux-gnu
```

### Manual install (equivalent)

**Amazon Linux 2023:**

```bash
sudo dnf install -y gcc gcc-c++ make gmp-devel cmake git curl

# LLVM 21.1.0 prebuilt @ 3623fe661ae35c6c80ac221f14d85be76aa870f1
curl -LO https://github.com/llvm/llvm-project/releases/download/llvmorg-21.1.0/LLVM-21.1.0-Linux-X64.tar.xz
tar xf LLVM-21.1.0-Linux-X64.tar.xz -C "$HOME"
export PATH="$HOME/LLVM-21.1.0-Linux-X64/bin:$PATH"
export CC=clang CXX=clang++

# google/benchmark @ 04ccbd86038796c319ea19987457e651a24f6b44
git clone https://github.com/google/benchmark.git "$HOME/benchmark"
git -C "$HOME/benchmark" checkout 04ccbd86038796c319ea19987457e651a24f6b44
cmake -B "$HOME/benchmark/build2" -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER=clang++ -DBENCHMARK_ENABLE_GTEST_TESTS=OFF
cmake --build "$HOME/benchmark/build2"
export BENCHMARK_LIBS="$HOME/benchmark/build2/src/libbenchmark.a -lpthread"
export BENCHMARK_INC="-I$HOME/benchmark/include"
```

**Ubuntu** (also supported for builds; paper timings still use clang 21 + source-built benchmark):

```bash
sudo apt install build-essential libgmp-dev cmake git curl
./scripts/setup_toolchain.sh
```

Check IFMA support:

```bash
grep -q avx512ifma /proc/cpuinfo && echo "IFMA ok" || echo "Use fallback build (see below)"
```

## 3. Build

With clang 21 and a source-built Google Benchmark (see §2):

```bash
eval "$(./scripts/setup_toolchain.sh --print)"

cd blst && make gbench && cd ..
cd src && make bench-bls12-381
```

Default (`-lbenchmark` from the system package) works on Ubuntu when the distro library matches your compiler.

Without AVX-512 IFMA (correctness only, slow compile):

```bash
cd src && make -f Makefile.fallback
```

## 4. Minimal working examples

```bash
make -C examples
./examples/01_parameter_setup
./examples/02_singular_modmul
./examples/03_sum_of_products
./examples/04_difference_of_products
./examples/05_product_of_sums
```

Use `make -C examples FALLBACK=1` on machines without IFMA.

## 5. Smoke test (recommended verification)

From the repo root:

```bash
./scripts/smoke_test.sh
```

Runs BLST build, core unit tests, examples, and one field-multiply regression.

## 6. Full test and benchmark matrix

From `src/` (AVX-512 IFMA unless using `Makefile.fallback`):

| Command | Purpose |
|---------|---------|
| `make test` | Elementwise reduce + inversion |
| `make test-bls12-381-field-mul` | BLS12-381 `a*b mod p` |
| `make test-bls12-381-pairing` | BLS12-381 pairing |
| `make test-pairing` | Pairing (BLST baseline + RNS path) |
| `make test-ec` | EC point arithmetic |
| `make test-ec-bn254` | BN254 EC |
| `make test-scalar-mult` | Scalar multiplication |
| `make bench-bls12-381` | Paper CPU benchmarks (BLS12-381, Matrix + MatrixNoK) |
| `make bench-bls12-381` | BLS12-381 component benchmarks |
| `make bench-ec-bn254` | BN254 EC benchmarks |

## 7. Reproduce paper CPU numbers

**Reference machine** (reported numbers):

| | |
|--|--|
| Cloud | AWS **c7i.metal-24xl** |
| OS | **Amazon Linux 2023** (`amzn2023` kernel 6.1.x) |
| CPU | Intel **Xeon Platinum 8488C** (`avx512ifma` in `/proc/cpuinfo`) |
| Compiler | **clang 21.1.0** — LLVM prebuilt at `~/LLVM-21.1.0-Linux-X64/bin` (`CXX=clang++`) |

Development / smaller instances: **c7i-flex.large** on the same AMI family is sufficient for iteration; use metal for final timings.

**One command** (build, benchmark JSON, parsed table, resource log):

```bash
chmod +x scripts/setup_toolchain.sh scripts/reproduce_cpu_bench.sh
eval "$(./scripts/setup_toolchain.sh --print)"
./scripts/reproduce_cpu_bench.sh
```

Outputs under `results/`:

| File | Contents |
|------|----------|
| `bench_bls12_381.json` | Raw Google Benchmark JSON (RNS) |
| `bench_blst_complete.json` | Raw Google Benchmark JSON (BLST baseline) |
| `bench_bls12_381_table.txt` | Parsed RNS ns (Matrix + MatrixNoK); amortized ns/mul for batches |
| `bench_bls12_381_vs_blst.txt` | Aligned pairs + speedup ratio (`blst_ns / rns_ns`; >1 means RNS is faster) |
| `bench_bls12_381_resources.txt` | Wall time and peak RSS per phase + totals |

**Resource measurement (reference host, `results/bench_bls12_381_resources.txt`):**

Measured on AWS **c7i** / Amazon Linux 2023 (2026-06-17). Approximate resources **per paper claim**:

| Claim | Command | Wall time | Peak RAM | Disk (outputs + bench binaries) |
|-------|---------|-----------|----------|----------------------------------|
| **CPU benchmark table** (main) | `./scripts/reproduce_cpu_bench.sh` | **~2 minutes** (126 s end-to-end) | **~650 MB** (compile peak) | **~8 MB** (`du -sh results/`) |
| Benchmark run only (already built) | run `bench_bls12_381` + `bench_blst_complete` + parse | **~2 minutes** (RNS 86 s + BLST 22 s) | **~7 MB** | **~8 MB** (`results/`) |
| Build `bench_bls12_381` only | `make -C src bench_bls12_381` | **~15 seconds** | **~650 MB** | — |

Phase breakdown (full reproduce): build blst 2 s; build RNS bench 15 s; RNS benchmark 86 s; BLST benchmark 22 s; parse under 1 s.

Raw numbers are in `bench_bls12_381_resources.txt` (`total_elapsed_s=125.7`, `peak_rss_kb=663668`). Disk: **`du -sh results/` → 7.8 MB** on the reference host (JSON, parsed tables, resource log; the script’s `artifact_disk_kb` field sums only selected binaries and understates this).

RNS↔BLST name mapping lives in `scripts/parse_bench_json.py` (`RNS_BLST_MAP`). Unmapped entries are listed at the bottom of the comparison file for manual pairing.

BLST upstream revision: **supranational/blst @ `8065152`** (modified vendored copy in `blst/`; see §8 commit table).

Without IFMA: `FALLBACK=1 ./scripts/reproduce_cpu_bench.sh` (correctness only, not paper timings).

Parsed tables include **Matrix** and **MatrixNoK** rows by default. Use `--exclude-nok` on `reproduce_cpu_bench.sh` or `parse_bench_json.py` for Matrix only.

```bash
./scripts/reproduce_cpu_bench.sh --exclude-nok
```

Manual equivalent:

```bash
cd src
CXX=clang++ make bench-bls12-381
./bench_bls12_381 --benchmark_out=../results/bench_bls12_381.json --benchmark_out_format=json
python3 scripts/parse_bench_json.py results/bench_bls12_381.json --blst ../results/bench_blst_complete.json
```

Compare table to the sample in the root `README.md` (§ Running). `BM_BatchModMul_N` amortized time = reported cpu_ns ÷ N.

## 8. External CPU baselines (optional)

Third-party libraries for **optional** comparison: arkworks, zkcrypto, zksync-crypto, gnark-crypto. Vendored under **`baselines/`**; full details in **`baselines/README.md`**.

The **primary** paper CPU baseline is **BLST** (§7 above). External baselines are not required for the functional badge.

### Pinned upstream commits (paper bibliography)

| Dependency | Repository | Commit |
|------------|------------|--------|
| BLST (primary CPU baseline) | [supranational/blst](https://github.com/supranational/blst) | `8065152` |
| arkworks | [arkworks-rs/algebra](https://github.com/arkworks-rs/algebra) | `598a5fba` |
| zkcrypto | [zkcrypto/bls12_381](https://github.com/zkcrypto/bls12_381) | `6bb9695` |
| zksync | [matter-labs/zksync-crypto](https://github.com/matter-labs/zksync-crypto) | `e770ffd` |
| gnark | [Consensys/gnark-crypto](https://github.com/Consensys/gnark-crypto) | `d1dece6` |

`blst/` is a **modified vendored fork** based on `8065152` (benchmark wrappers, `extern "C"`, etc.) — use the in-tree copy for reproduction. Optional baselines are refreshed via `baselines/clone.sh`, which checks out the four external repos at the commits above and applies `baselines/overlays/`.

**Dependencies:** Rust/cargo, Go.

```bash
chmod +x baselines/reproduce_baselines.sh
./baselines/reproduce_baselines.sh          # BLS12-381 + BN254
./baselines/reproduce_baselines.sh bls12    # BLS12-381 pairing libs only
```

Outputs under `results/`:

| File | Contents |
|------|----------|
| `baselines_bls12_resources.txt` | Per-phase wall time + peak RSS + totals |
| `baselines_bls12_bench.log` | Raw Criterion / Go bench output |
| `baselines_bn254_resources.txt` | BN254 phases + totals |
| `baselines_bn254_bench.log` | Raw BN254 bench output |

**Resource measurement:** same badge pattern as §7 — one totals block per optional claim (first full run is dominated by Rust/Go compile; expect tens of minutes and several GiB peak RSS). There is no JSON parser aligning these to VROOM benchmark names; compare manually from the bench logs.

## 9. Parameter generation

### Precomputed (default)

Ready-to-use IntRNS4 exports for the paper curves are already wired in `src/`:

- **BLS12-381:** `bls12_381_ring_params.hpp`, `bls12_381_intrns4_system.hpp`, `bls12_381_perm*.hpp`
- **BN254:** `bn254_ring_params.hpp`, `bn254_intrns4_system.hpp`, `bn254_perm*.hpp`

Generated constants also live in `scripts/bls12_381_intrns4_params.hpp` and `scripts/bn254_intrns4_params.hpp`.

### Regenerating from Python

`scripts/rns_secp256k1.py` computes RNS/CRNS matrices and **writes C++ header files** (`print_intrns4_params`, `print_bls12_381_intrns4_params`, `print_bn254_intrns4_exports`). Run from `scripts/`:

```bash
python rns_secp256k1.py
```

### Advanced parameter generation (QR “no k” / MatrixNoK)

The quadratic-residue **“no k”** branch (`ChangeBaseVariation::MatrixNoK`, r1 cyclic permutation without k-correction) is **advanced parameter generation**:

- Runtime is **exponential in the number of radix words** — practical only for **small moduli** (the shipped BLS12-381 and BN254 configs).
- Higher redundance requires retuning: the MatrixNoK path uses **50-bit** (not 52-bit) words for the radix part of integer↔RNS conversion (`CONVERT_TO_WORD_BITS` in `bounded_ring.hpp`).

Helper scripts:

- **`scripts/gen_qr.py`** — parallel brute-force search for QR decomposition parameters.
- **`scripts/rns_secp256k1.py`** — full pipeline; emits `.hpp` parameter files consumed by `src/*_intrns4_system.hpp`.

After regenerating headers, rebuild `src/` and re-run `scripts/smoke_test.sh`.

## 10. Repository layout

| Path | Contents |
|------|----------|
| `src/` | `BoundedRing`, pairing/EC stack, tests, benchmarks |
| `cpu/` | AVX vector CRNS, Montgomery reduction, precompute |
| `examples/` | Minimal working examples |
| `scripts/` | Parameter generation, optional GPU benchmarks |
| `test_data/` | Reference test vectors (`test_data/README.md`) |
| `blst/` | BLST fork based on supranational/blst @ `8065152` (primary CPU baseline, FP12 reference) |
| `baselines/` | Optional external CPU baselines (`baselines/README.md`) |
| `gpu/` | CUDA kernels (optional; **separate GPU host** — see **ARTIFACT_GPU.md**) |

Module map and paper terminology: **`src/README.md`**.

## 11. Known limitations

- Research-grade code; not audited for production cryptography.
- Integer fallback (`Makefile.fallback`) validates correctness only; not for performance measurement.
- GPU path requires a separate CUDA/CGBN setup on an NVIDIA host — **[ARTIFACT_GPU.md](ARTIFACT_GPU.md)** (`scripts/setup_gpu_deps.sh`, `latency.py`).
