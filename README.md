# RNS Montgomery multiplication with AVX512

Standalone C++ template/header code. Benchmarks use Google Benchmark. Inversion with BLST requires BLST with extern "C" modification.

**Artifact reviewers:** start with **[ARTIFACT.md](ARTIFACT.md)** (environment, smoke test, parameter generation).

## Repository layout

| Path | Role |
|------|------|
| [`src/`](src/) | RNS rings, pairing/EC code, tests, benchmarks |
| [`cpu/`](cpu/) | AVX vector CRNS, Montgomery reduction, precompute |
| [`examples/`](examples/) | Minimal working examples |
| [`scripts/`](scripts/) | Parameter generation (`rns_secp256k1.py`), optional GPU benchmarks |
| [`test_data/`](test_data/) | Reference test vectors |
| [`blst/`](blst/) | BLST fork — primary CPU baseline for paper numbers |
| [`baselines/`](baselines/) | Optional external CPU baselines (arkworks, gnark, zksync, dalek) |
| [`gpu/`](gpu/) | CUDA kernels (**optional**, not required for CPU artifact) |

## Supported environment

| | |
|--|--|
| OS | Linux x86_64 — **Amazon Linux 2023** on AWS for paper numbers; Ubuntu also works |
| Cloud | AWS **c7i** (Intel Xeon Platinum **8488C**, Sapphire Rapids) |
| CPU | AVX-512 **IFMA** required for production build (`grep avx512ifma /proc/cpuinfo`) |
| Compiler | **clang 21.1.0** (LLVM 21 prebuilt) on reference host — see [ARTIFACT.md](ARTIFACT.md) §2 |
| Packages | C++ toolchain, GMP, Google Benchmark — **build from source with clang 21 on AL2023** ([ARTIFACT.md](ARTIFACT.md) §2) |

Verify IFMA on a fresh instance:

```bash
cat /etc/os-release && uname -r
lscpu | grep -E 'Model name|Flags' | head -5
grep -q avx512ifma /proc/cpuinfo && echo "IFMA ok" || echo "Use fallback build"
```

GPU benchmarks are optional; see [`scripts/README.md`](scripts/README.md).

## Building & testing

Install dependencies — **Amazon Linux 2023** (paper host):

```bash
sudo dnf install -y gcc gcc-c++ make gmp-devel cmake git
# Google Benchmark: build from source with clang++ (see ARTIFACT.md §2)
export BENCHMARK_LIBS="$HOME/benchmark/build2/src/libbenchmark.a -lpthread"
export BENCHMARK_INC="-I$HOME/benchmark/include"
```

**Ubuntu:**

```bash
sudo apt install build-essential libgmp-dev libbenchmark-dev
```

Build:

```bash
cd blst
make
cd ../src
make
```

Requires **AVX-512 IFMA** (`-mavx512ifma`). For integer-simulated vector ops (debugging / machines without IFMA), use `make -f Makefile.fallback` in `src/` — compiles slowly and is not used for performance work.

Quick verification:

```bash
chmod +x scripts/smoke_test.sh
./scripts/smoke_test.sh
```

Use `FALLBACK=1 ./scripts/smoke_test.sh` without IFMA.

### Minimal working examples

```bash
make -C examples
./examples/01_parameter_setup
./examples/02_singular_modmul
./examples/03_sum_of_products
./examples/04_difference_of_products
./examples/05_product_of_sums
```

See `examples/README.md` for details. Use `make -C examples FALLBACK=1` without IFMA.

### Tests & benchmarks

| Target | Purpose |
|--------|---------|
| `make test` | Elementwise reduce + inversion |
| `make test-bls12-381-field-mul` | BLS12-381 field multiply |
| `make test-bls12-381-pairing` | BLS12-381 pairing |
| `make test-pairing` | Pairing (AVX512 IFMA) |
| `make test-ec` | EC point arithmetic |
| `make test-ec-bn254` | BN254 EC |
| `make test-scalar-mult` | Scalar multiplication |
| `make bench-bls12-381` | Paper CPU benchmarks (BLS12-381, Matrix + MatrixNoK) |
| `make bench-bls12-381` | BLS12-381 component benchmarks |
| `make bench-ec-bn254` | BN254 EC benchmarks |
| `make -f Makefile.fallback test-pairing` | Pairing test with integer fallback (correctness only) |

Module map and parameter files: [`src/README.md`](src/README.md). Parameter regeneration: [`ARTIFACT.md`](ARTIFACT.md) § Parameter generation.

### External CPU baselines (optional)

VROOM vs BLST is automated in `scripts/reproduce_cpu_bench.sh`. For arkworks, gnark, zksync, zkcrypto, and curve25519-dalek comparisons, see [`baselines/README.md`](baselines/README.md) and `./baselines/reproduce_baselines.sh`.

## Running

**Paper numbers:** AWS **c7i.metal-24xl**, **Amazon Linux 2023**, Intel **Xeon Platinum 8488C**, **clang 21.1.0** (`~/LLVM-21.1.0-Linux-X64/bin/clang++`).

**Development:** AWS **c7i-flex.large** on the same AMI family is enough for iteration.

GCC is significantly slower than clang on this codebase; clang **21.x** (we used **21.1.0**) is the reference for paper timings.

- **`./bench_bls12_381`** — paper CPU suite (`BM_*_Matrix` and `BM_*_MatrixNoK` variants). Parsed tables include both by default; pass `--exclude-nok` to `parse_bench_json.py` for Matrix only.

Sample output on `c7i.metal-24xl` with **clang 21.1.0** (Matrix rows; abbreviated):

```
Running ./bench_bls12_381
...
BM_ModMul_Matrix                         28.4 ns         28.4 ns     ... ModMul_Matrix
BM_BatchModMul_Matrix_1                  28.3 ns         28.3 ns     ... BatchModMul_Matrix_1
...
BM_MillerLoop_Matrix                  58985 ns        58985 ns        ... MillerLoop_Matrix
BM_Pairing_RNS_BLST_Inverter_Matrix  151537 ns       151534 ns        ... Pairing_RNS_BLST_Inverter_Matrix
```

- **BM_BatchModMul_N** — batch Montgomery multiplication, `a*b mod p`.
- **RNS_BLST_Inverter** — uses RNS for FP12 inversion and BLST for FP inversion.
- **BM_Inverse** — uses an addition chain of length 425 to invert in RNS without BLST dependency (compare to **BM_FP12_Inverse_BLST** and  **BM_FP12_Inverse_RNS_BLST**)

## GPU (optional)

CPU artifact evaluation does **not** require GPU setup. For paper GPU figures, see [`scripts/README.md`](scripts/README.md) (CUDA toolkit, CGBN, `latency.py`, `gpugraph.py`).

## License

MIT. The **blst/** subtree keeps its original Apache license and copyright (see `blst/LICENSE`).

## Security

This is research-grade code, demonstrating the speedups of combining AVX512-IFMA, RNS Montgomery, and the BLS12-381 curve.
It has not been vetted or audited for security.