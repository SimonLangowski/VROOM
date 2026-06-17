# scripts

Python tooling for **RNS parameter generation** and optional **GPU benchmarks**. The CPU artifact does not require Python unless you regenerate parameters or run GPU experiments.

## CPU toolchain (clang 21 + Google Benchmark)

Paper CPU timings use **clang 21.1.0** and a **clang-built** `libbenchmark.a`. Automated setup (pinned commits, same as the reference host):

```bash
chmod +x scripts/setup_toolchain.sh scripts/smoke_test.sh scripts/reproduce_cpu_bench.sh
./scripts/setup_toolchain.sh
eval "$(./scripts/setup_toolchain.sh --print)"
```

| Component | Pin |
|-----------|-----|
| LLVM 21.1.0 prebuilt | llvm-project @ `3623fe661ae35c6c80ac221f14d85be76aa870f1` |
| google/benchmark | `04ccbd86038796c319ea19987457e651a24f6b44` |

See **ARTIFACT.md** §2 for manual steps and verification.

## Parameter generation (CPU artifact)

| Script | Role |
|--------|------|
| **`rns_secp256k1.py`** | Main pipeline: IntRNS2/IntRNS4 matrices, Montgomery factors, convert_to/from. **Prints C++ header files** (e.g. `bls12_381_intrns4_params.hpp`, `bn254_intrns4_params.hpp`) and can refresh `test_data/test_values_*.hpp`. |
| **`gen_qr.py`** | Helper for **advanced parameter generation**: parallel brute-force search for quadratic-residue (“QR no k”) decomposition used by the MatrixNoK change-base path. |
| `rns_helpers.py`, `precompute_matrix.py`, `precompute_montgomery.py` | Shared linear-algebra and Montgomery precompute routines. |

### Precomputed parameters (default)

You do **not** need to run Python to build the shipped artifact. Precomputed IntRNS4 exports for **BLS12-381** and **BN254** are already in `src/`:

- `bls12_381_ring_params.hpp`, `bls12_381_intrns4_system.hpp`, `bls12_381_perm*.hpp`
- `bn254_ring_params.hpp`, `bn254_intrns4_system.hpp`, `bn254_perm*.hpp`

Raw generated constants: `scripts/bls12_381_intrns4_params.hpp`, `scripts/bn254_intrns4_params.hpp`.

### Regenerating headers

```bash
cd scripts
python rns_secp256k1.py
cd ../src && make
```

### Advanced parameter generation (QR “no k” / MatrixNoK)

The QR **“no k”** method (`ChangeBaseVariation::MatrixNoK`) is labeled **advanced parameter generation** in this artifact:

- Search cost is **exponential in the number of radix words** — feasible only for **small moduli** (the included BLS12-381 and BN254 configs).
- Higher redundance requires different tuning; MatrixNoK uses **50-bit** (not 52-bit) radix words for the integer↔RNS conversion path.

Workflow: run `gen_qr.py` (parallel brute force) to find QR candidates, then `rns_secp256k1.py` to emit `.hpp` files and wire them through `src/*_intrns4_system.hpp`.

See **ARTIFACT.md** § Parameter generation for the full reviewer workflow.

### External CPU baselines (optional)

Vendored third-party trees under **`../baselines/`** (arkworks, gnark, zksync, zkcrypto). Not required for the primary CPU artifact.

```bash
./baselines/reproduce_baselines.sh
```

Writes `results/baselines_*_resources.txt` (per-phase time/RSS) and `results/baselines_*_bench.log`. See **`baselines/README.md`**.

---

## GPU benchmarks (optional, separate host)

GPU evaluation uses a **different machine** than the CPU (AVX) artifact. Full reviewer instructions: **[ARTIFACT_GPU.md](../ARTIFACT_GPU.md)**.

| Script | Role |
|--------|------|
| **`setup_gpu_deps.sh`** | Clone NVlabs CGBN + `pip install -r requirements-gpu.txt` |
| **`smoke_test_gpu.sh`** | Minimal working example (one CGBN + one RNS point) |
| **`reproduce_gpu_latency.sh`** | Full paper latency sweep → `results/gpu_latency.*` |
| **`latency.py`** | Sweep 64–4096 bit moduli: CGBN vs **MontRNSReduceVec** |
| **`bench_cgbn.py`** | CGBN baseline hooks (`gen`, `result`, `correct`) |
| **`gpu.py`** | RNS Montgomery GPU path (same hook pattern) |
| **`basic_benchmarker.py`** | JIT `nvcc` compile, run kernels, check correctness, time |
| **`create_graph.py`** | Parse CSV and plot (`gpu_latency_plot`) |
| **`gpugraph.py`** | Supplementary paper figures from named CSVs in `data/` |

### Quick start

```bash
./scripts/setup_gpu_deps.sh
./scripts/smoke_test_gpu.sh
./scripts/reproduce_gpu_latency.sh   # ~8 min on RTX 3090; runs: cd scripts && python3 latency.py
```

Run from `scripts/` when invoking Python directly (`python3 latency.py`). Requires `nvidia-smi`, `nvcc`, and `cuda-python`.

Sample paper-host CSV: `scripts/data/1748900378.1196628.csv`.
