# GPU CUDA kernels

CUDA implementation of **RNS Montgomery multiplication** and baselines used in the paper’s GPU latency experiments. Terminology matches the paper and `scripts/gpu.py` (`MontRNSReduceVec`, `Mont32`, TPI, stride).

CPU (AVX-512 IFMA) code lives under `cpu/` and `src/` — see **[ARTIFACT.md](../ARTIFACT.md)**. GPU evaluation uses a **separate machine**; see **[ARTIFACT_GPU.md](../ARTIFACT_GPU.md)**.

## Layout

| Path | Role |
|------|------|
| `benching/` | Benchmark harness: `loadingvec.cuh` loads operands and runs `bench_modmult2`; `include.cuh` pulls in algorithm headers |
| `algorithms/baselines/` | **CGBN** baseline (`cgbn_mul.cuh`) and redundant Montgomery reference (`montgomery.cuh`) |
| `algorithms/reductions/` | Elementwise **Montgomery reduction** (`Mont32`, `MontRedundant30`) on 32-bit RNS limbs |
| `algorithms/matrix/` | RNS matrix multiply fragments, epilogue, and `mod_multiply_rns` in `parts.cuh` |
| `algorithms/vector/` | **MontRNSReduceVec** — warp-level accumulation and change-base (`vec_acc.cuh`) |
| `crypto/` | `bench_modmult2` microbenchmark (repeated `mulRed`) |
| `debug.cuh` | Optional device-side tracing when `DEBUG` is set |

## Paper ↔ code mapping

| Paper term | Code |
|------------|------|
| CGBN baseline | `algorithms/baselines/cgbn_mul.cuh` + NVlabs CGBN (`setup_gpu_deps.sh`) |
| MontRNSReduceVec | `variation == 4` in `scripts/gpu.py`; kernel template in `benching/loadingvec.cuh` |
| Mont32 reduction | `algorithms/reductions/montgomery.cuh` (`montgomery_reduce32`) |
| TPI (threads per instance) | `PYREPLACE_TPI` / `--btpi` — each big integer spread across a power-of-2 thread group |
| stride | `PYREPLACE_STRIDE` / `--stride` — over-parallelization factor (1 or 2) |

Kernels are **not built with a standalone Makefile**. `scripts/basic_benchmarker.py` substitutes `PYREPLACE_*` placeholders into a `.cuh` file, writes `sourceCode.cu`, invokes `nvcc --cubin`, and loads the module at runtime.

## Non-obvious design notes

- **Include paths in `.cuh` files** are written for generated `scripts/sourceCode.cu` (repo root relative), not for in-tree compilation of the `.cuh` alone.
- **MontRNSReduceVec** packs precomputed change-base matrices (`pre_ptr1`, `pre_ptr2`) and moduli per thread block; Python (`scripts/gpu.py`) computes layouts via `precompute_matrix.py` / `python_reduction.py`.
- **Correctness**: each benchmark run checks GPU output against Python big-int reference (`correct_func_*` in `scripts/gpu.py`, `scripts/bench_cgbn.py`).
- **CGBN** is an external dependency cloned by `scripts/setup_gpu_deps.sh` into `algorithms/baselines/CGBN/` (not vendored in git).

## Quick start

```bash
./scripts/setup_gpu_deps.sh
./scripts/smoke_test_gpu.sh          # one CGBN + one RNS point
./scripts/reproduce_gpu_latency.sh   # full paper latency sweep → results/
```

See **ARTIFACT_GPU.md** for environment, claim-to-script mapping, and resource estimates.
