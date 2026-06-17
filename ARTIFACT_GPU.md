# GPU artifact guide

Step-by-step instructions for reproducing the **GPU latency** experiments from the paper. This path is **independent** of the CPU (AVX-512 IFMA) artifact — use a CUDA-capable Linux host with an NVIDIA GPU.

**CPU artifact:** [ARTIFACT.md](ARTIFACT.md) (Amazon Linux 2023, AWS c7i, clang 21).

## 1. Supported environment

| Requirement | Details |
|-------------|---------|
| OS | Linux x86_64 — **Ubuntu 20.04 LTS** on the paper GPU host; other recent distros with NVIDIA drivers also work |
| GPU | NVIDIA GPU with **compute capability ≥ 7.0** (paper: **GeForce RTX 3090**, sm_86) |
| Driver | NVIDIA driver with CUDA 12.x support (paper host: **570.x**) |
| CUDA toolkit | **nvcc** on `PATH` (paper host: CUDA **12.2** toolkit; driver reports CUDA 12.8) |
| Python | **3.10+** with `numpy`, `matplotlib`, `cuda-python` |
| RAM | ≥ 8 GiB |
| Disk | ~500 MB for CGBN clone + benchmark CSV/plots |

**Not required:** AVX-512 IFMA, clang 21, Google Benchmark, or GMP on the GPU host (GMP is only needed if you rebuild CGBN samples).

### Verify your host

Paper GPU numbers were collected on a workstation with **RTX 3090**. On the reference host:

```bash
cat /etc/os-release && uname -r
nvidia-smi --query-gpu=name,driver_version,compute_cap --format=csv,noheader
nvcc --version
python3 -c "from cuda import cuda; import numpy, matplotlib; print('deps ok')"
```

Expected GPU line (from benchmark CSV metadata):

```
NVIDIA GeForce RTX 3090, 570.x, 8.6
```

Select a non-default GPU with `GPU_ID=1` (see §5).

## 2. Install dependencies

One-shot setup (clones [NVlabs/CGBN](https://github.com/NVlabs/CGBN) and installs Python packages):

```bash
chmod +x scripts/setup_gpu_deps.sh
./scripts/setup_gpu_deps.sh
```

Manual equivalent:

```bash
# NVIDIA driver + CUDA toolkit (distro packages or conda)
nvidia-smi
nvcc --version

git clone --depth 1 https://github.com/NVlabs/CGBN.git gpu/algorithms/baselines/CGBN

python3 -m pip install -r scripts/requirements-gpu.txt
```

**Conda alternative** (CUDA toolkit + cuda-python):

```bash
conda install -c conda-forge cuda-python
conda install -c "nvidia/label/cuda-12.2.1" cuda-toolkit
```

Set `CUDA_HOME` if `nvcc` is not on `PATH` (e.g. `/usr/local/cuda`).

## 3. Repository layout (GPU)

| Path | Contents |
|------|----------|
| `gpu/` | CUDA kernels — see **`gpu/README.md`** for module map |
| `scripts/` | Python drivers: `latency.py`, `gpu.py`, `bench_cgbn.py`, `basic_benchmarker.py` |
| `scripts/data/` | Sample benchmark CSVs from the paper host; new runs write `data/<timestamp>.csv` |
| `scripts/create_graph.py` | Parse CSV and plot latency curves |
| `results/` | Outputs from `reproduce_gpu_latency.sh` (`gpu_latency.csv`, `gpu_latency.png`, resources log) |

Parameter precompute for GPU moduli reuses `scripts/python_reduction.py` and `scripts/precompute_matrix.py` (invoked from `scripts/gpu.py`).

## 4. Minimal working example

After §2:

```bash
./scripts/smoke_test_gpu.sh
```

Runs one **CGBN** point (64-bit modulus) and one **MontRNSReduceVec** point (60-bit, stride 1) with correctness checks. Expect ~5 seconds and &lt; 200 MB RAM.

## 5. Configuration

| Variable / flag | Purpose |
|-----------------|--------|
| `GPU_ID` | CUDA device index (default `0`). Passed through `basic_benchmarker.py --gpu`. |
| `CUDA_HOME` | CUDA toolkit root for `nvcc` includes (default `/usr/local/cuda`). |
| `CGBN_REPO`, `CGBN_REF` | Override CGBN git URL/branch in `setup_gpu_deps.sh`. |
| `SKIP_SMOKE=1` | Skip smoke test in `reproduce_gpu_latency.sh`. |
| `RESULTS_DIR` | Output directory (default `results/`). |

No IP addresses or host-specific paths are hardcoded. Benchmarks use random moduli from `--modbits` unless `--modulus` is passed on the CLI.

## 6. Usage

### Full paper latency sweep

From the repo root:

```bash
chmod +x scripts/reproduce_gpu_latency.sh
./scripts/reproduce_gpu_latency.sh
```

This runs `python3 latency.py` inside `scripts/`, which:

1. Sweeps modulus sizes **64–4096 bits** (step 32).
2. For each size, benchmarks **CGBN** (`bench_cgbn.py`) and **MontRNSReduceVec** (`gpu.py` / `run_rns_latency_benchmark`) at stride 1 (and stride 2 for RNS).
3. Writes CSV to `scripts/data/<timestamp>.csv`, copies to `results/gpu_latency.csv`.
4. Plots **`results/gpu_latency.png`** (ns per modular multiplication vs modulus size).

### Manual equivalent

```bash
cd scripts
export MPLBACKEND=Agg
python3 latency.py
```

### Additional paper figures (optional)

`scripts/gpugraph.py` builds supplementary plots (CGBN vs RNS acceleration, TPI / warp-shuffle sweeps) from named CSVs under `scripts/data/`. See comments at the top of `gpugraph.py` for which sample files map to which figure.

```bash
cd scripts
python3 gpugraph.py
```

## 7. Reproduce paper GPU claims

### Reference machine (reported numbers)

| | |
|--|--|
| OS | **Ubuntu 20.04.6 LTS** (Focal) |
| GPU | **NVIDIA GeForce RTX 3090** (82 SMs, **sm_86**) |
| Driver | **570.207** (CUDA 12.8 reported by `nvidia-smi`) |
| Toolkit | CUDA **12.2** / `nvcc` 12.x |
| Python | 3.12, `numpy` 2.2, `matplotlib` 3.10, `cuda-python` 12.8 |

This is a **different machine** from the CPU artifact (AWS c7i / Amazon Linux 2023 / Xeon 8488C).

### Claim → script mapping

| Paper claim | Script / command | Output |
|-------------|------------------|--------|
| **GPU latency**: CGBN vs **MontRNSReduceVec** (modular multiplication time vs modulus bits, 64–4096) | `./scripts/reproduce_gpu_latency.sh` or `cd scripts && python3 latency.py` | `results/gpu_latency.csv`, `results/gpu_latency.png` |
| **Correctness** of GPU kernels (per configuration) | automatic in `basic_benchmarker.run_benchmark` (`correct` column in CSV) | `correct=True` rows in CSV |
| **Supplementary GPU figures** (baseline %, TPI / warp-shuffle) | `cd scripts && python3 gpugraph.py` | `*.png` in `scripts/` (see `gpugraph.py` for names) |

The latency CSV documents every configuration: columns `args_modbits`, `args_variation_name`, `args_stride`, `args_tpi`, `time`, `correct`, and `device` (GPU model string from CUDA).

### Resource measurement

Approximate resources **per claim** on the reference RTX 3090 host:

| Claim | Command | Wall time | Peak RAM | Disk |
|-------|---------|-----------|----------|------|
| **GPU latency sweep** | `./scripts/reproduce_gpu_latency.sh` | **~8 min** (492 s on RTX 3090; many `nvcc` kernel builds) | **~1.8 GB** | **~250 KB** CSV + PNG in `results/` |
| **Minimal GPU check** | `./scripts/smoke_test_gpu.sh` | **~5 s** | **~200 MB** | negligible |
| **Setup (first time)** | `./scripts/setup_gpu_deps.sh` | **~1 min** | — | **~50 MB** (CGBN clone) |

After a full run, exact numbers are in `results/gpu_latency_resources.txt` (`total_elapsed_s=491.68`, `peak_rss_kb≈1.8e6`, `sweep_rows=262` on the reference RTX 3090 host).

Sample paper-host CSV (shipped for comparison): `scripts/data/1748900378.1196628.csv`.

### Compare to reference output

```bash
cd scripts
python3 -c "
from create_graph import read_latency_csv
outx, outy = read_latency_csv('data/1748900378.1196628.csv')
print('reference series:', sorted(outx.keys()))
"
# Replot reference:
python3 -c "from create_graph import gpu_latency_plot; gpu_latency_plot('data/1748900378.1196628.csv')"
```

Timings vary by GPU model and driver; shapes (RNS faster than CGBN at midsize moduli) should match.

## 8. Smoke test (recommended verification)

```bash
./scripts/smoke_test_gpu.sh
```

## 9. Known limitations

- Each benchmark point **JIT-compiles** a specialized CUDA kernel via `nvcc`; the full sweep is slow but requires no separate build step.
- `nvcc` may warn that `-std=c++20` is ignored with the host GCC — benchmarks still run.
- CGBN is cloned at setup time (not vendored); pin with `CGBN_REF` if reproducibility across CGBN versions matters.
- Supplementary `gpugraph.py` figures depend on specific historical CSV filenames; regenerate those sweeps or edit paths in `gpugraph.py` if files are missing.

## 10. Related documentation

- **`gpu/README.md`** — CUDA module organization and paper terminology
- **`scripts/README.md`** — script index (CPU parameter generation + GPU)
- **`ARTIFACT.md`** — CPU / AVX artifact (separate host)
