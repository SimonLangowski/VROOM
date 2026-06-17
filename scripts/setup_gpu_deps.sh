#!/usr/bin/env bash
# Install GPU benchmark dependencies: CGBN (NVlabs) and Python packages.
#
# Usage (from repo root):
#   ./scripts/setup_gpu_deps.sh
#
# Environment overrides:
#   CGBN_REPO   — git URL (default: https://github.com/NVIDIA/CGBN.git)
#   CGBN_REF    — branch or tag (default: master)
#   SKIP_PIP=1  — skip pip install

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
CGBN_DIR="$ROOT/gpu/algorithms/baselines/CGBN"
CGBN_REPO="${CGBN_REPO:-https://github.com/NVlabs/CGBN.git}"
CGBN_REF="${CGBN_REF:-master}"

echo "== Verify NVIDIA GPU + CUDA =="
if ! command -v nvidia-smi >/dev/null 2>&1; then
  echo "nvidia-smi not found. Install NVIDIA drivers before GPU benchmarks." >&2
  exit 1
fi
nvidia-smi --query-gpu=name,driver_version,compute_cap --format=csv,noheader

if ! command -v nvcc >/dev/null 2>&1; then
  echo "nvcc not found. Install the CUDA toolkit (nvcc on PATH or set CUDA_HOME)." >&2
  exit 1
fi
nvcc --version | head -1

echo ""
echo "== CGBN (GPU baseline library) =="
if [[ -f "$CGBN_DIR/include/cgbn/cgbn.h" ]]; then
  echo "CGBN already present at $CGBN_DIR"
else
  echo "Cloning $CGBN_REPO (ref $CGBN_REF) -> $CGBN_DIR"
  mkdir -p "$(dirname "$CGBN_DIR")"
  git clone --depth 1 --branch "$CGBN_REF" "$CGBN_REPO" "$CGBN_DIR"
fi

if [[ "${SKIP_PIP:-0}" != "1" ]]; then
  echo ""
  echo "== Python packages =="
  python3 -m pip install --user -r "$ROOT/scripts/requirements-gpu.txt"
fi

echo ""
echo "GPU dependencies ready. Next:"
echo "  cd scripts && ./smoke_test_gpu.sh"
echo "  cd scripts && ./reproduce_gpu_latency.sh"
