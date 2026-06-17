#!/usr/bin/env bash
# Reproduce GPU latency figure (paper): CGBN baseline vs MontRNSReduceVec sweep.
#
# Usage (from repo root):
#   ./scripts/reproduce_gpu_latency.sh
#   SKIP_SMOKE=1 ./scripts/reproduce_gpu_latency.sh
#   GPU_ID=1 ./scripts/reproduce_gpu_latency.sh
#
# Outputs under results/:
#   gpu_latency.csv          — timing CSV (all sweep points)
#   gpu_latency.png          — latency comparison plot
#   gpu_latency_resources.txt — wall time, peak RSS, disk

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SCRIPT_DIR="$ROOT/scripts"
RESULTS_DIR="${RESULTS_DIR:-$ROOT/results}"
CSV_OUT="$RESULTS_DIR/gpu_latency.csv"
PNG_OUT="$RESULTS_DIR/gpu_latency.png"
RESOURCES_OUT="$RESULTS_DIR/gpu_latency_resources.txt"
CGBN_HDR="$ROOT/gpu/algorithms/baselines/CGBN/include/cgbn/cgbn.h"

mkdir -p "$RESULTS_DIR"
export MPLBACKEND=Agg
export GPU_ID="${GPU_ID:-0}"

if [[ ! -f "$CGBN_HDR" ]]; then
  echo "CGBN not found — running setup_gpu_deps.sh"
  "$ROOT/scripts/setup_gpu_deps.sh"
fi

if [[ "${SKIP_SMOKE:-0}" != "1" ]]; then
  echo "== Phase: GPU smoke test =="
  "$SCRIPT_DIR/smoke_test_gpu.sh"
fi

run_timed() {
  local label="$1"
  shift
  local tf
  tf="$(mktemp)"
  local t0 t1
  t0=$(date +%s.%N)
  /usr/bin/time -f '%M' -o "$tf" "$@"
  t1=$(date +%s.%N)
  PHASE_RSS_KB=$(cat "$tf")
  rm -f "$tf"
  PHASE_ELAPSED=$(python3 -c "print(round($t1 - $t0, 2))")
  echo "  [$label] elapsed=${PHASE_ELAPSED}s  max_rss=${PHASE_RSS_KB}KB"
}

SCRIPT_START=$(date +%s.%N)
TOTAL_RSS_KB=0

{
  echo "GPU latency reproduction (latency.py)"
  echo "date: $(date -Is)"
  echo "host: $(uname -n)"
  echo "gpu_id: $GPU_ID"
  nvidia-smi --query-gpu=name,driver_version,compute_cap --format=csv,noheader 2>/dev/null || true
  echo ""
} | tee "$RESOURCES_OUT"

echo "== Phase: full latency sweep (CGBN + MontRNSReduceVec) =="
echo "This compiles many CUDA kernels; expect ~8–15 minutes on an RTX 3090-class GPU."
cd "$SCRIPT_DIR"
run_timed "latency.py" python3 latency.py
echo "latency_sweep_elapsed_s=$PHASE_ELAPSED" >> "$RESOURCES_OUT"
echo "latency_sweep_max_rss_kb=$PHASE_RSS_KB" >> "$RESOURCES_OUT"
TOTAL_RSS_KB=$(( TOTAL_RSS_KB > PHASE_RSS_KB ? TOTAL_RSS_KB : PHASE_RSS_KB ))

# latency.py writes data/<timestamp>.csv and latency.png in scripts/
LATEST_CSV="$(ls -t "$SCRIPT_DIR/data/"*.csv 2>/dev/null | head -1 || true)"
if [[ -z "$LATEST_CSV" || ! -f "$LATEST_CSV" ]]; then
  echo "No CSV produced by latency.py" >&2
  exit 1
fi
cp "$LATEST_CSV" "$CSV_OUT"
if [[ -f "$SCRIPT_DIR/latency.png" ]]; then
  cp "$SCRIPT_DIR/latency.png" "$PNG_OUT"
fi
echo "latency_csv=$CSV_OUT" >> "$RESOURCES_OUT"
echo "latency_png=$PNG_OUT" >> "$RESOURCES_OUT"
echo "source_csv=$LATEST_CSV" >> "$RESOURCES_OUT"

SCRIPT_END=$(date +%s.%N)
TOTAL_ELAPSED=$(python3 -c "print(round($SCRIPT_END - $SCRIPT_START, 2))")
DISK_KB=$(du -sk "$RESULTS_DIR" "$SCRIPT_DIR/data" 2>/dev/null | awk '{s+=$1} END {print s}')

{
  echo ""
  echo "# Totals (for ARTIFACT_GPU.md — paper GPU latency claim)"
  echo "total_elapsed_s=$TOTAL_ELAPSED"
  echo "peak_rss_kb=$TOTAL_RSS_KB"
  echo "artifact_disk_kb=$DISK_KB"
} | tee -a "$RESOURCES_OUT"

echo ""
echo "Wrote: $CSV_OUT"
echo "       $PNG_OUT"
echo "       $RESOURCES_OUT"
