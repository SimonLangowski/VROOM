# Shared wall-time / peak-RSS helpers for baselines/bench*.sh
# Source from bench scripts; do not execute directly.

bench_init() {
  local tag="$1"
  BASELINES_DIR="$(cd "$(dirname "${BASH_SOURCE[1]}")" && pwd)"
  ROOT="$(cd "$BASELINES_DIR/.." && pwd)"
  RESULTS_DIR="${RESULTS_DIR:-$ROOT/results}"
  mkdir -p "$RESULTS_DIR"
  RESOURCES_OUT="${RESOURCES_OUT:-$RESULTS_DIR/baselines_${tag}_resources.txt}"
  BENCH_LOG="${BENCH_LOG:-$RESULTS_DIR/baselines_${tag}_bench.log}"
  TOTAL_RSS_KB=0
  SCRIPT_START=$(date +%s.%N)

  {
    echo "Baselines benchmark ($tag)"
    echo "date: $(date -Is)"
    echo "host: $(uname -n)"
    echo "rustc: $(rustc --version 2>/dev/null || echo missing)"
    echo "cargo: $(cargo --version 2>/dev/null || echo missing)"
    echo "go: $(go version 2>/dev/null || echo missing)"
    echo "bench_log=$BENCH_LOG"
    echo ""
  } | tee "$RESOURCES_OUT"
  : >"$BENCH_LOG"
}

run_timed() {
  local label="$1"
  shift
  local key
  key="$(echo "$label" | tr '[:upper:]' '[:lower:]' | tr ' /' '__' | tr -cd '[:alnum:]_')"
  local tf
  tf="$(mktemp)"
  local t0 t1
  t0=$(date +%s.%N)
  set +e
  /usr/bin/time -f '%M' -o "$tf" "$@" 2>&1 | tee -a "$BENCH_LOG"
  local exit_code=${PIPESTATUS[0]}
  set -e
  t1=$(date +%s.%N)
  PHASE_RSS_KB=$(cat "$tf")
  rm -f "$tf"
  PHASE_ELAPSED=$(python3 -c "print(round($t1 - $t0, 2))")
  echo "  [$label] elapsed=${PHASE_ELAPSED}s  max_rss=${PHASE_RSS_KB}KB" | tee -a "$RESOURCES_OUT"
  echo "${key}_elapsed_s=$PHASE_ELAPSED" >>"$RESOURCES_OUT"
  echo "${key}_max_rss_kb=$PHASE_RSS_KB" >>"$RESOURCES_OUT"
  TOTAL_RSS_KB=$((TOTAL_RSS_KB > PHASE_RSS_KB ? TOTAL_RSS_KB : PHASE_RSS_KB))
  return "$exit_code"
}

bench_finish() {
  local tag="$1"
  SCRIPT_END=$(date +%s.%N)
  TOTAL_ELAPSED=$(python3 -c "print(round($SCRIPT_END - $SCRIPT_START, 2))")
  local disk_kb
  disk_kb=$(du -sk "$BASELINES_DIR" 2>/dev/null | awk '{print $1}')
  {
    echo ""
    echo "# Totals (for ARTIFACT.md — optional baseline reproduction claim)"
    echo "total_elapsed_s=$TOTAL_ELAPSED"
    echo "peak_rss_kb=$TOTAL_RSS_KB"
    echo "baselines_disk_kb=$disk_kb"
    echo "bench_log=$BENCH_LOG"
  } | tee -a "$RESOURCES_OUT"
  echo ""
  echo "Wrote: $RESOURCES_OUT"
  echo "       $BENCH_LOG"
}
