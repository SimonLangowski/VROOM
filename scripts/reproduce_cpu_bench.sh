#!/usr/bin/env bash
# Reproduce CPU benchmark table (paper numbers): build, run bench_bls12_381 +
# bench_blst_complete, write JSON, parse with RNS/BLST speedup ratios.
#
# Usage:
#   ./scripts/reproduce_cpu_bench.sh              # AVX-512 IFMA (production)
#   FALLBACK=1 ./scripts/reproduce_cpu_bench.sh   # integer fallback (slow compile)
#   BENCHMARK_LIBS=/path/to/libbenchmark.a\ -lpthread ./scripts/reproduce_cpu_bench.sh
#   ./scripts/reproduce_cpu_bench.sh --exclude-nok  # Matrix rows only in parsed table

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

FALLBACK="${FALLBACK:-0}"
EXCLUDE_NOK=0
BENCH_BASE="${BENCH_BASE:-bench_bls12_381}"
BLST_BENCH_BIN="${BLST_BENCH_BIN:-bench_blst_complete}"
MIN_TIME="${BENCH_MIN_TIME:-0.5s}"
RESULTS_DIR="${RESULTS_DIR:-$ROOT/results}"

if [[ "$FALLBACK" == "1" ]]; then
  BENCH_TARGET="${BENCH_TARGET:-${BENCH_BASE}_fallback}"
  BENCH_BIN="${BENCH_BIN:-${BENCH_BASE}_fallback}"
else
  BENCH_TARGET="${BENCH_TARGET:-$BENCH_BASE}"
  BENCH_BIN="${BENCH_BIN:-$BENCH_BASE}"
fi

JSON_OUT="$RESULTS_DIR/${BENCH_BIN}.json"
BLST_JSON_OUT="$RESULTS_DIR/${BLST_BENCH_BIN}.json"
TABLE_OUT="$RESULTS_DIR/${BENCH_BIN}_table.txt"
COMPARE_OUT="$RESULTS_DIR/${BENCH_BIN}_vs_blst.txt"
RESOURCES_OUT="$RESULTS_DIR/${BENCH_BIN}_resources.txt"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --exclude-nok) EXCLUDE_NOK=1; shift ;;
    --bench)
      BENCH_BASE="$2"
      BENCH_TARGET="$2"
      BENCH_BIN="$2"
      JSON_OUT="$RESULTS_DIR/${BENCH_BIN}.json"
      TABLE_OUT="$RESULTS_DIR/${BENCH_BIN}_table.txt"
      COMPARE_OUT="$RESULTS_DIR/${BENCH_BIN}_vs_blst.txt"
      RESOURCES_OUT="$RESULTS_DIR/${BENCH_BIN}_resources.txt"
      shift 2
      ;;
    -h|--help)
      sed -n '2,12p' "$0"
      exit 0
      ;;
    *) echo "Unknown option: $1" >&2; exit 2 ;;
  esac
done

MAKE_BENCH_EXTRA=()
if [[ -n "${BENCHMARK_LIBS:-}" ]]; then
  MAKE_BENCH_EXTRA+=(BENCHMARK_LIBS="$BENCHMARK_LIBS")
fi
if [[ -n "${BENCHMARK_INC:-}" ]]; then
  MAKE_BENCH_EXTRA+=(BENCHMARK_INC="$BENCHMARK_INC")
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
  echo "CPU benchmark reproduction"
  echo "date: $(date -Is)"
  echo "host: $(uname -n)"
  echo "fallback: $FALLBACK"
  echo "rns_bench: $BENCH_BIN"
  echo "blst_bench: $BLST_BENCH_BIN"
  echo ""
} | tee "$RESOURCES_OUT"

echo "== Phase: build blst (+ bench_blst_complete) =="
run_timed "build blst" make -C blst gbench "${MAKE_BENCH_EXTRA[@]}"
echo "build_blst_elapsed_s=$PHASE_ELAPSED" >> "$RESOURCES_OUT"
echo "build_blst_max_rss_kb=$PHASE_RSS_KB" >> "$RESOURCES_OUT"
TOTAL_RSS_KB=$(( TOTAL_RSS_KB > PHASE_RSS_KB ? TOTAL_RSS_KB : PHASE_RSS_KB ))

echo "== Phase: build RNS benchmark =="
if [[ "$FALLBACK" == "1" ]]; then
  run_timed "build bench (fallback)" make -C src -f Makefile.fallback "${MAKE_BENCH_EXTRA[@]}" "$BENCH_TARGET"
else
  run_timed "build bench" make -C src "${MAKE_BENCH_EXTRA[@]}" "$BENCH_TARGET"
fi
echo "build_bench_elapsed_s=$PHASE_ELAPSED" >> "$RESOURCES_OUT"
echo "build_bench_max_rss_kb=$PHASE_RSS_KB" >> "$RESOURCES_OUT"
TOTAL_RSS_KB=$(( TOTAL_RSS_KB > PHASE_RSS_KB ? TOTAL_RSS_KB : PHASE_RSS_KB ))

BENCH_PATH="$ROOT/src/$BENCH_BIN"
BLST_BENCH_PATH="$ROOT/blst/$BLST_BENCH_BIN"
if [[ ! -x "$BENCH_PATH" ]]; then
  echo "RNS benchmark binary not found: $BENCH_PATH" >&2
  exit 1
fi
if [[ ! -x "$BLST_BENCH_PATH" ]]; then
  echo "BLST benchmark binary not found: $BLST_BENCH_PATH" >&2
  exit 1
fi

echo "== Phase: run RNS benchmark (JSON) =="
run_timed "rns benchmark" "$BENCH_PATH" \
  --benchmark_out="$JSON_OUT" \
  --benchmark_out_format=json \
  --benchmark_min_time="$MIN_TIME"
echo "rns_bench_elapsed_s=$PHASE_ELAPSED" >> "$RESOURCES_OUT"
echo "rns_bench_max_rss_kb=$PHASE_RSS_KB" >> "$RESOURCES_OUT"
TOTAL_RSS_KB=$(( TOTAL_RSS_KB > PHASE_RSS_KB ? TOTAL_RSS_KB : PHASE_RSS_KB ))
echo "rns_json=$JSON_OUT" >> "$RESOURCES_OUT"

echo "== Phase: run BLST benchmark (JSON) =="
run_timed "blst benchmark" "$BLST_BENCH_PATH" \
  --benchmark_out="$BLST_JSON_OUT" \
  --benchmark_out_format=json \
  --benchmark_min_time="$MIN_TIME"
echo "blst_bench_elapsed_s=$PHASE_ELAPSED" >> "$RESOURCES_OUT"
echo "blst_bench_max_rss_kb=$PHASE_RSS_KB" >> "$RESOURCES_OUT"
TOTAL_RSS_KB=$(( TOTAL_RSS_KB > PHASE_RSS_KB ? TOTAL_RSS_KB : PHASE_RSS_KB ))
echo "blst_json=$BLST_JSON_OUT" >> "$RESOURCES_OUT"

echo "== Phase: parse + compare =="
PARSE_ARGS=(
  "$ROOT/scripts/parse_bench_json.py" "$JSON_OUT"
  --blst "$BLST_JSON_OUT"
  -o "$TABLE_OUT"
  --compare-output "$COMPARE_OUT"
)
if [[ "$EXCLUDE_NOK" == "1" ]]; then
  PARSE_ARGS+=(--exclude-nok)
fi
run_timed "parse" python3 "${PARSE_ARGS[@]}"
echo "parse_elapsed_s=$PHASE_ELAPSED" >> "$RESOURCES_OUT"
echo "table_out=$TABLE_OUT" >> "$RESOURCES_OUT"
echo "compare_out=$COMPARE_OUT" >> "$RESOURCES_OUT"

SCRIPT_END=$(date +%s.%N)
TOTAL_ELAPSED=$(python3 -c "print(round($SCRIPT_END - $SCRIPT_START, 2))")
DISK_KB=$(du -sk "$RESULTS_DIR" "$ROOT/src/$BENCH_BIN" "$ROOT/blst/$BLST_BENCH_BIN" "$ROOT/blst/libblst_small.a" 2>/dev/null | awk '{s+=$1} END {print s}')

{
  echo ""
  echo "# Totals (for ARTIFACT.md — one row per paper claim)"
  echo "total_elapsed_s=$TOTAL_ELAPSED"
  echo "peak_rss_kb=$TOTAL_RSS_KB"
  echo "artifact_disk_kb=$DISK_KB"
} | tee -a "$RESOURCES_OUT"

echo ""
echo "Wrote: $JSON_OUT"
echo "       $BLST_JSON_OUT"
echo "       $TABLE_OUT"
echo "       $COMPARE_OUT"
echo "       $RESOURCES_OUT"
