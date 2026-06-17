#!/usr/bin/env bash
# Optional: run external CPU baselines with per-phase wall time and peak RSS.
#
# Usage:
#   ./baselines/reproduce_baselines.sh           # BLS12-381 + BN254
#   ./baselines/reproduce_baselines.sh bls12     # BLS12-381 only
#   ./baselines/reproduce_baselines.sh bn254     # BN254 only
#
# Outputs under results/:
#   baselines_bls12_resources.txt, baselines_bls12_bench.log
#   baselines_bn254_resources.txt, baselines_bn254_bench.log
#
# Not required for the primary CPU artifact claim (see scripts/reproduce_cpu_bench.sh).

set -euo pipefail

BASE="$(cd "$(dirname "$0")" && pwd)"
CURVE="${1:-all}"

case "$CURVE" in
  bls12|bls12-381|bls)
    "$BASE/bench.sh"
    ;;
  bn254|bn)
    "$BASE/bench_bn254.sh"
    ;;
  all)
    "$BASE/bench.sh"
    "$BASE/bench_bn254.sh"
    ;;
  -h|--help)
    sed -n '2,14p' "$0"
    exit 0
    ;;
  *)
    echo "Unknown curve: $CURVE (use bls12, bn254, or all)" >&2
    exit 2
    ;;
esac
