#!/usr/bin/env bash
# Minimal GPU verification: one CGBN + one RNS latency point with correctness checks.
# Exit non-zero on first failure.
#
# Usage (from repo root or scripts/):
#   ./scripts/smoke_test_gpu.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$SCRIPT_DIR"

CGBN_HDR="$ROOT/gpu/algorithms/baselines/CGBN/include/cgbn/cgbn.h"
if [[ ! -f "$CGBN_HDR" ]]; then
  echo "CGBN not found. Run: $ROOT/scripts/setup_gpu_deps.sh" >&2
  exit 1
fi

export MPLBACKEND=Agg

echo "== CGBN baseline (64-bit modulus, tpi=2) =="
python3 -c "
from bench_cgbn import run_cgbn_benchmark
from run_on_params import params_to_args
args = {
    'blocksize': 32, 'blocksper': 1, 'total': 16, 'btpi': 2,
    'modbits': 64, 'pow': 64, 'm': 8,
}
r = run_cgbn_benchmark(params_to_args(args), True)
assert r is not None and r.get('correct') is True, r
print('CGBN smoke OK:', r['time'], 'ms')
"

echo "== RNS Montgomery GPU (MontRNSReduceVec, 60-bit, stride=1) =="
python3 -c "
from gpu import run_rns_latency_benchmark, VECTOR
from run_on_params import params_to_args
args = {
    'blocksize': 32, 'blocksper': 1, 'total': 8, 'btpi': 4,
    'modbits': 60, 'pow': 60, 'm': 1, 'variation': VECTOR, 'stride': 1,
}
r = run_rns_latency_benchmark(params_to_args(args), True)
assert r is not None and r.get('correct') is True, r
print('RNS GPU smoke OK:', r['time'], 'ms')
"

echo "== smoke_test_gpu OK =="
