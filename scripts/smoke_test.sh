#!/usr/bin/env bash
# Quick artifact verification. Exit non-zero on first failure.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

# Paper toolchain (clang + source-built benchmark). Auto-load if not already exported.
if [[ -z "${BENCHMARK_INC:-}" || -z "${BENCHMARK_LIBS:-}" ]]; then
  eval "$(bash "$ROOT/scripts/setup_toolchain.sh" --print)"
fi
bench_inc_dir="${BENCHMARK_INC#-I}"
if [[ ! -f "$bench_inc_dir/benchmark/benchmark.h" ]]; then
  echo "Google Benchmark not found at $bench_inc_dir" >&2
  echo "Run: ./scripts/setup_toolchain.sh" >&2
  exit 1
fi

FALLBACK="${FALLBACK:-0}"

echo "== BLST =="
make -C blst

echo "== src core tests =="
if [[ "$FALLBACK" == "1" ]]; then
  make -C src -f Makefile.fallback test
else
  make -C src test
  make -C src test-bls12-381-field-mul
fi

echo "== examples =="
if [[ "$FALLBACK" == "1" ]]; then
  make -C examples FALLBACK=1
else
  make -C examples
fi

echo "== run examples =="
./examples/01_parameter_setup | head -20
./examples/02_singular_modmul
./examples/03_sum_of_products
./examples/04_difference_of_products
./examples/05_product_of_sums

echo "== smoke test OK =="
