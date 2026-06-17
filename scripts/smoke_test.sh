#!/usr/bin/env bash
# Quick artifact verification. Exit non-zero on first failure.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

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

echo "== smoke test OK =="
