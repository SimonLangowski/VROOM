#!/usr/bin/env bash
# Install paper CPU toolchain: clang 21.1.0 + Google Benchmark (clang-built).
# Tested on Amazon Linux 2023, Intel Xeon Platinum 8488C (AVX-512 IFMA).
#
# Usage:
#   ./scripts/setup_toolchain.sh           # install / refresh toolchain
#   ./scripts/setup_toolchain.sh --print   # print export lines only (no install)
#   ./scripts/setup_toolchain.sh --help
#
# After running, source the exports (or add them to ~/.bashrc):
#   eval "$(./scripts/setup_toolchain.sh --print)"

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"

LLVM_VERSION=21.1.0
LLVM_COMMIT=3623fe661ae35c6c80ac221f14d85be76aa870f1
LLVM_TARBALL="LLVM-${LLVM_VERSION}-Linux-X64.tar.xz"
LLVM_URL="https://github.com/llvm/llvm-project/releases/download/llvmorg-${LLVM_VERSION}/${LLVM_TARBALL}"
LLVM_SHA256=4a8c4b07646a4a2eb76ccf1d73522c3e13519745b72ef09d631c09b7577b0ed2
LLVM_PREFIX="${LLVM_PREFIX:-$HOME/LLVM-${LLVM_VERSION}-Linux-X64}"

BENCHMARK_REPO=https://github.com/google/benchmark.git
BENCHMARK_COMMIT=04ccbd86038796c319ea19987457e651a24f6b44
BENCHMARK_DIR="${BENCHMARK_DIR:-$HOME/benchmark}"
BENCHMARK_BUILD="${BENCHMARK_BUILD:-$BENCHMARK_DIR/build2}"

print_exports() {
  cat <<EOF
export PATH="$LLVM_PREFIX/bin:\$PATH"
export CC=clang
export CXX=clang++
export BENCHMARK_LIBS="$BENCHMARK_BUILD/src/libbenchmark.a -lpthread"
export BENCHMARK_INC="-I$BENCHMARK_DIR/include"
EOF
}

usage() {
  cat <<'EOF'
Install paper CPU toolchain: clang 21.1.0 + Google Benchmark (clang-built).

Usage:
  ./scripts/setup_toolchain.sh           # install / refresh toolchain
  ./scripts/setup_toolchain.sh --print   # print export lines only (no install)
  ./scripts/setup_toolchain.sh --help

After running:
  eval "$(./scripts/setup_toolchain.sh --print)"
EOF
  echo ""
  echo "Pinned versions (reference host, Amazon Linux 2023):"
  echo "  LLVM ${LLVM_VERSION} prebuilt  llvm-project @ ${LLVM_COMMIT}"
  echo "  google/benchmark             @ ${BENCHMARK_COMMIT}  (v1.9.4-89-g04ccbd8)"
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi

if [[ "${1:-}" == "--print" ]]; then
  print_exports
  exit 0
fi

if [[ "${1:-}" != "" ]]; then
  echo "Unknown option: $1" >&2
  usage >&2
  exit 2
fi

echo "== System packages (Amazon Linux 2023) =="
if command -v dnf >/dev/null 2>&1; then
  sudo dnf install -y gcc gcc-c++ make gmp-devel cmake git curl
elif command -v apt-get >/dev/null 2>&1; then
  sudo apt-get update
  sudo apt-get install -y build-essential libgmp-dev cmake git curl
else
  echo "    (no dnf/apt — install gcc, g++, make, gmp-devel, cmake, git, curl manually)"
fi

echo "== LLVM ${LLVM_VERSION} prebuilt @ ${LLVM_COMMIT} =="
if [[ -x "$LLVM_PREFIX/bin/clang++" ]]; then
  echo "    already installed at $LLVM_PREFIX"
else
  tmpdir="$(mktemp -d)"
  trap 'rm -rf "$tmpdir"' EXIT
  curl -fsSL -o "$tmpdir/$LLVM_TARBALL" "$LLVM_URL"
  echo "    sha256 $LLVM_SHA256  $LLVM_TARBALL"
  echo "$LLVM_SHA256  $tmpdir/$LLVM_TARBALL" | sha256sum -c -
  tar xf "$tmpdir/$LLVM_TARBALL" -C "$HOME"
  echo "    extracted to $LLVM_PREFIX"
fi
export PATH="$LLVM_PREFIX/bin:$PATH"
export CC=clang
export CXX=clang++
"$CXX" --version | head -1

echo "== google/benchmark @ ${BENCHMARK_COMMIT} =="
if [[ ! -d "$BENCHMARK_DIR/.git" ]]; then
  git clone "$BENCHMARK_REPO" "$BENCHMARK_DIR"
fi
git -C "$BENCHMARK_DIR" fetch --tags origin
git -C "$BENCHMARK_DIR" checkout "$BENCHMARK_COMMIT"
echo "    checked out $(git -C "$BENCHMARK_DIR" rev-parse --short HEAD) \
($(git -C "$BENCHMARK_DIR" describe --tags --always))"

echo "== Build libbenchmark.a with clang++ =="
cmake -B "$BENCHMARK_BUILD" -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER="$CXX" \
  -DBENCHMARK_ENABLE_GTEST_TESTS=OFF
cmake --build "$BENCHMARK_BUILD"
test -f "$BENCHMARK_BUILD/src/libbenchmark.a"

echo ""
echo "Toolchain ready. Add to your shell (or ~/.bashrc):"
print_exports
echo ""
echo "Verify from repo root ($ROOT):"
echo "  eval \"\$(./scripts/setup_toolchain.sh --print)\""
echo "  ./scripts/smoke_test.sh"
echo "  ./scripts/reproduce_cpu_bench.sh"
