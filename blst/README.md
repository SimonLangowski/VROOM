# blst

Modified copy of the BLST code (Supranational; Apache license, allowing modification — see this folder).

**Upstream baseline (paper):** [supranational/blst](https://github.com/supranational/blst) @ **`8065152`**. The vendored tree in this directory includes local modifications on top of that revision; use it for benchmarks (do not replace with a fresh upstream clone).

- **Modifications:** `extern "C"` visibility for C++ compatibility; extra wrappers for finer-grained benchmarking.
- **`empirical_basis.cpp`** — derives constants for “flattened” operations (useful in RNS), effectively removing Karatsuba (desirable for RNS).
- **`flattener.py`** — turns output polynomials into C++ code.
