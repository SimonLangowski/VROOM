# blst

Modified copy of the BLST code (Supranational; Apache license, allowing modification — see this folder).

- **Modifications:** `extern "C"` visibility for C++ compatibility; extra wrappers for finer-grained benchmarking.
- **`empirical_basis.cpp`** — derives constants for “flattened” operations (useful in RNS), effectively removing Karatsuba (desirable for RNS).
- **`flattener.py`** — turns output polynomials into C++ code.
