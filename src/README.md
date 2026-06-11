# src

Main C++ code (templates/headers, tests, benchmarks). Build with `make`; see repo root **README.md** for targets.

Update for camera ready: added increased performance via quadratic residue decomposition.
Decomposition is computed in python `rns_secp256k1.py` and then parameters are output to a .hpp file, and imported in the corresponding int4rnssystem, perm_batch, perm_params, ring_params, and perm files.

- **`*.hpp`** — RNS Montgomery, pairing, EC, scalar mult.
- **`test_*.cpp`** — unit tests.
- **`bench_*.cpp`** — Google Benchmark harnesses.
- **`gmp_wrapper.hpp`** — pulls in **`../cpu/precompute/gmp_wrapper.hpp`** and **`../test_data/test_values.hpp`**.
