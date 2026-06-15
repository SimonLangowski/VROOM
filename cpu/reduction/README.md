# cpu/reduction

Per-lane Montgomery reduction used by the vector CRNS layer.

- **`montgomery.hpp`** — `MontgomeryReduce` / reduction helpers included from `cpu/vector/conversion.hpp` and `cpu/vector/multiplication.hpp`.

The standalone reduction templates here support RNS matrix operations; the hot path is selected via `cpu/vector/vector_impl.hpp` (AVX-512 IFMA or integer fallback).
