# RNS Montgomery multiplication with AVX512

Standalone C++ template/header code. Benchmarks use Google Benchmark; Python scripts use numpy, sympy, sagemath, and matplotlib.

## Building & testing

```bash
cd src
make
```
Note that making the fallback version takes a while.  The AVX version is much faster to compile.

- **`make bench_pairing_50bit`** — pairing benchmark with 50-bit parameters (BLS12-381, AVX).
- **`make test_pairing`** — pairing test (integer fallback).
- **`make test_pairing_avx`** — pairing test (AVX).

## License

MIT. The **blst/** subtree keeps its original Apache license and copyright (see `blst/LICENSE`).
