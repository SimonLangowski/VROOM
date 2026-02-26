# RNS Montgomery multiplication with AVX512

Standalone C++ template/header code. Benchmarks use Google Benchmark.  Inversion with BLST requires BLST with extern "C" modification.

## Building & testing

```bash
cd src
make
```
Note that making the fallback version takes a while.  The AVX version is much faster to compile.

- **`make bench_pairing_50bit`** — pairing benchmark with 50-bit parameters (BLS12-381, AVX).
- **`make test_pairing_avx`** — pairing test (AVX512).
- **`make test_pairing`** — pairing test (integer fallback, for machines without AVX IFMA, testing correctness only).

## License

MIT. The **blst/** subtree keeps its original Apache license and copyright (see `blst/LICENSE`).

## Security

This is research-grade code, demonstrating the speedups of combining AVX512-IFMA, RNS Montgomery, and the BLS12-381 curve.  
It has not been vetted or audited for security.