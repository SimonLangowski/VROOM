# RNS Montgomery multiplication with AVX512

Standalone C++ template/header code. Benchmarks use Google Benchmark.  Inversion with BLST requires BLST with extern "C" modification.

## Building & testing

```bash
cd blst
make
cd ../src
make
```
Note that making the fallback version takes a while.  The AVX version is much faster to compile.

- **`make bench_pairing_50bit`** — pairing benchmark with 50-bit parameters (BLS12-381, AVX512 IFMA).
- **`make test_pairing_avx`** — pairing test (AVX512 IFMA).
- **`make test_pairing`** — pairing test (integer fallback, for machines without AVX512 IFMA, testing correctness only).

## Running
To see speedups and develop, we used the AWS `c7i-flex.large` machine.
To match the paper numbers, use the AWS `c7i.metal-24xl` machine, and the `clang-21` compiler.
GCC seems to be significantly slower than clang, and most modern clang versions perform similarly.
- **./bench_pairing_50bit**
Sample output on `c7i.metal-24xl` with `clang-21`:
```


```

## License

MIT. The **blst/** subtree keeps its original Apache license and copyright (see `blst/LICENSE`).

## Security

This is research-grade code, demonstrating the speedups of combining AVX512-IFMA, RNS Montgomery, and the BLS12-381 curve.  
It has not been vetted or audited for security.