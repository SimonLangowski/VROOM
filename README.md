# RNS Montgomery multiplication with AVX512

Standalone C++ template/header code. Benchmarks use Google Benchmark. Inversion with BLST requires BLST with extern "C" modification.

## Building & testing

```bash
cd blst
make
cd ../src
make
```
Note that making the fallback version takes a while. The AVX version is much faster to compile.

- **`make bench_pairing_50bit`** — pairing benchmark with 50-bit parameters (BLS12-381, AVX512 IFMA).
- **`make test_pairing_avx`** — pairing test (AVX512 IFMA).
- **`make test_pairing`** — pairing test (integer fallback, for machines without AVX512 IFMA, testing correctness only).

## Running

To see speedups and develop, we used the AWS `c7i-flex.large` machine.
To match the paper numbers, use the AWS `c7i.metal-24xl` machine, and the `clang-21` compiler.
GCC seems to be significantly slower than clang, and most modern clang versions perform similarly.

- **`./bench_pairing_50bit`**

Sample output on `c7i.metal-24xl` with `clang-21`:

```
Running ./bench_pairing_50bit
Run on (96 X 800 MHz CPU s)
CPU Caches:
  L1 Data 48 KiB (x48)
  L1 Instruction 32 KiB (x48)
  L2 Unified 2048 KiB (x48)
  L3 Unified 107520 KiB (x1)
Load Average: 0.68, 0.45, 0.18
-----------------------------------------------------------------------------
Benchmark                                   Time             CPU   Iterations
-----------------------------------------------------------------------------
BM_ModMul                                28.4 ns         28.4 ns     24694133 ModMul
BM_BatchModMul_1                         28.3 ns         28.3 ns     24767606 BatchModMul_1
BM_BatchModMul_2                         46.6 ns         46.6 ns     15026360 BatchModMul_2
BM_BatchModMul_3                         62.1 ns         62.1 ns     11180624 BatchModMul_3
BM_BatchModMul_4                         77.4 ns         77.4 ns      9037047 BatchModMul_4
BM_BatchModMul_5                         88.3 ns         88.3 ns      7922070 BatchModMul_5
BM_BatchModMul_6                          109 ns          109 ns      6391533 BatchModMul_6
BM_BatchModMul_7                          115 ns          115 ns      6062577 BatchModMul_7
BM_BatchModMul_8                          132 ns          132 ns      5286519 BatchModMul_8
BM_BatchModMul_10                         158 ns          158 ns      4416745 BatchModMul_10
BM_BatchModMul_12                         207 ns          207 ns      3378190 BatchModMul_12
BM_BatchReduceExpand_1                   20.9 ns         20.9 ns     33528312 BatchReduceExpand_1
BM_BatchReduceExpand_2                   36.2 ns         36.2 ns     19324062 BatchReduceExpand_2
BM_BatchReduceExpand_3                   50.4 ns         50.4 ns     13850523 BatchReduceExpand_3
BM_BatchReduceExpand_4                   65.6 ns         65.6 ns     10681823 BatchReduceExpand_4
BM_BatchReduceExpand_5                   78.1 ns         78.1 ns      8946299 BatchReduceExpand_5
BM_BatchReduceExpand_6                   94.3 ns         94.3 ns      7411180 BatchReduceExpand_6
BM_BatchReduceExpand_12                   196 ns          196 ns      3620070 BatchReduceExpand_12
BM_FP2_Mul                               46.9 ns         46.9 ns     14260145 FP2_Mul_50bit
BM_FP2_Sqr                               49.3 ns         49.3 ns     14191739 FP2_Sqr_50bit
BM_G1_Add                                 174 ns          174 ns      4016688 G1_Add
BM_G1_MixedAdd                            168 ns          168 ns      4170086 G1_MixedAdd
BM_G1_Double                              169 ns          169 ns      4142371 G1_Double
BM_G2_Add                                 356 ns          356 ns      1965909 G2_Add
BM_G2_MixedAdd                            311 ns          311 ns      2250196 G2_MixedAdd
BM_G2_Double                              259 ns          259 ns      2700671 G2_Double
BM_G1_ScalarMult_GLV_50bit              35893 ns        35893 ns        19501 G1_ScalarMult_GLV_50bit
BM_G2_ScalarMult_GLS_50bit              48104 ns        48104 ns        14553 G2_ScalarMult_GLS_50bit
BM_FP12_Mul                               299 ns          299 ns      2338161 FP12_Mul
BM_FP12_Sqr                               265 ns          265 ns      2644300 FP12_Sqr
BM_FP12_CyclotomicSqr                     238 ns          238 ns      2939282 FP12_CyclotomicSqr
BM_FP12_Conjugate                        8.85 ns         8.85 ns     79075566 FP12_Conjugate
BM_FP12_FrobeniusMap1                     161 ns          161 ns      4335378 FP12_FrobeniusMap1
BM_FP12_FrobeniusMap2                     133 ns          133 ns      5269124 FP12_FrobeniusMap2
BM_FP12_FrobeniusMap3                     104 ns          104 ns      6728849 FP12_FrobeniusMap3
BM_FP12_MulByXy0z00                       253 ns          253 ns      2766511 FP12_MulByXy0z00
BM_FP12_FinalExpAfterInverse            88630 ns        88629 ns         7891 FP12_FinalExpAfterInverse
BM_Inverse                              20856 ns        20856 ns        33563 Inverse
BM_MillerLoop_DoubleStep                  622 ns          622 ns      1125633 MillerLoop_DoubleStep
BM_MillerLoop_AddStep                     688 ns          688 ns      1017913 MillerLoop_AddStep
BM_MillerLoop                           58985 ns        58985 ns        11864 MillerLoop
BM_FP12_Inverse_BLST                     7807 ns         7807 ns        89491 FP12_Inverse_BLST_50bit
BM_FP12_Inverse_RNS_BLST                 3578 ns         3578 ns       195749 FP12_Inverse_RNS_BLST_50bit
BM_FP12_FinalExp_BLST_Inverter          96778 ns        96778 ns         7247 FP12_FinalExp_BLST_Inverter_50bit
BM_FP12_FinalExp_RNS_BLST_Inverter      92103 ns        92104 ns         7596 FP12_FinalExp_RNS_BLST_Inverter_50bit
BM_Pairing_BLST_Inverter               156004 ns       156005 ns         4485 Pairing_BLST_Inverter_50bit
BM_Pairing_RNS_BLST_Inverter           151537 ns       151534 ns         4620 Pairing_RNS_BLST_Inverter_50bit
```

- **BM_BatchModMul_N** — batch Montgomery multiplication, `a*b mod p`.
- **RNS_BLST_Inverter** — uses RNS for FP12 inversion and BLST for FP inversion.
- **BM_Inverse** — uses an addition chain of length 425 to invert in RNS without BLST dependency (compare to **BM_FP12_Inverse_BLST** and  **BM_FP12_Inverse_RNS_BLST**)

## License

MIT. The **blst/** subtree keeps its original Apache license and copyright (see `blst/LICENSE`).

## Security

This is research-grade code, demonstrating the speedups of combining AVX512-IFMA, RNS Montgomery, and the BLS12-381 curve.
It has not been vetted or audited for security.