// BLS12-381 benchmarks (8×50-bit RNS): full bench_pairing_50bit suite, Matrix vs MatrixNoK.
#include <benchmark/benchmark.h>
#include "bench_bls12_381_impl.hpp"
#include "bls12_381_ring_params.hpp"

using namespace bls12_381_params;
using namespace bench_bls12_381;

#define BENCH_RING(name, impl)                          \
    static void BM_##name##_Matrix(benchmark::State &s) { \
        impl<BlsRing>(s);                               \
        s.SetLabel(#name "_Matrix");                    \
    }                                                   \
    static void BM_##name##_MatrixNoK(benchmark::State &s) { \
        impl<BlsRingMatrixNoK>(s);                      \
        s.SetLabel(#name "_MatrixNoK");                 \
    }                                                   \
    BENCHMARK(BM_##name##_Matrix);                      \
    BENCHMARK(BM_##name##_MatrixNoK)

#define BENCH_RING_BATCH(name, impl, n)                 \
    static void BM_##name##_Matrix_##n(benchmark::State &s) { \
        impl<n, BlsRing>(s);                            \
        s.SetLabel(#name "_Matrix_" #n);                \
    }                                                   \
    static void BM_##name##_MatrixNoK_##n(benchmark::State &s) { \
        impl<n, BlsRingMatrixNoK>(s);                   \
        s.SetLabel(#name "_MatrixNoK_" #n);             \
    }                                                   \
    BENCHMARK(BM_##name##_Matrix_##n);                  \
    BENCHMARK(BM_##name##_MatrixNoK_##n)

#define BENCH_RING_EC(name, impl)                       \
    static void BM_##name##_Matrix(benchmark::State &s) { \
        impl<BlsRingEC>(s);                             \
        s.SetLabel(#name "_Matrix");                    \
    }                                                   \
    static void BM_##name##_MatrixNoK(benchmark::State &s) { \
        impl<BlsRingMatrixNoKEC>(s);                      \
        s.SetLabel(#name "_MatrixNoK");                 \
    }                                                   \
    BENCHMARK(BM_##name##_Matrix);                      \
    BENCHMARK(BM_##name##_MatrixNoK)

BENCH_RING(ModMul, BM_ModMul_Impl);
BENCH_RING_BATCH(BatchModMul, BM_BatchModMul_Impl, 1);
BENCH_RING_BATCH(BatchModMul, BM_BatchModMul_Impl, 2);
BENCH_RING_BATCH(BatchModMul, BM_BatchModMul_Impl, 3);
BENCH_RING_BATCH(BatchModMul, BM_BatchModMul_Impl, 4);
BENCH_RING_BATCH(BatchModMul, BM_BatchModMul_Impl, 5);
BENCH_RING_BATCH(BatchModMul, BM_BatchModMul_Impl, 6);
BENCH_RING_BATCH(BatchModMul, BM_BatchModMul_Impl, 7);
BENCH_RING_BATCH(BatchModMul, BM_BatchModMul_Impl, 8);
BENCH_RING_BATCH(BatchModMul, BM_BatchModMul_Impl, 10);
BENCH_RING_BATCH(BatchModMul, BM_BatchModMul_Impl, 12);

BENCH_RING_BATCH(BatchReduceExpand, BM_BatchReduceExpand_Impl, 1);
BENCH_RING_BATCH(BatchReduceExpand, BM_BatchReduceExpand_Impl, 2);
BENCH_RING_BATCH(BatchReduceExpand, BM_BatchReduceExpand_Impl, 3);
BENCH_RING_BATCH(BatchReduceExpand, BM_BatchReduceExpand_Impl, 4);
BENCH_RING_BATCH(BatchReduceExpand, BM_BatchReduceExpand_Impl, 5);
BENCH_RING_BATCH(BatchReduceExpand, BM_BatchReduceExpand_Impl, 6);
BENCH_RING_BATCH(BatchReduceExpand, BM_BatchReduceExpand_Impl, 12);

BENCH_RING(FP2_Mul, BM_FP2_Mul_Impl);
BENCH_RING(FP2_Sqr, BM_FP2_Sqr_Impl);

BENCH_RING_EC(G1_Add, BM_G1_Add_Impl);
BENCH_RING_EC(G1_MixedAdd, BM_G1_MixedAdd_Impl);
BENCH_RING_EC(G1_Double, BM_G1_Double_Impl);
BENCH_RING_EC(G2_Add, BM_G2_Add_Impl);
BENCH_RING_EC(G2_MixedAdd, BM_G2_MixedAdd_Impl);
BENCH_RING_EC(G2_Double, BM_G2_Double_Impl);

BENCH_RING_EC(G1_ScalarMult_GLV, BM_G1_ScalarMult_GLV_Impl);
BENCH_RING_EC(G2_ScalarMult_GLS, BM_G2_ScalarMult_GLS_Impl);

BENCH_RING(FP12_Mul, BM_FP12_Mul_Impl);
BENCH_RING(FP12_Sqr, BM_FP12_Sqr_Impl);
BENCH_RING(FP12_CyclotomicSqr, BM_FP12_CyclotomicSqr_Impl);
BENCH_RING(FP12_Conjugate, BM_FP12_Conjugate_Impl);
BENCH_RING(FP12_FrobeniusMap1, BM_FP12_FrobeniusMap1_Impl);
BENCH_RING(FP12_FrobeniusMap2, BM_FP12_FrobeniusMap2_Impl);
BENCH_RING(FP12_FrobeniusMap3, BM_FP12_FrobeniusMap3_Impl);
BENCH_RING(FP12_MulByXy0z00, BM_FP12_MulByXy0z00_Impl);
BENCH_RING(FP12_FinalExpAfterInverse, BM_FP12_FinalExpAfterInverse_Impl);

BENCH_RING(Inverse, BM_Inverse_Impl);

BENCH_RING(MillerLoop_DoubleStep, BM_MillerLoop_DoubleStep_Impl);
BENCH_RING(MillerLoop_AddStep, BM_MillerLoop_AddStep_Impl);
BENCH_RING(MillerLoop, BM_MillerLoop_Impl);

BENCH_RING(FP12_Inverse_BLST, BM_FP12_Inverse_BLST_Impl);
BENCH_RING(FP12_Inverse_RNS_BLST, BM_FP12_Inverse_RNS_BLST_Impl);
BENCH_RING(FP12_FinalExp_BLST_Inverter, BM_FP12_FinalExp_BLST_Inverter_Impl);
BENCH_RING(FP12_FinalExp_RNS_BLST_Inverter, BM_FP12_FinalExp_RNS_BLST_Inverter_Impl);
BENCH_RING(Pairing_BLST_Inverter, BM_Pairing_BLST_Inverter_Impl);
BENCH_RING(Pairing_RNS_BLST_Inverter, BM_Pairing_RNS_BLST_Inverter_Impl);

BENCHMARK_MAIN();
