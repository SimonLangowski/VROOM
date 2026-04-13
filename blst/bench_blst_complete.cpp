#include <benchmark/benchmark.h>
#include <cstring>

extern "C" {
#include "bindings/blst.h"
}

// ============================================================================
// G1 Operations
// ============================================================================

static void BM_G1_Add(benchmark::State& state) {
    blst_p1 p1, p2, result;
    const blst_p1 *gen = blst_p1_generator();
    
    blst_p1_add(&p1, gen, gen);  // 2*G
    blst_p1_add(&p2, gen, &p1);   // 3*G
    
    for (auto _ : state) {
        blst_p1_add(&result, &p1, &p2);
        benchmark::DoNotOptimize(result);
    }
}

static void BM_G1_Double(benchmark::State& state) {
    blst_p1 p1, result;
    const blst_p1 *gen = blst_p1_generator();
    
    blst_p1_add(&p1, gen, gen);  // 2*G
    
    for (auto _ : state) {
        blst_p1_double(&result, &p1);
        benchmark::DoNotOptimize(result);
    }
}

static void BM_G1_MixedAdd(benchmark::State& state) {
    blst_p1 p1, result;
    blst_p1_affine p2_affine;
    const blst_p1 *gen = blst_p1_generator();
    
    blst_p1_add(&p1, gen, gen);  // 2*G
    blst_p1_to_affine(&p2_affine, gen);  // Convert generator to affine
    
    for (auto _ : state) {
        blst_p1_add_affine(&result, &p1, &p2_affine);
        benchmark::DoNotOptimize(result);
    }
}

static void BM_G1_ScalarMul(benchmark::State& state) {
    blst_p1 p1, result;
    const blst_p1 *gen = blst_p1_generator();
    byte scalar[32];
    
    memset(scalar, 0x42, sizeof(scalar));
    scalar[0] = 0x01;
    
    blst_p1_add(&p1, gen, gen);  // 2*G
    
    for (auto _ : state) {
        scalar[31] = state.iterations() & 0xff;
        blst_p1_mult(&result, &p1, scalar, 256);
        benchmark::DoNotOptimize(result);
    }
}

// ============================================================================
// G2 Operations
// ============================================================================

static void BM_G2_Add(benchmark::State& state) {
    blst_p2 p1, p2, result;
    const blst_p2 *gen = blst_p2_generator();
    
    blst_p2_add(&p1, gen, gen);  // 2*G
    blst_p2_add(&p2, gen, &p1);   // 3*G
    
    for (auto _ : state) {
        blst_p2_add(&result, &p1, &p2);
        benchmark::DoNotOptimize(result);
    }
}

static void BM_G2_Double(benchmark::State& state) {
    blst_p2 p1, result;
    const blst_p2 *gen = blst_p2_generator();
    
    blst_p2_add(&p1, gen, gen);  // 2*G
    
    for (auto _ : state) {
        blst_p2_double(&result, &p1);
        benchmark::DoNotOptimize(result);
    }
}

static void BM_G2_MixedAdd(benchmark::State& state) {
    blst_p2 p1, result;
    blst_p2_affine p2_affine;
    const blst_p2 *gen = blst_p2_generator();
    
    blst_p2_add(&p1, gen, gen);  // 2*G
    blst_p2_to_affine(&p2_affine, gen);  // Convert generator to affine
    
    for (auto _ : state) {
        blst_p2_add_affine(&result, &p1, &p2_affine);
        benchmark::DoNotOptimize(result);
    }
}

static void BM_G2_ScalarMul(benchmark::State& state) {
    blst_p2 p1, result;
    const blst_p2 *gen = blst_p2_generator();
    byte scalar[32];
    
    memset(scalar, 0x42, sizeof(scalar));
    scalar[0] = 0x01;
    
    blst_p2_add(&p1, gen, gen);  // 2*G
    
    for (auto _ : state) {
        scalar[31] = state.iterations() & 0xff;
        blst_p2_mult(&result, &p1, scalar, 256);
        benchmark::DoNotOptimize(result);
    }
}

// ============================================================================
// Pairing Operations
// ============================================================================

static void BM_MillerLoop(benchmark::State& state) {
    blst_p1_affine p1_affine;
    blst_p2_affine p2_affine;
    blst_fp12 result;
    const blst_p1 *gen1 = blst_p1_generator();
    const blst_p2 *gen2 = blst_p2_generator();
    
    blst_p1_to_affine(&p1_affine, gen1);
    blst_p2_to_affine(&p2_affine, gen2);
    
    for (auto _ : state) {
        blst_miller_loop(&result, &p2_affine, &p1_affine);
        benchmark::DoNotOptimize(result);
    }
}

static void BM_FinalExp(benchmark::State& state) {
    blst_p1_affine p1_affine;
    blst_p2_affine p2_affine;
    blst_fp12 miller_result, result;
    const blst_p1 *gen1 = blst_p1_generator();
    const blst_p2 *gen2 = blst_p2_generator();
    
    blst_p1_to_affine(&p1_affine, gen1);
    blst_p2_to_affine(&p2_affine, gen2);
    blst_miller_loop(&miller_result, &p2_affine, &p1_affine);
    
    for (auto _ : state) {
        blst_final_exp(&result, &miller_result);
        benchmark::DoNotOptimize(result);
    }
}

static void BM_FullPairing(benchmark::State& state) {
    blst_p1_affine p1_affine;
    blst_p2_affine p2_affine;
    blst_fp12 miller_result, result;
    const blst_p1 *gen1 = blst_p1_generator();
    const blst_p2 *gen2 = blst_p2_generator();
    
    blst_p1_to_affine(&p1_affine, gen1);
    blst_p2_to_affine(&p2_affine, gen2);
    
    for (auto _ : state) {
        blst_miller_loop(&miller_result, &p2_affine, &p1_affine);
        blst_final_exp(&result, &miller_result);
        benchmark::DoNotOptimize(result);
    }
}

// ============================================================================
// Fp Operations
// ============================================================================

static void BM_Fp_Mul(benchmark::State& state) {
    blst_p1 p;
    blst_fp a, b, result;
    const blst_p1 *gen = blst_p1_generator();
    
    blst_p1_add(&p, gen, gen);  // 2*G
    a = p.x;
    b = p.y;
    
    for (auto _ : state) {
        blst_fp_mul(&result, &a, &b);
        benchmark::DoNotOptimize(result);
    }
}

static void BM_Fp_Sqr(benchmark::State& state) {
    blst_p1 p;
    blst_fp a, result;
    const blst_p1 *gen = blst_p1_generator();
    
    blst_p1_add(&p, gen, gen);  // 2*G
    a = p.x;
    
    for (auto _ : state) {
        blst_fp_sqr(&result, &a);
        benchmark::DoNotOptimize(result);
    }
}

// ============================================================================
// Fp2 Operations
// ============================================================================

static void BM_Fp2_Mul(benchmark::State& state) {
    blst_p2 p;
    blst_fp2 a, b, result;
    const blst_p2 *gen = blst_p2_generator();
    
    blst_p2_add(&p, gen, gen);  // 2*G
    a = p.x;
    b = p.y;
    
    for (auto _ : state) {
        blst_fp2_mul(&result, &a, &b);
        benchmark::DoNotOptimize(result);
    }
}

static void BM_Fp2_Sqr(benchmark::State& state) {
    blst_p2 p;
    blst_fp2 a, result;
    const blst_p2 *gen = blst_p2_generator();
    
    blst_p2_add(&p, gen, gen);  // 2*G
    a = p.x;
    
    for (auto _ : state) {
        blst_fp2_sqr(&result, &a);
        benchmark::DoNotOptimize(result);
    }
}

// ============================================================================
// Fp12 Operations
// ============================================================================

static void BM_Fp12_Sqr(benchmark::State& state) {
    blst_p1_affine p1_affine;
    blst_p2_affine p2_affine;
    blst_fp12 miller_result, result;
    const blst_p1 *gen1 = blst_p1_generator();
    const blst_p2 *gen2 = blst_p2_generator();
    
    blst_p1_to_affine(&p1_affine, gen1);
    blst_p2_to_affine(&p2_affine, gen2);
    blst_miller_loop(&miller_result, &p2_affine, &p1_affine);
    
    for (auto _ : state) {
        blst_fp12_sqr(&result, &miller_result);
        benchmark::DoNotOptimize(result);
    }
}

static void BM_Fp12_CyclotomicSqr(benchmark::State& state) {
    blst_p1_affine p1_affine;
    blst_p2_affine p2_affine;
    blst_fp12 miller_result, result;
    const blst_p1 *gen1 = blst_p1_generator();
    const blst_p2 *gen2 = blst_p2_generator();
    
    blst_p1_to_affine(&p1_affine, gen1);
    blst_p2_to_affine(&p2_affine, gen2);
    blst_miller_loop(&miller_result, &p2_affine, &p1_affine);
    
    for (auto _ : state) {
        blst_fp12_cyclotomic_sqr(&result, &miller_result);
        benchmark::DoNotOptimize(result);
    }
}

static void BM_Fp12_Mul(benchmark::State& state) {
    blst_p1_affine p1_affine;
    blst_p2_affine p2_affine;
    blst_fp12 a, b, result;
    const blst_p1 *gen1 = blst_p1_generator();
    const blst_p2 *gen2 = blst_p2_generator();
    
    blst_p1_to_affine(&p1_affine, gen1);
    blst_p2_to_affine(&p2_affine, gen2);
    blst_miller_loop(&a, &p2_affine, &p1_affine);
    blst_fp12_sqr(&b, &a);
    
    for (auto _ : state) {
        blst_fp12_mul(&result, &a, &b);
        benchmark::DoNotOptimize(result);
    }
}

static void BM_Fp12_MulByXy00z0(benchmark::State& state) {
    blst_p1_affine p1_affine;
    blst_p2_affine p2_affine;
    blst_fp12 a, result;
    blst_fp6 xy00z0;
    const blst_p1 *gen1 = blst_p1_generator();
    const blst_p2 *gen2 = blst_p2_generator();
    
    blst_p1_to_affine(&p1_affine, gen1);
    blst_p2_to_affine(&p2_affine, gen2);
    blst_miller_loop(&a, &p2_affine, &p1_affine);
    
    xy00z0.fp2[0] = p2_affine.x;
    xy00z0.fp2[1] = p2_affine.y;
    xy00z0.fp2[2] = p2_affine.x;
    
    for (auto _ : state) {
        blst_fp12_mul_by_xy00z0(&result, &a, &xy00z0);
        benchmark::DoNotOptimize(result);
    }
}

static void BM_Fp12_Conjugate(benchmark::State& state) {
    blst_p1_affine p1_affine;
    blst_p2_affine p2_affine;
    blst_fp12 a;
    const blst_p1 *gen1 = blst_p1_generator();
    const blst_p2 *gen2 = blst_p2_generator();
    
    blst_p1_to_affine(&p1_affine, gen1);
    blst_p2_to_affine(&p2_affine, gen2);
    blst_miller_loop(&a, &p2_affine, &p1_affine);
    
    for (auto _ : state) {
        blst_fp12_conjugate(&a);
        benchmark::DoNotOptimize(a);
    }
}

static void BM_Fp12_Inverse(benchmark::State& state) {
    blst_p1_affine p1_affine;
    blst_p2_affine p2_affine;
    blst_fp12 miller_result, result;
    const blst_p1 *gen1 = blst_p1_generator();
    const blst_p2 *gen2 = blst_p2_generator();
    
    blst_p1_to_affine(&p1_affine, gen1);
    blst_p2_to_affine(&p2_affine, gen2);
    blst_miller_loop(&miller_result, &p2_affine, &p1_affine);
    
    for (auto _ : state) {
        blst_fp12_inverse(&result, &miller_result);
        benchmark::DoNotOptimize(result);
    }
}

static void BM_Fp12_FrobeniusMap(benchmark::State& state) {
    blst_p1_affine p1_affine;
    blst_p2_affine p2_affine;
    blst_fp12 miller_result, result;
    const blst_p1 *gen1 = blst_p1_generator();
    const blst_p2 *gen2 = blst_p2_generator();
    
    blst_p1_to_affine(&p1_affine, gen1);
    blst_p2_to_affine(&p2_affine, gen2);
    blst_miller_loop(&miller_result, &p2_affine, &p1_affine);
    
    for (auto _ : state) {
        blst_fp12_frobenius_map(&result, &miller_result, 1);
        benchmark::DoNotOptimize(result);
    }
}

// ============================================================================
// Line Operations
// ============================================================================

static void BM_LineAdd(benchmark::State& state) {
    blst_p2 T, R;
    blst_p2_affine Q;
    blst_fp6 line;
    const blst_p2 *gen = blst_p2_generator();
    
    blst_p2_add(&R, gen, gen);  // 2*G
    blst_p2_add(&T, gen, &R);   // 3*G
    blst_p2_to_affine(&Q, gen);
    
    for (auto _ : state) {
        blst_line_add(&line, &T, &R, &Q);
        benchmark::DoNotOptimize(line);
    }
}

static void BM_LineDbl(benchmark::State& state) {
    blst_p2 T, Q;
    blst_fp6 line;
    const blst_p2 *gen = blst_p2_generator();
    
    blst_p2_add(&Q, gen, gen);  // 2*G
    blst_p2_add(&T, gen, &Q);   // 3*G
    
    for (auto _ : state) {
        blst_line_dbl(&line, &T, &Q);
        benchmark::DoNotOptimize(line);
    }
}

// Helper function to multiply line by Px2 coordinates
static void line_by_Px2(blst_fp6 *line, const blst_p1_affine *Px2) {
    blst_fp_mul(&line->fp2[1].fp[0], &line->fp2[1].fp[0], &Px2->x);
    blst_fp_mul(&line->fp2[1].fp[1], &line->fp2[1].fp[1], &Px2->x);
    blst_fp_mul(&line->fp2[2].fp[0], &line->fp2[2].fp[0], &Px2->y);
    blst_fp_mul(&line->fp2[2].fp[1], &line->fp2[2].fp[1], &Px2->y);
}

static void BM_LineAddByPx2(benchmark::State& state) {
    blst_p2 T, R;
    blst_p2_affine Q;
    blst_p1_affine P, Px2;
    blst_fp6 line;
    const blst_p2 *gen2 = blst_p2_generator();
    const blst_p1 *gen1 = blst_p1_generator();
    
    blst_p2_add(&R, gen2, gen2);  // 2*G
    blst_p2_add(&T, gen2, &R);    // 3*G
    blst_p2_to_affine(&Q, gen2);
    blst_p1_to_affine(&P, gen1);
    
    blst_fp tmp;
    blst_fp_add(&tmp, &P.x, &P.x);
    blst_fp_cneg(&Px2.x, &tmp, 1);
    blst_fp_add(&Px2.y, &P.y, &P.y);
    
    for (auto _ : state) {
        blst_line_add(&line, &T, &R, &Q);
        line_by_Px2(&line, &Px2);
        benchmark::DoNotOptimize(line);
    }
}

static void BM_LineDblByPx2(benchmark::State& state) {
    blst_p2 T, Q;
    blst_p1_affine P, Px2;
    blst_fp6 line;
    const blst_p2 *gen2 = blst_p2_generator();
    const blst_p1 *gen1 = blst_p1_generator();
    
    blst_p2_add(&Q, gen2, gen2);  // 2*G
    blst_p2_add(&T, gen2, &Q);    // 3*G
    blst_p1_to_affine(&P, gen1);
    
    blst_fp tmp;
    blst_fp_add(&tmp, &P.x, &P.x);
    blst_fp_cneg(&Px2.x, &tmp, 1);
    blst_fp_add(&Px2.y, &P.y, &P.y);
    
    for (auto _ : state) {
        blst_line_dbl(&line, &T, &Q);
        line_by_Px2(&line, &Px2);
        benchmark::DoNotOptimize(line);
    }
}

// ============================================================================
// Register Benchmarks
// ============================================================================

// G1 Operations
BENCHMARK(BM_G1_Add);
BENCHMARK(BM_G1_Double);
BENCHMARK(BM_G1_MixedAdd);
BENCHMARK(BM_G1_ScalarMul);

// G2 Operations
BENCHMARK(BM_G2_Add);
BENCHMARK(BM_G2_Double);
BENCHMARK(BM_G2_MixedAdd);
BENCHMARK(BM_G2_ScalarMul);

// Pairing Operations
BENCHMARK(BM_MillerLoop);
BENCHMARK(BM_FinalExp);
BENCHMARK(BM_FullPairing);

// Fp Operations
BENCHMARK(BM_Fp_Mul);
BENCHMARK(BM_Fp_Sqr);

// Fp2 Operations
BENCHMARK(BM_Fp2_Mul);
BENCHMARK(BM_Fp2_Sqr);

// Fp12 Operations
BENCHMARK(BM_Fp12_Sqr);
BENCHMARK(BM_Fp12_CyclotomicSqr);
BENCHMARK(BM_Fp12_Mul);
BENCHMARK(BM_Fp12_MulByXy00z0);
BENCHMARK(BM_Fp12_Conjugate);
BENCHMARK(BM_Fp12_Inverse);
BENCHMARK(BM_Fp12_FrobeniusMap);

// Line Operations
BENCHMARK(BM_LineAdd);
BENCHMARK(BM_LineDbl);
BENCHMARK(BM_LineAddByPx2);
BENCHMARK(BM_LineDblByPx2);

BENCHMARK_MAIN();

