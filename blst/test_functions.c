/*
 * Test program for blst library operations (src_small version)
 * Uses fixed seeded random inputs to ensure deterministic output
 */

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "point.h"
#include "fields.h"

// External declarations
extern const POINTonE1 BLS12_381_G1;
extern const POINTonE2 BLS12_381_G2;

// Function declarations
void line_add(vec384fp6 line, POINTonE2 *T, const POINTonE2 *R,
              const POINTonE2_affine *Q);
void line_dbl(vec384fp6 line, POINTonE2 *T, const POINTonE2 *Q);
void line_by_Px2(vec384fp6 line, const POINTonE1_affine *Px2);
void blst_miller_loop(vec384fp12 ret, const POINTonE2_affine *Q,
                      const POINTonE1_affine *P);
void blst_final_exp(vec384fp12 ret, const vec384fp12 f);

void blst_p1_add(POINTonE1 *out, const POINTonE1 *a, const POINTonE1 *b);
void blst_p1_double(POINTonE1 *out, const POINTonE1 *a);
void blst_p1_mult(POINTonE1 *out, const POINTonE1 *a, const byte *scalar, size_t nbits);
void blst_p1_to_affine(POINTonE1_affine *out, const POINTonE1 *a);

void blst_p2_add(POINTonE2 *out, const POINTonE2 *a, const POINTonE2 *b);
void blst_p2_double(POINTonE2 *out, const POINTonE2 *a);
void blst_p2_mult(POINTonE2 *out, const POINTonE2 *a, const byte *scalar, size_t nbits);
void blst_p2_to_affine(POINTonE2_affine *out, const POINTonE2 *a);

void blst_fp12_sqr(vec384fp12 ret, const vec384fp12 a);
void blst_fp12_cyclotomic_sqr(vec384fp12 ret, const vec384fp12 a);
void blst_fp12_mul(vec384fp12 ret, const vec384fp12 a, const vec384fp12 b);
void blst_fp12_inverse(vec384fp12 ret, const vec384fp12 a);
void blst_fp12_frobenius_map(vec384fp12 ret, const vec384fp12 a, size_t n);

// Simple seeded PRNG (LCG)
static uint64_t seed = 12345;
static uint32_t lcg_next(void) {
    seed = seed * 1103515245 + 12345;
    return (uint32_t)(seed >> 32);
}

static void fill_random_bytes(byte *buf, size_t len) {
    for (size_t i = 0; i < len; i++) {
        buf[i] = (byte)(lcg_next() & 0xff);
    }
}

// Helper to print vec384 in hex
static void print_vec384(const char *label, const vec384 v) {
    printf("%s: ", label);
    for (int i = 5; i >= 0; i--) {
        printf("%016lx", v[i]);
    }
    printf("\n");
}

// Helper to print vec384fp12 in hex
static void print_vec384fp12(const char *label, const vec384fp12 v) {
    printf("%s: ", label);
    for (int i = 1; i >= 0; i--) {
        for (int j = 2; j >= 0; j--) {
            for (int k = 1; k >= 0; k--) {
                for (int l = 5; l >= 0; l--) {
                    printf("%016lx", v[i][j][k][l]);
                }
            }
        }
    }
    printf("\n");
}

static void test_g1_add(void) {
    POINTonE1 p1, p2, result;
    const POINTonE1 *gen = &BLS12_381_G1;
    
    blst_p1_add(&p1, gen, gen);
    blst_p1_add(&p2, gen, &p1);
    blst_p1_add(&result, &p1, &p2);
    
    printf("=== G1 Addition ===\n");
    print_vec384("result.X", result.X);
    print_vec384("result.Y", result.Y);
    print_vec384("result.Z", result.Z);
}

static void test_g1_double(void) {
    POINTonE1 p1, result;
    const POINTonE1 *gen = &BLS12_381_G1;
    
    blst_p1_add(&p1, gen, gen);
    blst_p1_double(&result, &p1);
    
    printf("=== G1 Double ===\n");
    print_vec384("result.X", result.X);
    print_vec384("result.Y", result.Y);
    print_vec384("result.Z", result.Z);
}

static void test_g1_scalar_mul(void) {
    POINTonE1 p1, result;
    const POINTonE1 *gen = &BLS12_381_G1;
    byte scalar[32];
    
    seed = 12345;
    fill_random_bytes(scalar, sizeof(scalar));
    scalar[0] |= 0x01; // Ensure non-zero
    
    blst_p1_add(&p1, gen, gen);
    blst_p1_mult(&result, &p1, scalar, 256);
    
    printf("=== G1 Scalar Mul ===\n");
    print_vec384("result.X", result.X);
    print_vec384("result.Y", result.Y);
    print_vec384("result.Z", result.Z);
}

static void test_g2_add(void) {
    POINTonE2 p1, p2, result;
    const POINTonE2 *gen = &BLS12_381_G2;
    
    blst_p2_add(&p1, gen, gen);
    blst_p2_add(&p2, gen, &p1);
    blst_p2_add(&result, &p1, &p2);
    
    printf("=== G2 Addition ===\n");
    print_vec384("result.X[0]", result.X[0]);
    print_vec384("result.X[1]", result.X[1]);
    print_vec384("result.Y[0]", result.Y[0]);
    print_vec384("result.Y[1]", result.Y[1]);
    print_vec384("result.Z[0]", result.Z[0]);
    print_vec384("result.Z[1]", result.Z[1]);
}

static void test_g2_double(void) {
    POINTonE2 p1, result;
    const POINTonE2 *gen = &BLS12_381_G2;
    
    blst_p2_add(&p1, gen, gen);
    blst_p2_double(&result, &p1);
    
    printf("=== G2 Double ===\n");
    print_vec384("result.X[0]", result.X[0]);
    print_vec384("result.X[1]", result.X[1]);
    print_vec384("result.Y[0]", result.Y[0]);
    print_vec384("result.Y[1]", result.Y[1]);
    print_vec384("result.Z[0]", result.Z[0]);
    print_vec384("result.Z[1]", result.Z[1]);
}

static void test_g2_scalar_mul(void) {
    POINTonE2 p1, result;
    const POINTonE2 *gen = &BLS12_381_G2;
    byte scalar[32];
    
    seed = 12345;
    fill_random_bytes(scalar, sizeof(scalar));
    scalar[0] |= 0x01; // Ensure non-zero
    
    blst_p2_add(&p1, gen, gen);
    blst_p2_mult(&result, &p1, scalar, 256);
    
    printf("=== G2 Scalar Mul ===\n");
    print_vec384("result.X[0]", result.X[0]);
    print_vec384("result.X[1]", result.X[1]);
    print_vec384("result.Y[0]", result.Y[0]);
    print_vec384("result.Y[1]", result.Y[1]);
    print_vec384("result.Z[0]", result.Z[0]);
    print_vec384("result.Z[1]", result.Z[1]);
}

static void test_fp12_sqr(void) {
    POINTonE1_affine p1_affine;
    POINTonE2_affine p2_affine;
    vec384fp12 miller_result, result;
    const POINTonE1 *gen1 = &BLS12_381_G1;
    const POINTonE2 *gen2 = &BLS12_381_G2;
    
    blst_p1_to_affine(&p1_affine, gen1);
    blst_p2_to_affine(&p2_affine, gen2);
    blst_miller_loop(miller_result, &p2_affine, &p1_affine);
    blst_fp12_sqr(result, miller_result);
    
    printf("=== Fp12 Square ===\n");
    print_vec384fp12("result", result);
}

static void test_fp12_mul(void) {
    POINTonE1_affine p1_affine;
    POINTonE2_affine p2_affine;
    vec384fp12 a, b, result;
    const POINTonE1 *gen1 = &BLS12_381_G1;
    const POINTonE2 *gen2 = &BLS12_381_G2;
    
    blst_p1_to_affine(&p1_affine, gen1);
    blst_p2_to_affine(&p2_affine, gen2);
    blst_miller_loop(a, &p2_affine, &p1_affine);
    blst_fp12_sqr(b, a);
    blst_fp12_mul(result, a, b);
    
    printf("=== Fp12 Multiply ===\n");
    print_vec384fp12("result", result);
}

static void test_fp12_frobenius(void) {
    POINTonE1_affine p1_affine;
    POINTonE2_affine p2_affine;
    vec384fp12 miller_result, result;
    const POINTonE1 *gen1 = &BLS12_381_G1;
    const POINTonE2 *gen2 = &BLS12_381_G2;
    
    blst_p1_to_affine(&p1_affine, gen1);
    blst_p2_to_affine(&p2_affine, gen2);
    blst_miller_loop(miller_result, &p2_affine, &p1_affine);
    blst_fp12_frobenius_map(result, miller_result, 1);
    
    printf("=== Fp12 Frobenius ===\n");
    print_vec384fp12("result", result);
}

static void test_fp12_inverse(void) {
    POINTonE1_affine p1_affine;
    POINTonE2_affine p2_affine;
    vec384fp12 miller_result, result;
    const POINTonE1 *gen1 = &BLS12_381_G1;
    const POINTonE2 *gen2 = &BLS12_381_G2;
    
    blst_p1_to_affine(&p1_affine, gen1);
    blst_p2_to_affine(&p2_affine, gen2);
    blst_miller_loop(miller_result, &p2_affine, &p1_affine);
    blst_fp12_inverse(result, miller_result);
    
    printf("=== Fp12 Inverse ===\n");
    print_vec384fp12("result", result);
}

static void test_cyclotomic_sqr_fp12(void) {
    POINTonE1_affine p1_affine;
    POINTonE2_affine p2_affine;
    vec384fp12 miller_result, result;
    const POINTonE1 *gen1 = &BLS12_381_G1;
    const POINTonE2 *gen2 = &BLS12_381_G2;
    
    blst_p1_to_affine(&p1_affine, gen1);
    blst_p2_to_affine(&p2_affine, gen2);
    blst_miller_loop(miller_result, &p2_affine, &p1_affine);
    blst_fp12_cyclotomic_sqr(result, miller_result);
    
    printf("=== Cyclotomic Sqr ===\n");
    print_vec384fp12("result", result);
}

static void test_line_add(void) {
    POINTonE2 T, R;
    POINTonE2_affine Q;
    vec384fp6 line;
    const POINTonE2 *gen = &BLS12_381_G2;
    
    blst_p2_add(&R, gen, gen);
    blst_p2_add(&T, gen, &R);
    blst_p2_to_affine(&Q, gen);
    line_add(line, &T, &R, &Q);
    
    printf("=== Line Add ===\n");
    print_vec384("line[0][0]", line[0][0]);
    print_vec384("line[0][1]", line[0][1]);
    print_vec384("line[1][0]", line[1][0]);
    print_vec384("line[1][1]", line[1][1]);
    print_vec384("line[2][0]", line[2][0]);
    print_vec384("line[2][1]", line[2][1]);
}

static void test_line_dbl(void) {
    POINTonE2 T, Q;
    vec384fp6 line;
    const POINTonE2 *gen = &BLS12_381_G2;
    
    blst_p2_add(&Q, gen, gen);
    blst_p2_add(&T, gen, &Q);
    line_dbl(line, &T, &Q);
    
    printf("=== Line Dbl ===\n");
    print_vec384("line[0][0]", line[0][0]);
    print_vec384("line[0][1]", line[0][1]);
    print_vec384("line[1][0]", line[1][0]);
    print_vec384("line[1][1]", line[1][1]);
    print_vec384("line[2][0]", line[2][0]);
    print_vec384("line[2][1]", line[2][1]);
}

static void test_line_by_Px2(void) {
    POINTonE2 T, R;
    POINTonE2_affine Q;
    POINTonE1_affine Px2;
    vec384fp6 line;
    const POINTonE2 *gen2 = &BLS12_381_G2;
    const POINTonE1 *gen1 = &BLS12_381_G1;
    
    blst_p2_add(&R, gen2, gen2);
    blst_p2_add(&T, gen2, &R);
    blst_p2_to_affine(&Q, gen2);
    blst_p1_to_affine(&Px2, gen1);
    line_add(line, &T, &R, &Q);
    line_by_Px2(line, &Px2);
    
    printf("=== Line by Px2 ===\n");
    print_vec384("line[0][0]", line[0][0]);
    print_vec384("line[0][1]", line[0][1]);
    print_vec384("line[1][0]", line[1][0]);
    print_vec384("line[1][1]", line[1][1]);
    print_vec384("line[2][0]", line[2][0]);
    print_vec384("line[2][1]", line[2][1]);
}

static void test_miller_loop(void) {
    POINTonE1_affine p1_affine;
    POINTonE2_affine p2_affine;
    vec384fp12 result;
    const POINTonE1 *gen1 = &BLS12_381_G1;
    const POINTonE2 *gen2 = &BLS12_381_G2;
    
    blst_p1_to_affine(&p1_affine, gen1);
    blst_p2_to_affine(&p2_affine, gen2);
    blst_miller_loop(result, &p2_affine, &p1_affine);
    
    printf("=== Miller Loop ===\n");
    print_vec384fp12("result", result);
}

static void test_final_exp(void) {
    POINTonE1_affine p1_affine;
    POINTonE2_affine p2_affine;
    vec384fp12 miller_result, result;
    const POINTonE1 *gen1 = &BLS12_381_G1;
    const POINTonE2 *gen2 = &BLS12_381_G2;
    
    blst_p1_to_affine(&p1_affine, gen1);
    blst_p2_to_affine(&p2_affine, gen2);
    blst_miller_loop(miller_result, &p2_affine, &p1_affine);
    blst_final_exp(result, miller_result);
    
    printf("=== Final Exp ===\n");
    print_vec384fp12("result", result);
}

int main(void) {
    printf("BLST Function Tests (src_small)\n");
    printf("===============================\n\n");
    
    test_g1_add();
    test_g1_double();
    test_g1_scalar_mul();
    test_g2_add();
    test_g2_double();
    test_g2_scalar_mul();
    test_fp12_sqr();
    test_fp12_mul();
    test_fp12_frobenius();
    test_fp12_inverse();
    test_cyclotomic_sqr_fp12();
    test_line_add();
    test_line_dbl();
    test_line_by_Px2();
    test_miller_loop();
    test_final_exp();
    
    return 0;
}

#ifdef __cplusplus
}
#endif














