#include <iostream>
#include <cstring>
#include <cassert>

extern "C" {
#include "fields.h"
#include "point.h"
void blst_fp12_inverse(vec384fp12 ret, const vec384fp12 a);
void mul_by_u_plus_1_fp2_export(vec384x ret, const vec384x a);
void blst_fp12_mul(vec384fp12 ret, const vec384fp12 a, const vec384fp12 b);
int blst_fp12_is_one(const vec384fp12 a);
const vec384fp12 *blst_fp12_one(void);
void blst_p1_to_affine(POINTonE1_affine *out, const POINTonE1 *a);
void blst_p2_to_affine(POINTonE2_affine *out, const POINTonE2 *a);
void blst_miller_loop(vec384fp12 ret, const POINTonE2_affine *Q, const POINTonE1_affine *P);
extern const POINTonE1 BLS12_381_G1;
extern const POINTonE2 BLS12_381_G2;
}

#include "inversion_parts.hpp"

// Helper to check if two fp12 elements are equal
static bool fp12_equal(const vec384fp12 a, const vec384fp12 b) {
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 2; k++) {
                if (memcmp(a[i][j][k], b[i][j][k], sizeof(vec384)) != 0) {
                    return false;
                }
            }
        }
    }
    return true;
}

// Helper to print fp (vec384) as hex
static void print_fp(const char *name, const vec384 x) {
    printf("%s = ", name);
    for (int j = 5; j >= 0; j--) {
        printf("%016lx", x[j]);
    }
    printf("\n");
}

// Helper to print fp2 (two fp elements)
static void print_fp2(const char *name, const vec384fp2 x) {
    printf("%s[0] = ", name);
    for (int j = 5; j >= 0; j--) printf("%016lx", x[0][j]);
    printf("\n%s[1] = ", name);
    for (int j = 5; j >= 0; j--) printf("%016lx", x[1][j]);
    printf("\n");
}

// Helper to print fp6 (three fp2 components)
static void print_fp6(const char *name, const vec384fp6 x) {
    for (int i = 0; i < 3; i++) {
        printf("%s[%d][0] = ", name, i);
        for (int j = 5; j >= 0; j--) printf("%016lx", x[i][0][j]);
        printf("\n%s[%d][1] = ", name, i);
        for (int j = 5; j >= 0; j--) printf("%016lx", x[i][1][j]);
        printf("\n");
    }
}

// Helper to print fp12 element (first element only for brevity)
static void print_fp12_element(const char *name, const vec384fp12 x) {
    printf("%s[0][0][0] = ", name);
    for (int j = 5; j >= 0; j--) {
        printf("%016lx", x[0][0][0][j]);
    }
    printf("\n");
}

// Print full fp12 (all 2*3*2 = 12 fp components)
static void print_fp12_full(const char *name, const vec384fp12 x) {
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 3; j++)
            for (int k = 0; k < 2; k++) {
                printf("%s[%d][%d][%d] = ", name, i, j, k);
                for (int l = 5; l >= 0; l--) printf("%016lx", x[i][j][k][l]);
                printf("\n");
            }
}

// Run one end-to-end inversion and print input/output of each part with explicit labels.
// Each part lists INPUT(s) and OUTPUT(s) by name and count so test_fp12_flat_inverse can align.
static void print_inversion_parts_trace(const vec384fp12 a) {
    printf("========== Inversion parts trace (one end-to-end input) ==========\n\n");

    printf("--- Part 0 (initial): INPUT a (fp12, 12 fp) ---\n");
    print_fp12_full("a", a);
    printf("\n");

    // Part 1: fp12 -> fp6. INPUT: a (12). OUTPUT: t0_fp6 (6).
    vec384fp6 t0_fp6;
    inverse_fp12_to_fp6_wrapper_export(t0_fp6, a);
    printf("--- Part 1: inverse_fp12_to_fp6 ---\n");
    printf("  INPUT:  a (fp12, 12 fp)\n");
    printf("  OUTPUT: t0_fp6 (fp6, 6 fp)\n");
    print_fp6("t0_fp6", t0_fp6);
    printf("\n");

    // Part 2: fp6 -> ret_fp2, c0, c1, c2. INPUT: t0_fp6 (6). OUTPUT: ret_fp2 (2), c0 (2), c1 (2), c2 (2).
    vec384fp2 ret_fp2, c0, c1, c2;
    inverse_fp6_to_fp2_wrapper_export(ret_fp2, c0, c1, c2, t0_fp6);
    printf("--- Part 2: inverse_fp6_to_fp2 ---\n");
    printf("  INPUT:  t0_fp6 (fp6, 6 fp)\n");
    printf("  OUTPUT: ret_fp2 (fp2, 2 fp), c0 (fp2, 2 fp), c1 (fp2, 2 fp), c2 (fp2, 2 fp)\n");
    print_fp2("ret_fp2", ret_fp2);
    print_fp2("c0", c0);
    print_fp2("c1", c1);
    print_fp2("c2", c2);
    printf("\n");

    // Part 3: fp2 -> fp (norm). INPUT: ret_fp2 (2). OUTPUT: norm_fp (1).
    vec384 norm_fp;
    inverse_fp2_to_fp_wrapper_export(norm_fp, ret_fp2);
    printf("--- Part 3: inverse_fp2_to_fp ---\n");
    printf("  INPUT:  ret_fp2 (fp2, 2 fp)\n");
    printf("  OUTPUT: norm_fp (fp, 1 fp) = a^2 + b^2\n");
    print_fp("norm_fp", norm_fp);
    printf("\n");

    // Part 4: fp inverse. INPUT: norm_fp (1). OUTPUT: r_fp (1).
    vec384 r_fp;
    blst_fp_inverse(r_fp, norm_fp);
    printf("--- Part 4: blst_fp_inverse ---\n");
    printf("  INPUT:  norm_fp (fp, 1 fp)\n");
    printf("  OUTPUT: r_fp (fp, 1 fp) = 1/norm_fp\n");
    print_fp("r_fp", r_fp);
    printf("\n");

    // Part 5: fp -> fp2. INPUT: r_fp (1), ret_fp2 (2). OUTPUT: t1_fp2 (2).
    vec384fp2 t1_fp2;
    inverse_fp_to_fp2_wrapper_export(t1_fp2, r_fp, ret_fp2);
    printf("--- Part 5: inverse_fp_to_fp2 ---\n");
    printf("  INPUT:  r_fp (fp, 1 fp), ret_fp2 (fp2, 2 fp)\n");
    printf("  OUTPUT: t1_fp2 (fp2, 2 fp)\n");
    print_fp2("t1_fp2", t1_fp2);
    printf("\n");

    // Part 6: fp2 -> fp6. INPUT: t1_fp2 (2), c0 (2), c1 (2), c2 (2). OUTPUT: t1_fp6 (6).
    vec384fp6 t1_fp6;
    inverse_fp2_to_fp6_wrapper_export(t1_fp6, t1_fp2, c0, c1, c2);
    printf("--- Part 6: inverse_fp2_to_fp6 ---\n");
    printf("  INPUT:  t1_fp2 (fp2, 2 fp), c0 (fp2, 2 fp), c1 (fp2, 2 fp), c2 (fp2, 2 fp)\n");
    printf("  OUTPUT: t1_fp6 (fp6, 6 fp)\n");
    print_fp6("t1_fp6", t1_fp6);
    printf("\n");

    // Part 7: fp6 -> fp12. INPUT: t1_fp6 (6), a (12). OUTPUT: result (12).
    vec384fp12 result_parts;
    inverse_fp6_to_fp12_wrapper_export(result_parts, t1_fp6, a);
    printf("--- Part 7: inverse_fp6_to_fp12 ---\n");
    printf("  INPUT:  t1_fp6 (fp6, 6 fp), a (fp12, 12 fp)\n");
    printf("  OUTPUT: result (fp12, 12 fp) = inverse of a\n");
    print_fp12_full("result", result_parts);
    printf("\n========== End of trace ==========\n");
}

// Test function that chains all inversion parts together
static bool test_inversion_parts_chain(const vec384fp12 a) {
    vec384fp12 result_parts, result_blst;
    
    // Step 1: fp12 -> fp6
    vec384fp6 t0_fp6;
    inverse_fp12_to_fp6_wrapper_export(t0_fp6, a);
    
    // Step 2: fp6 -> fp2 (computes ret, c0, c1, c2)
    vec384fp2 ret_fp2, c0, c1, c2;
    inverse_fp6_to_fp2_wrapper_export(ret_fp2, c0, c1, c2, t0_fp6);
    
    // Step 3: fp2 -> fp (compute norm = a^2 + b^2)
    vec384 norm_fp;
    inverse_fp2_to_fp_wrapper_export(norm_fp, ret_fp2);
    
    // Step 4: fp -> fp (compute reciprocal of norm)
    vec384 r_fp;
    blst_fp_inverse(r_fp, norm_fp);
    
    // Step 5: fp -> fp2 (multiply ret_fp2 by reciprocal)
    vec384fp2 t1_fp2;
    inverse_fp_to_fp2_wrapper_export(t1_fp2, r_fp, ret_fp2);
    
    // Step 6: fp2 -> fp6 (multiply c0, c1, c2 by t1_fp2)
    vec384fp6 t1_fp6;
    inverse_fp2_to_fp6_wrapper_export(t1_fp6, t1_fp2, c0, c1, c2);
    
    // Step 7: fp6 -> fp12 (multiply a[0] and a[1] by t1_fp6, negate a[1])
    inverse_fp6_to_fp12_wrapper_export(result_parts, t1_fp6, a);
    
    // Compare with blst_fp12_inverse
    blst_fp12_inverse(result_blst, a);
    
    bool match = fp12_equal(result_parts, result_blst);
    
    if (!match) {
        printf("MISMATCH detected!\n");
        print_fp12_element("result_parts", result_parts);
        print_fp12_element("result_blst", result_blst);
    }
    
    return match;
}

/*
 * How inverse_fp6_to_fp2 outputs 0 and 1 (ret_fp2) are computed in BLST:
 *
 * Layout: vec384fp6 a has a[0]=(x0,x1), a[1]=(x2,x3), a[2]=(x4,x5) (each fp2: [0]=real, [1]=imag).
 * Wrapper (empirical_basis): flat a[i] -> a_fp6[i/2][i%2], so a_fp6[0][0]=x0, a_fp6[0][1]=x1, etc.
 *
 * Steps in inverse_fp6_to_fp2:
 *   c0 = a[0]^2 - (a[1]*a[2])*(u+1)
 *   c1 = a[2]^2*(u+1) - (a[0]*a[1])
 *   c2 = a[1]^2 - a[0]*a[2]
 *   ret = (a[2]*c1 + a[1]*c2)*(u+1) + a[0]*c0
 *
 * In Fp2, (u+1)*(a,b) = (a-b, a+b) (from fp12_tower.c mul_by_u_plus_1_fp2x2).
 * So ret[0] = (a2*c1 + a1*c2)_real - (a2*c1 + a1*c2)_imag + (a0*c0)_real
 *    ret[1] = (a2*c1 + a1*c2)_real + (a2*c1 + a1*c2)_imag + (a0*c0)_imag
 *
 * The empirical basis fits polynomials using inverse_fp6_to_fp2_wrapper, which returns
 * from_mont(ret[0]) and from_mont(ret[1]) as int64_t. from_mont() only uses the low limb,
 * so outputs 0 and 1 are truncated to 64 bits when fitting. The fitted Poly 0 and Poly 1
 * can therefore be wrong for full 384-bit values. Outputs 2-7 (c0, c1, c2) often stay
 * small for the probe inputs 0,1,-1 so those formulas match.
 *
 * Correct way to get output 0 and 1: use the BLST sequence above (compute c0,c1,c2 then
 * ret = (a2*c1 + a1*c2)*(u+1) + a0*c0), or run a 384-bit-aware basis generator.
 */
// Compute ret (outputs 0 and 1) exactly as BLST does: c0,c1,c2 then ret = (a2*c1 + a1*c2)*(u+1) + a0*c0.
static void compute_ret_fp2_from_blst_steps(vec384fp2 ret, const vec384fp2 c0, const vec384fp2 c1, const vec384fp2 c2, const vec384fp6 a) {
    vec384x t0, t1;
    mul_fp2(t0, c1, a[2]);   // a[2]*c1
    mul_fp2(t1, c2, a[1]);   // a[1]*c2
    add_fp2(t0, t0, t1);     // a[2]*c1 + a[1]*c2
    mul_by_u_plus_1_fp2_export(t0, t0);  // (u+1)*(...)
    mul_fp2(t1, c0, a[0]);   // a[0]*c0
    add_fp2(t0, t0, t1);     // ret
    vec_copy(ret[0], t0[0], sizeof(t0[0]));
    vec_copy(ret[1], t0[1], sizeof(t0[1]));
}

// Evaluate the stated fp6->fp2 polynomial formulas (in Fp) and compare to inverse_fp6_to_fp2 output.
// Flat layout: x0=a[0][0], x1=a[0][1], x2=a[1][0], x3=a[1][1], x4=a[2][0], x5=a[2][1].
// Output order: ret[0], ret[1], c0[0], c0[1], c1[0], c1[1], c2[0], c2[1] = poly0..poly7.
static bool test_fp6_to_fp2_polynomial_formulas(const vec384fp6 a) {
    const vec384 &x0 = a[0][0], &x1 = a[0][1], &x2 = a[1][0], &x3 = a[1][1], &x4 = a[2][0], &x5 = a[2][1];
    vec384 t, t2, t6;
    vec384 expected[8];

    // Poly 0: 1*x0 + 1*x2 + 1*x3 + 2*x5 - 3*x0*x1 - 6*x2*x3 - 6*x4*x5
    mul_fp(t, x0, x1);      mul_by_3_fp(t, t);   // 3*x0*x1
    mul_fp(t2, x2, x3);    mul_by_3_fp(t6, t2); add_fp(t6, t6, t6);  // 6*x2*x3
    mul_fp(t2, x4, x5);    mul_by_3_fp(t2, t2); add_fp(t2, t2, t2);  // 6*x4*x5
    add_fp(expected[0], x0, x2);
    add_fp(expected[0], expected[0], x3);
    add_fp(expected[0], expected[0], x5);
    add_fp(expected[0], expected[0], x5);       // + 2*x5
    sub_fp(expected[0], expected[0], t);
    sub_fp(expected[0], expected[0], t6);
    sub_fp(expected[0], expected[0], t2);

    // Poly 1: -1*x1 + 1*x2 - 1*x3 + 2*x4 + 3*x0*x1 - 6*x4*x5
    mul_fp(t, x0, x1);     mul_by_3_fp(t, t);   // 3*x0*x1
    mul_fp(t2, x4, x5);    mul_by_3_fp(t6, t2); add_fp(t6, t6, t6);  // 6*x4*x5
    neg_fp(expected[1], x1);
    add_fp(expected[1], expected[1], x2);
    sub_fp(expected[1], expected[1], x3);
    add_fp(t2, x4, x4);    add_fp(expected[1], expected[1], t2);
    add_fp(expected[1], expected[1], t);
    sub_fp(expected[1], expected[1], t6);

    // Poly 2: 1*x0*x0 - 1*x1*x1 - 1*x2*x4 + 1*x2*x5 + 1*x3*x4 + 1*x3*x5
    sqr_fp(expected[2], x0);
    sqr_fp(t, x1);         sub_fp(expected[2], expected[2], t);
    mul_fp(t, x2, x4);     sub_fp(expected[2], expected[2], t);
    mul_fp(t, x2, x5);     add_fp(expected[2], expected[2], t);
    mul_fp(t, x3, x4);     add_fp(expected[2], expected[2], t);
    mul_fp(t, x3, x5);     add_fp(expected[2], expected[2], t);

    // Poly 3: 2*x0*x1 - 1*x2*x4 - 1*x2*x5 - 1*x3*x4 + 1*x3*x5
    mul_fp(expected[3], x0, x1);  add_fp(expected[3], expected[3], expected[3]);
    mul_fp(t, x2, x4);     sub_fp(expected[3], expected[3], t);
    mul_fp(t, x2, x5);     sub_fp(expected[3], expected[3], t);
    mul_fp(t, x3, x4);     sub_fp(expected[3], expected[3], t);
    mul_fp(t, x3, x5);     add_fp(expected[3], expected[3], t);

    // Poly 4: 1*x4*x4 - 1*x5*x5 - 1*x0*x2 + 1*x1*x3 - 2*x4*x5
    sqr_fp(expected[4], x4);
    sqr_fp(t, x5);         sub_fp(expected[4], expected[4], t);
    mul_fp(t, x0, x2);     sub_fp(expected[4], expected[4], t);
    mul_fp(t, x1, x3);     add_fp(expected[4], expected[4], t);
    mul_fp(t, x4, x5);     add_fp(t, t, t);      sub_fp(expected[4], expected[4], t);

    // Poly 5: 1*x4*x4 - 1*x5*x5 - 1*x0*x3 - 1*x1*x2 + 2*x4*x5
    sqr_fp(expected[5], x4);
    sqr_fp(t, x5);         sub_fp(expected[5], expected[5], t);
    mul_fp(t, x0, x3);     sub_fp(expected[5], expected[5], t);
    mul_fp(t, x1, x2);     sub_fp(expected[5], expected[5], t);
    mul_fp(t, x4, x5);     add_fp(t, t, t);      add_fp(expected[5], expected[5], t);

    // Poly 6: 1*x2*x2 - 1*x3*x3 - 1*x0*x4 + 1*x1*x5
    sqr_fp(expected[6], x2);
    sqr_fp(t, x3);         sub_fp(expected[6], expected[6], t);
    mul_fp(t, x0, x4);     sub_fp(expected[6], expected[6], t);
    mul_fp(t, x1, x5);     add_fp(expected[6], expected[6], t);

    // Poly 7: -1*x0*x5 - 1*x1*x4 + 2*x2*x3
    mul_fp(expected[7], x0, x5);  neg_fp(expected[7], expected[7]);
    mul_fp(t, x1, x4);     sub_fp(expected[7], expected[7], t);
    mul_fp(t, x2, x3);     add_fp(t, t, t);      add_fp(expected[7], expected[7], t);

    vec384fp2 ret, c0, c1, c2;
    inverse_fp6_to_fp2_wrapper_export(ret, c0, c1, c2, a);

    bool all_ok = true;
    if (memcmp(expected[0], ret[0], sizeof(vec384)) != 0) { printf("fp6_to_fp2 formula mismatch at output 0 (poly 0)\n"); all_ok = false; }
    if (memcmp(expected[1], ret[1], sizeof(vec384)) != 0) { printf("fp6_to_fp2 formula mismatch at output 1 (poly 1)\n"); all_ok = false; }
    if (memcmp(expected[2], c0[0], sizeof(vec384)) != 0) { printf("fp6_to_fp2 formula mismatch at output 2 (poly 2)\n"); all_ok = false; }
    if (memcmp(expected[3], c0[1], sizeof(vec384)) != 0) { printf("fp6_to_fp2 formula mismatch at output 3 (poly 3)\n"); all_ok = false; }
    if (memcmp(expected[4], c1[0], sizeof(vec384)) != 0) { printf("fp6_to_fp2 formula mismatch at output 4 (poly 4)\n"); all_ok = false; }
    if (memcmp(expected[5], c1[1], sizeof(vec384)) != 0) { printf("fp6_to_fp2 formula mismatch at output 5 (poly 5)\n"); all_ok = false; }
    if (memcmp(expected[6], c2[0], sizeof(vec384)) != 0) { printf("fp6_to_fp2 formula mismatch at output 6 (poly 6)\n"); all_ok = false; }
    if (memcmp(expected[7], c2[1], sizeof(vec384)) != 0) { printf("fp6_to_fp2 formula mismatch at output 7 (poly 7)\n"); all_ok = false; }
    return all_ok;
}

// Test with multiple random inputs
static bool test_inversion_parts_multiple() {
    printf("Testing inversion parts chain against blst_fp12_inverse...\n\n");
    
    // Test 1: Use a known fp12 element (from miller loop)
    POINTonE1_affine p1_affine;
    POINTonE2_affine p2_affine;
    vec384fp12 test_input;
    
    blst_p1_to_affine(&p1_affine, &BLS12_381_G1);
    blst_p2_to_affine(&p2_affine, &BLS12_381_G2);
    blst_miller_loop(test_input, &p2_affine, &p1_affine);
    
    // Print input/output of each inversion part for this one end-to-end input
    print_inversion_parts_trace(test_input);

    printf("Test 1b: fp6_to_fp2 polynomial formulas (all 8 outputs; focus 0 and 1)\n");
    vec384fp6 t0_fp6_1b;
    inverse_fp12_to_fp6_wrapper_export(t0_fp6_1b, test_input);
    if (!test_fp6_to_fp2_polynomial_formulas(t0_fp6_1b)) {
        printf("FAILED: Test 1b - formula vs implementation mismatch\n");
        return false;
    }
    printf("PASSED: Test 1b (all 8 formulas match inverse_fp6_to_fp2 output)\n\n");

    printf("Test 1: Using miller loop result as input\n");
    if (!test_inversion_parts_chain(test_input)) {
        printf("FAILED: Test 1\n");
        return false;
    }
    printf("PASSED: Test 1\n\n");
    
    // Test 2: Use the result from test 1 as input (chain test)
    printf("Test 2: Using previous result as input (chain test)\n");
    vec384fp12 test_input2;
    blst_fp12_inverse(test_input2, test_input);
    if (!test_inversion_parts_chain(test_input2)) {
        printf("FAILED: Test 2\n");
        return false;
    }
    printf("PASSED: Test 2\n\n");

    printf("Test 2b: fp6_to_fp2 formulas on second fp6 (methodology check)\n");
    vec384fp6 t0_fp6_2;
    inverse_fp12_to_fp6_wrapper_export(t0_fp6_2, test_input2);
    if (!test_fp6_to_fp2_polynomial_formulas(t0_fp6_2)) {
        printf("FAILED: Test 2b\n");
        return false;
    }
    printf("PASSED: Test 2b\n\n");
    
    // Test 3: Verify that a * inverse(a) = 1
    printf("Test 3: Verifying a * inverse(a) = 1\n");
    vec384fp12 inverse_result, product;
    
    // Compute inverse using parts
    vec384fp6 t0_fp6;
    inverse_fp12_to_fp6_wrapper_export(t0_fp6, test_input);
    vec384fp2 ret_fp2, c0, c1, c2;
    inverse_fp6_to_fp2_wrapper_export(ret_fp2, c0, c1, c2, t0_fp6);
    vec384 norm_fp;
    inverse_fp2_to_fp_wrapper_export(norm_fp, ret_fp2);
    vec384 r_fp;
    blst_fp_inverse(r_fp, norm_fp);
    vec384fp2 t1_fp2;
    inverse_fp_to_fp2_wrapper_export(t1_fp2, r_fp, ret_fp2);
    vec384fp6 t1_fp6;
    inverse_fp2_to_fp6_wrapper_export(t1_fp6, t1_fp2, c0, c1, c2);
    inverse_fp6_to_fp12_wrapper_export(inverse_result, t1_fp6, test_input);
    
    // Multiply a * inverse(a) using parts result
    blst_fp12_mul(product, test_input, inverse_result);
    
    // Check if product is 1
    if (!blst_fp12_is_one(product)) {
        printf("FAILED: Test 3 - a * inverse(a) is not 1\n");
        const vec384fp12 *one = blst_fp12_one();
        print_fp12_element("expected (one)", *one);
        print_fp12_element("got (product)", product);
        return false;
    }
    printf("PASSED: Test 3\n\n");
    
    // Test 4: Verify our result matches blst's result
    printf("Test 4: Verifying parts result matches blst_fp12_inverse\n");
    vec384fp12 inverse_blst;
    blst_fp12_inverse(inverse_blst, test_input);
    if (!fp12_equal(inverse_result, inverse_blst)) {
        printf("FAILED: Test 4 - inverse doesn't match blst\n");
        print_fp12_element("parts result", inverse_result);
        print_fp12_element("blst result", inverse_blst);
        return false;
    }
    printf("PASSED: Test 4\n\n");
    
    printf("=== ALL TESTS PASSED ===\n");
    return true;
}

int main() {
    if (test_inversion_parts_multiple()) {
        return 0;
    } else {
        return 1;
    }
}
