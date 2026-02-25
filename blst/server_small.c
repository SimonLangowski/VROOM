/*
 * Single compilation unit for src_small - includes all source files
 * This matches how the original build.sh compiles everything together
 * 
 * Order matters: consts.c must come first to define constants used by others
 */

#ifdef __cplusplus
extern "C" {
#endif

#include "consts.c"
#include "vect.c"
#include "recip.c"
#include "fp12_tower.c"
#include "e1.c"
#include "e2.c"
#include "bulk_addition.c"
#include "multi_scalar.c"
#include "pairing.c"

/* Wrappers for fp12 functions to make them accessible from C++ */
void mul_fp12_export(vec384fp12 ret, const vec384fp12 a, const vec384fp12 b) {
    mul_fp12(ret, a, b);
}

void sqr_fp12_export(vec384fp12 ret, const vec384fp12 a) {
    sqr_fp12(ret, a);
}

void cyclotomic_sqr_fp12_export(vec384fp12 ret, const vec384fp12 a) {
    cyclotomic_sqr_fp12(ret, a);
}

void conjugate_fp12_export(vec384fp12 a) {
    conjugate_fp12(a);
}

void frobenius_map_fp12_export(vec384fp12 ret, const vec384fp12 a, size_t n) {
    frobenius_map_fp12(ret, a, n);
}

void mul_by_xy00z0_fp12_export(vec384fp12 ret, const vec384fp12 a, const vec384fp6 xy00z0) {
    mul_by_xy00z0_fp12(ret, a, xy00z0);
}

void mul_fp6_export(vec384fp6 ret, const vec384fp6 a, const vec384fp6 b) {
    mul_fp6(ret, a, b);
}

void sqr_fp4_export(vec384fp4 ret, const vec384x a0, const vec384x a1) {
    sqr_fp4(ret, a0, a1);
}

void sqr_fp6_export(vec384fp6 ret, const vec384fp6 a) {
    sqr_fp6(ret, a);
}

void mul_by_u_plus_1_fp2_export(vec384x ret, const vec384x a) {
    mul_by_u_plus_1_fp2(ret, a);
}

void neg_fp6_export(vec384fp6 ret, const vec384fp6 a) {
    neg_fp6(ret, a);
}

void inverse_fp6_export(vec384fp6 ret, const vec384fp6 a) {
    inverse_fp6(ret, a);
}/* Wrapper for static inline from_fp to make it accessible from C++ */
void blst_from_fp(vec384 ret, const vec384 a) {
    from_mont_384(ret, a, BLS12_381_P, BLS12_381_p0);
}

/* Public wrappers for div_by_zz and div_by_z */
void blst_div_by_zz(limb_t val[]) {
    div_by_zz(val);
}

void blst_div_by_z(limb_t val[]) {
    div_by_z(val);
}

/* Public wrapper for POINTonE1_mult_glv */
void blst_p1_mult_glv(POINTonE1 *out, const POINTonE1 *in, const pow256 SK) {
    POINTonE1_mult_glv(out, in, SK);
}

/* Public wrappers for precompute and gather_booth */
void blst_p1_precompute_w5(POINTonE1 table[16], const POINTonE1 *point) {
    POINTonE1_precompute_w5(table, point);
}

bool_t blst_p1_gather_booth_w5(POINTonE1 *p, const POINTonE1 table[16], limb_t booth_idx) {
    return POINTonE1_gather_booth_w5(p, table, booth_idx);
}

#ifdef __cplusplus
}
#endif
