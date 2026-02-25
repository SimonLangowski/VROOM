
#ifndef __INVERSION_PARTS_HPP__
#define __INVERSION_PARTS_HPP__

extern "C" {
#include "fields.h"
}

extern "C" {
void sqr_fp6_export(vec384fp6 ret, const vec384fp6 a);
// mul_fp6_export is already declared in fields.h
void neg_fp6_export(vec384fp6 ret, const vec384fp6 a);
void mul_by_u_plus_1_fp2_export(vec384x ret, const vec384x a);
void blst_fp_inverse(vec384 out, const vec384 inp);
}

static void inverse_fp12_to_fp6(vec384fp6 ret, const vec384fp12 a) {
    vec384fp6 t0, t1;
    sqr_fp6_export(t0, a[0]);
    sqr_fp6_export(t1, a[1]);
    mul_by_u_plus_1_fp2_export(t1[2], t1[2]);
    sub_fp2(t0[0], t0[0], t1[2]);
    sub_fp2(t0[1], t0[1], t1[0]);
    sub_fp2(t0[2], t0[2], t1[1]);
    vec_copy(ret[0], t0[0], sizeof(t0[0]));
    vec_copy(ret[1], t0[1], sizeof(t0[1]));
    vec_copy(ret[2], t0[2], sizeof(t0[2]));
}

static void inverse_fp6_to_fp12(vec384fp12 ret, const vec384fp6 t1, const vec384fp12 a) {
    mul_fp6_export(ret[0], a[0], t1);
    mul_fp6_export(ret[1], a[1], t1);
    neg_fp6_export(ret[1], ret[1]);
}

static void inverse_fp6_to_fp2(vec384fp2 ret, vec384fp2 c0_out, vec384fp2 c1_out, vec384fp2 c2_out, const vec384fp6 a) {
    vec384x c0, c1, c2, t0, t1;

    /* c0 = a0^2 - (a1*a2)*(u+1) */
    sqr_fp2(c0, a[0]);
    mul_fp2(t0, a[1], a[2]);
    mul_by_u_plus_1_fp2_export(t0, t0);
    sub_fp2(c0, c0, t0);

    /* c1 = a2^2*(u+1) - (a0*a1) */
    sqr_fp2(c1, a[2]);
    mul_by_u_plus_1_fp2_export(c1, c1);
    mul_fp2(t0, a[0], a[1]);
    sub_fp2(c1, c1, t0);

    /* c2 = a1^2 - a0*a2 */
    sqr_fp2(c2, a[1]);
    mul_fp2(t0, a[0], a[2]);
    sub_fp2(c2, c2, t0);

    /* Output c0, c1, c2 */
    vec_copy(c0_out[0], c0[0], sizeof(c0[0]));
    vec_copy(c0_out[1], c0[1], sizeof(c0[1]));
    vec_copy(c1_out[0], c1[0], sizeof(c1[0]));
    vec_copy(c1_out[1], c1[1], sizeof(c1[1]));
    vec_copy(c2_out[0], c2[0], sizeof(c2[0]));
    vec_copy(c2_out[1], c2[1], sizeof(c2[1]));

    /* (a2*c1 + a1*c2)*(u+1) + a0*c0 */
    mul_fp2(t0, c1, a[2]);
    mul_fp2(t1, c2, a[1]);
    add_fp2(t0, t0, t1);
    mul_by_u_plus_1_fp2_export(t0, t0);
    mul_fp2(t1, c0, a[0]);
    add_fp2(t0, t0, t1);
    vec_copy(ret[0], t0[0], sizeof(t0[0]));
    vec_copy(ret[1], t0[1], sizeof(t0[1]));
}

static void inverse_fp2_to_fp6(vec384fp6 ret, const vec384fp2 t1, const vec384fp2 c0, const vec384fp2 c1, const vec384fp2 c2) {
    mul_fp2(ret[0], c0, t1);
    mul_fp2(ret[1], c1, t1);
    mul_fp2(ret[2], c2, t1);
}

// inverse_fp2_to_fp: takes fp2 (2 elements) -> outputs norm (1 element) = 1 output
// Computes a^2 + b^2 (the norm, without the reciprocal since that's not quadratic)
static void inverse_fp2_to_fp(vec384 norm_out, const vec384fp2 a) {
    vec384 t0, t1;
    sqr_fp(t0, a[0]);
    sqr_fp(t1, a[1]);
    add_fp(norm_out, t0, t1);
}

// inverse_fp_to_fp2: takes fp r (1 element) and fp2 (a, b) (2 elements) = 3 inputs -> outputs fp2 (2 elements)
// Computes a*r and -b*r (multiplication parts, assuming r is already computed elsewhere)
static void inverse_fp_to_fp2(vec384fp2 ret, const vec384 r, const vec384fp2 a) {
    mul_fp(ret[0], a[0], r);
    mul_fp(ret[1], a[1], r);
    neg_fp(ret[1], ret[1]);
}

// Non-static wrapper functions for use in empirical_basis
void inverse_fp12_to_fp6_wrapper_export(vec384fp6 ret, const vec384fp12 a) {
    inverse_fp12_to_fp6(ret, a);
}

void inverse_fp6_to_fp12_wrapper_export(vec384fp12 ret, const vec384fp6 t1, const vec384fp12 a) {
    inverse_fp6_to_fp12(ret, t1, a);
}

void inverse_fp6_to_fp2_wrapper_export(vec384fp2 ret, vec384fp2 c0_out, vec384fp2 c1_out, vec384fp2 c2_out, const vec384fp6 a) {
    inverse_fp6_to_fp2(ret, c0_out, c1_out, c2_out, a);
}

/* Part 1: a (fp6) -> c0, c1, c2 (fp2 each). Quadratic in a. */
void inverse_fp6_to_fp2_c0_c1_c2_wrapper_export(vec384fp2 c0_out, vec384fp2 c1_out, vec384fp2 c2_out, const vec384fp6 a) {
    vec384x c0, c1, c2, t0;
    sqr_fp2(c0, a[0]);
    mul_fp2(t0, a[1], a[2]);
    mul_by_u_plus_1_fp2_export(t0, t0);
    sub_fp2(c0, c0, t0);
    sqr_fp2(c1, a[2]);
    mul_by_u_plus_1_fp2_export(c1, c1);
    mul_fp2(t0, a[0], a[1]);
    sub_fp2(c1, c1, t0);
    sqr_fp2(c2, a[1]);
    mul_fp2(t0, a[0], a[2]);
    sub_fp2(c2, c2, t0);
    vec_copy(c0_out[0], c0[0], sizeof(c0[0]));
    vec_copy(c0_out[1], c0[1], sizeof(c0[1]));
    vec_copy(c1_out[0], c1[0], sizeof(c1[0]));
    vec_copy(c1_out[1], c1[1], sizeof(c1[1]));
    vec_copy(c2_out[0], c2[0], sizeof(c2[0]));
    vec_copy(c2_out[1], c2[1], sizeof(c2[1]));
}

/* Part 2: c0, c1, c2, a (fp2 each + fp6) -> ret (fp2). Bilinear in (c,a), so quadratic for empirical basis. */
void inverse_fp6_to_fp2_ret_wrapper_export(vec384fp2 ret, const vec384fp2 c0, const vec384fp2 c1, const vec384fp2 c2, const vec384fp6 a) {
    vec384x t0, t1;
    mul_fp2(t0, c1, a[2]);
    mul_fp2(t1, c2, a[1]);
    add_fp2(t0, t0, t1);
    mul_by_u_plus_1_fp2_export(t0, t0);
    mul_fp2(t1, c0, a[0]);
    add_fp2(t0, t0, t1);
    vec_copy(ret[0], t0[0], sizeof(t0[0]));
    vec_copy(ret[1], t0[1], sizeof(t0[1]));
}

void inverse_fp2_to_fp6_wrapper_export(vec384fp6 ret, const vec384fp2 t1, const vec384fp2 c0, const vec384fp2 c1, const vec384fp2 c2) {
    inverse_fp2_to_fp6(ret, t1, c0, c1, c2);
}

void inverse_fp2_to_fp_wrapper_export(vec384 norm_out, const vec384fp2 a) {
    inverse_fp2_to_fp(norm_out, a);
}

void inverse_fp_to_fp2_wrapper_export(vec384fp2 ret, const vec384 r, const vec384fp2 a) {
    inverse_fp_to_fp2(ret, r, a);
}

#endif /* __INVERSION_PARTS_HPP__ */