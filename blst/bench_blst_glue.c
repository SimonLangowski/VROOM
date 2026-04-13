/*
 * Wrappers so bench_blst_complete.cpp can call blst.h Fp/Fp2 entry points while
 * linking libblst_small (field ops are static inline in fields.h there).
 *
 * Types match bindings/blst.h layout; only fields.h is included to avoid
 * limb_t conflicts with vect.h vs blst.h.
 */
#include <stdbool.h>
#include "fields.h"

typedef struct {
    limb_t l[NLIMBS(384)];
} blst_fp;

typedef struct {
    blst_fp fp[2];
} blst_fp2;

void blst_fp_add(blst_fp *ret, const blst_fp *a, const blst_fp *b)
{
    add_fp(ret->l, a->l, b->l);
}

void blst_fp_mul(blst_fp *ret, const blst_fp *a, const blst_fp *b)
{
    mul_fp(ret->l, a->l, b->l);
}

void blst_fp_sqr(blst_fp *ret, const blst_fp *a)
{
    sqr_fp(ret->l, a->l);
}

void blst_fp_cneg(blst_fp *ret, const blst_fp *a, bool flag)
{
    cneg_fp(ret->l, a->l, (bool_t)flag);
}

void blst_fp2_mul(blst_fp2 *ret, const blst_fp2 *a, const blst_fp2 *b)
{
    mul_fp2(*(vec384x *)((void *)ret), *(const vec384x *)((void *)a),
            *(const vec384x *)((void *)b));
}

void blst_fp2_sqr(blst_fp2 *ret, const blst_fp2 *a)
{
    sqr_fp2(*(vec384x *)((void *)ret), *(const vec384x *)((void *)a));
}
