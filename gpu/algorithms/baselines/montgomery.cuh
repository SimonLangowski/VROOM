#pragma once
#include <stdint.h>
typedef struct {
    uint32_t lo;
    uint64_t hi;
} u96;

__device__ __forceinline__ void mul_split_acc(const uint32_t x, const uint32_t y, u96 &acc) {
    uint64_t product = (uint64_t)(x) * (uint64_t)(y) + acc.lo; // can never overflow
    uint32_t hi = product >> 32;
    uint32_t lo = product;
    acc.lo = lo;
    acc.hi += hi;
}

template <int limbs>
__device__ void mont_mulreduce(const uint32_t (&a)[limbs], const uint32_t (&b)[limbs], uint32_t (&out)[limbs], const uint32_t mont_factor, const uint32_t (&p)[limbs]) {
    uint32_t u[limbs];
    u96 carry = {0, 0};
    #pragma unroll
    for (int k = 0; k < limbs; k++) {
        u96 ut = carry;
        #pragma unroll
        for (int i = 0; i < k; i++) {
            int j = k - i;
            mul_split_acc(a[i], b[j], ut);
            mul_split_acc(u[i], p[j], ut);
        }
        mul_split_acc(a[k], b[0], ut);
        u[k] = (uint32_t)(ut.lo) * mont_factor;
        mul_split_acc(u[k], p[0], ut);
        uint64_t ut_hi = ut.hi;
        carry.lo = (uint32_t)(ut_hi);
        carry.hi = ut_hi >> 32;
    }
    #pragma unroll
    for (int k = limbs; k < 2*limbs; k++) {
        u96 ut = carry;
        #pragma unroll
        for (int i = k - limbs + 1; i < limbs; i++) {
            int j = k - i;
            mul_split_acc(a[i], b[j], ut);
            mul_split_acc(u[i], p[j], ut);
        }
        out[k - limbs] = (uint32_t)(ut.lo);
        uint64_t ut_hi = ut.hi;
        carry.lo = (uint32_t)(ut_hi);
        carry.hi = ut_hi >> 32;
    }
}

template<int LIMBS>
class MontRedundantBaselineReduce {

    const uint32_t mont_factor;
    uint32_t modulus[LIMBS];

    public:
    __device__ MontRedundantBaselineReduce(const uint32_t* pre_ptr, const uint32_t* mod_ptr, const uint32_t* pre_ptr2, const int precision) : mont_factor(*pre_ptr) {
        #pragma unroll
        for (int i = 0; i < LIMBS; i++) {
            modulus[i] = mod_ptr[i];
        }
    }

    template<class Fragment>
    __device__ void mulRed(const Fragment &a, const Fragment &b, Fragment &c) const {
        static_assert(a.COL_LIMBS == LIMBS);
        #pragma unroll
        for (int i = 0; i < a.ROW_LIMBS; i++) {
            uint32_t al[LIMBS];
            uint32_t bl[LIMBS];
            #pragma unroll
            for (int j = 0; j < LIMBS; j++) {
                al[j] = a.get(i, j);
                bl[j] = b.get(i, j);
            }
            uint32_t out[LIMBS];
            mont_mulreduce(al, bl, out, mont_factor, modulus);
            #pragma unroll
            for (int j = 0; j < LIMBS; j++) {
                c.set(i, j, out[j]);
            }
        }
    }
};