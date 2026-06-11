#pragma once
#include "fragments.cuh"

// Elementwise RNS multiplication.
template<class Fragment, class ModReduce>
__device__ __forceinline__ void mod_multiply_rns(const Fragment &a, const Fragment &b, Fragment &out, const ModReduce &mod) {
    debug_printvec<DEBUG_TPI>("a", a);
    debug_printvec<DEBUG_TPI>("b", b);
    #pragma unroll
    for (int j = 0; j < Fragment::COL_LIMBS; j++) {
        auto m = mod.load_modulus(j);
        #pragma unroll
        for (int i = 0; i < Fragment::ROW_LIMBS; i++) {
            uint64_t c = (uint64_t)(a.get(i, j)) * (uint64_t)(b.get(i, j));
            out.set(i, j, mod.mod_reduce(c, m));
        }
    }
    debug_printvec<DEBUG_TPI>("product", out);
}
