#pragma once
#include "./fragments.cuh"
#include "../../debug.cuh"
// Alright, so later can do variation with q correction
// And variation with no correction? <- Let's make it a separate epilogue class
// Do as variation template if add
// make fragment 2d instead of 1d (hacky (q_correct + 1) is the dimension)

template<int offset, int cols>
class EpilogueExtraM {
    public:
        static constexpr int OFFSET_LIMBS = offset / 4;
        static constexpr int COLS = cols;

        __device__ __forceinline__ uint64_t extra(uint64_t w, int i, int j) const {
            return w;
        }
};

template<class Fragment, int offset, int cols>
class EpilogueExtraAcc {
    public:
        static constexpr int OFFSET_LIMBS = offset / 4;
        static constexpr int COLS = cols;

        const Fragment &am2;
        const Fragment &bm2;

        __device__ EpilogueExtraAcc(const Fragment &a, const Fragment &b) : am2(a), bm2(b) {}

        __device__ __forceinline__ uint64_t extra(uint64_t w, int i, int j) const {
            return w + ((uint64_t)am2.get(i, j) * (uint64_t)bm2.get(i, j));
        }
};

// Epilogue for matrix; unmerged method incurred many unnessecary local registers
template<template <int, int, int> class Fragment, uint32_t K_LIMBS_OUT, class ModReduce, int EVAR = 1>
class MCorrection {
    using CDataType = Fragment<EVAR, K_LIMBS_OUT, 4>;
    private:
        const CDataType correctionData;
        const int precision;
    public:
        static constexpr int E_COLS = EVAR;
        static constexpr int COLS = K_LIMBS_OUT;
        
        __device__ MCorrection(const uint32_t* correction_data, const int offset, const int precision): correctionData(correction_data, offset), precision(precision) {        
        }

        template<uint32_t M_LIMBS>
        class CorrectionData {
            public:
            uint32_t k_cors[M_LIMBS];
        };

        template<int p, uint32_t M_LIMBS, class C, class EE>
        __device__ __forceinline__ CorrectionData<M_LIMBS> compute_corrections(uint32_t (&out0)[M_LIMBS], const uint32_t (&out1)[M_LIMBS], int j, C &out_residues, const EE &ee, const ModReduce &mod) const {
            static_assert(COLS == EE::COLS);
            static_assert(out_residues.ROW_LIMBS == M_LIMBS);
            static_assert(out_residues.COL_LIMBS == CDataType::COL_LIMBS);
            // unpack silly compiler template values
            constexpr int is_tmp = p & 1;
            constexpr int is_half = (p >> 1) & 1;
            if (is_tmp) {
                j = 0; // moduli values will not be used; this prevents errors but the compiler should eliminate the statements as we do not use the results.
            }
            constexpr int k_index = is_half ? 1 : 3;
            auto m = mod.load_modulus(j + EE::OFFSET_LIMBS);
            uint32_t cor_value = correctionData.get(0, j);
            CorrectionData<M_LIMBS> c;
            for (int i = 0; i < M_LIMBS; i++) {
                uint32_t cl = out0[i];
                uint32_t ch;
                if (is_half) {
                    ch = __shfl_xor_sync((uint32_t)(-1), cl, 2, 4);
                } else {
                    ch = out1[i];
                }
                uint64_t w = ((uint64_t)(ch) << 16) + cl;
                
                debug_print1<k_index>("k_raw", w);
                uint32_t k = (w >> precision);
                // broadcast the quotient to all threads in this row
                c.k_cors[i] = __shfl_sync((uint32_t)(-1), k, k_index, 4);
                if(!is_tmp) {
                    w += (uint64_t)(c.k_cors[i]) * (uint64_t)(cor_value);
                    w = ee.extra(w, i, j);
                    out_residues.set(i, j, mod.mod_reduce(w, m));
                }
            }
            debug_print1("k", c.k_cors[0]);
            return c;
        }

        template<uint32_t M_LIMBS, class C, class EE>
        __device__ __forceinline__ void apply_corrections(const CorrectionData<M_LIMBS> &c, uint32_t (&out0)[M_LIMBS], const uint32_t (&out1)[M_LIMBS], int j, C &out_residues, const EE &ee, const ModReduce &mod) const {
            static_assert(out_residues.ROW_LIMBS == M_LIMBS);
            static_assert(out_residues.COL_LIMBS == CDataType::COL_LIMBS);
            auto m = mod.load_modulus(j + EE::OFFSET_LIMBS);
            uint32_t cor_value = correctionData.get(0, j);
            // Apply modular reduction
            #pragma unroll
            for (int i = 0; i < M_LIMBS; i++) {
                uint32_t cl = out0[i];
                uint32_t ch = out1[i];
                uint64_t w = ((uint64_t)(ch) << 16) + cl;
                w += (uint64_t)(c.k_cors[i]) * (uint64_t)(cor_value);
                w = ee.extra(w, i, j);
                out_residues.set(i, j, mod.mod_reduce(w, m));
            }
        }
};