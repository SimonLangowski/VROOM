#pragma once
#include <stdint.h>
#include "../matrix/fragments.cuh"
#include "../matrix/parts.cuh"
#include "../matrix/epilogue.cuh"
#include "../baselines/montgomery.cuh"

#define ALL_SHARED 0

constexpr __device__ __host__ __forceinline__ int log2(int x) {
    int result = 0;
    while (x >>= 1) {
        result++;
    }
    return result;
}

template<int modbits>
__device__ __forceinline__ uint64_t full_reduce_start(u96 &acc, uint32_t modulus) {
    constexpr int shift = (2*modbits - 32);
    constexpr uint32_t mask = ((uint64_t)(1) << shift) - 1;
    uint32_t upper = acc.hi >> shift;
    uint32_t rem = acc.hi & mask;

    // compiler can flip sign in madd so no-op?
    uint32_t t = -modulus; // 2^32 - (2^32 - t) = t
   
    // if 32 montgomery need upper conditional correction (however shift and mask will be free)
    if (modbits == 32) {
        // 33 bit number doesn't fit
        uint64_t tmp = (uint64_t)(rem) + (upper * t);
        tmp -= (tmp >= modulus) * modulus;
        rem = tmp;
        // uint32_t tmp = upper * modulus;
        // rem = rem - tmp + ((tmp > rem) * modulus); // correct underflow
    } else {
        // 28 bit rem, 29 bit sum, 61 bits in acceptable redundance range
        rem = rem + upper * t; // e.g 7 bits is 128 accumulations (e.g 4096 bit modulus) allows 21 bits for t. which is a million possible t
        // rem = rem - upper * modulus;
    }
    return ((uint64_t)(rem) << 32) + acc.lo;
}

// Stride columns as even, odd to avoid bank conflicts?, rather than first half/second half of column. 
// On the other hand, consecutive values allows larger shared memory read blocks
template <int MPT, int TPI, int BPT, int BROADCASTS, class AType>
__device__ __forceinline__ uint32_t broadcast_get(const AType &residues, int j) {
    // todo tpi == 2 shfl xor sync when > 1/2 and prearrange matrix
    // debug_print1<0>("tpi", TPI);
    // debug_print1<0>("value", j / TPI);
    // debug_print1<0>("thread", j % TPI);
    // debug_print1("BPT", BPT);
    constexpr int res_read_per_thread = BROADCASTS;
    int group = (threadIdx.x & (BPT - 1));
    if ((TPI <= 32) && (!ALL_SHARED)) {
        static_assert(AType::COL_LIMBS == 1);
        int idx = res_read_per_thread * group + j;
        if (MPT > 1) {
            idx = idx / TPI; // untested
        }
        idx *= BPT;
        debug_print1<0>("j", idx);
        debug_print1<1>("j", idx);
        return __shfl_sync(0xffffffff, residues.get(0, 0), idx, TPI);
    } else {
        extern __shared__ uint32_t shared[];
        // residues for each subset, add 1 to avoid conflicts
        // = BROADCASTS?
        int bank_offset = (res_read_per_thread + 4) * group;
        // int bank_offset = (res_read_per_thread) * group;
        int idx = bank_offset + j;
        debug_print1<0>("j", idx);
        debug_print1<1>("j", idx);
        return shared[idx];
    }
}

template <int MPT, int TPI, int BPT, int BROADCASTS, class AType>
__device__ __forceinline__ uint32_t broadcast_set(const AType &residues) {
    extern __shared__ uint32_t shared[];
    if ((TPI > 32) || (ALL_SHARED)) {
        static_assert(AType::COL_LIMBS == 1);
        // residues for each subset, add 1 to avoid conflicts
        // = BROADCASTS?
        constexpr int res_read_per_thread = BROADCASTS;
        if ((threadIdx.x % BPT == 0) && (threadIdx.x >= 2*BPT)) {
            int j_val = threadIdx.x / BPT;
            int group = (j_val - 2) / res_read_per_thread;
            // also add bank conflict destriding
            // do complicated stuff here so reading is sequential
            int idx = j_val + 4*group;
            // int idx = j_val;
            shared[idx] = residues.get(0, 0);
            debug_print1<4>("jwrite4", idx);
            debug_print1<6>("jwrite6", idx);
            debug_print1<8>("jwrite8", idx);
        }
        // shared[threadIdx.x] = residues.get(0, 0);
        __syncthreads();
    }
}

// Should I also make a u64 version?
template<int TPI, int iterations=0>
__device__ __forceinline__ u96 reduce_add(u96 v) {
    #pragma unroll
    for (int iter = 0; iter < iterations; iter++) {
        debug_print1<0>("vlo0", v.lo);
        debug_print1<0>("vhi0", v.hi);
        debug_print1<1>("vlo1", v.lo);
        debug_print1<1>("vhi1", v.hi);
        int shift = 1 << iter;
        uint32_t lo2 = __shfl_xor_sync(0xffffffff, v.lo, shift, TPI);
        uint64_t hi2 = __shfl_xor_sync(0xffffffff, v.hi, shift, TPI);
        uint64_t s1 = (uint64_t)(v.lo) + lo2;
        v.hi += hi2 + (s1 >> 32);
        v.lo = s1;
    }
    return v;
}

template<int T_IN, int T_OUT, int MPT_IN, int MPT_OUT, int TPI, int BPT, class ModType>
class RNSVecReduce {
    // Each thread has MPT (moduli per thread) columns of the B matrix
    // Each column has M1 rows (input # of moduli, but could be  halved to allow addition later)
    static constexpr int BROADCASTS = ceil_div(T_IN, BPT);
    using AType = RegisterFragmentInternal<1, T_IN, 1, MPT_IN, TPI>;
    using BType = RegisterFragmentInternal<BROADCASTS, T_OUT+2, BROADCASTS, MPT_OUT, TPI>;
    using CType = RegisterFragmentInternal<1, T_OUT, 1, MPT_OUT, TPI>;
    const BType precomputed;
    using EpilogueType = RegisterFragmentInternal<1,T_OUT,1,MPT_OUT, TPI>;
    const EpilogueType epilogue;

    public:
    __device__ RNSVecReduce(const uint32_t* pre_ptr) : precomputed(pre_ptr, padded(T_OUT, TPI)), epilogue(pre_ptr, 0) {
    }

    template<class EE, class ModReduction>
    __device__ __forceinline__ void reduce(const AType (&residues), CType (&out_residues), const EE &ee, const ModReduction &mod) const {
        u96 acc[MPT_OUT];

        broadcast_set<MPT_IN, TPI, BPT, BROADCASTS>(residues);
       
        #pragma unroll
        for (int j = 0; j < MPT_OUT; j++) {
            uint64_t w = ee.extra(0, 0, j);
            if (threadIdx.x % BPT != 0) {
                w = 0; // could be a 32 bit select before mult
            }
            acc[j].lo = w;
            acc[j].hi = w >> 32;
        }
    
        // every index required for column dot
        #pragma unroll
        for (int k = 0; k < BROADCASTS; k++) {
            // load a element
            uint32_t a = broadcast_get<MPT_IN, TPI, BPT, BROADCASTS>(residues, k + 2); // first two slots are reserved for k correction
            debug_print1<0>("a", a);
            debug_print1<1>("a", a);
            // for each column I'm responsible for
            #pragma unroll
            for (int j = 0; j < MPT_OUT; j++) {
                // load b element (should be stored in registers)
                uint32_t b_pre = precomputed.get(k, j);
                mul_split_acc(a, b_pre, acc[j]);
            }
        }
        // compute k
        // thread 1: low 32 * 32 bits
        // thread 0: high 32 * 32 bits

        debug_printvec<TPI>("bpre", precomputed);

        #pragma unroll
        for (int j = 0; j < MPT_OUT; j++) {
            acc[j] = reduce_add<TPI, log2(BPT)>(acc[j]);
        }
        
        debug_print1<0>("t0acc0lo", acc[0].lo);
        debug_print1<0>("t0acc0hi", acc[0].hi);
        debug_print1<BPT>("t1acc0lo", acc[0].lo);
        debug_print1<BPT>("t1acc0hi", acc[0].hi);
    
        // We will shift by 64; therefore only upper 64 bits of low part will be needed
        // always broadcast from thread 1, even with more than 32 threads.
        uint64_t low_carry = __shfl_sync(0xffffffff, acc[0].hi, BPT, TPI);
        // debug_print1<0>("t1acc0hi", low_carry);
        low_carry += acc[0].lo;
        debug_print1<0>("lc", low_carry);
        // compute k and finish shift
        uint64_t k = acc[0].hi + (low_carry >> 32);
        // broadcast k to all threads
        // broadcast_get(0);
        if ((TPI <= 32)  && (!ALL_SHARED)) {
            k = __shfl_sync(0xffffffff, k, 0, 32);
        } else {
            // an alternative is to duplicated the fixed point correction in every warp.
            extern __shared__ uint64_t shared64[]; // Hacky: use both 0 and 1 entry of uint32_t array but okay because both these threads only handle fixed point.
            if (threadIdx.x == 0) {
                shared64[0] = k; // only thread 0 should write
            }
            __syncthreads(); // ensure write complete
            k = shared64[0]; // all threads read value
        }
        debug_print1("k", k);

        uint32_t k_low = k;
        uint32_t k_high = k >> 32;
    
        // do corrections and reductions
        #pragma unroll
        for (int j = 0; j < MPT_OUT; j++) {
            // 64 by 32 bit product
            uint32_t correction_factor = epilogue.get(0, j);
            // debug_print1<45>("corr", correction_factor);
            // debug_print1<45>("kk", k);
            // debug_print1<45>("vv", acc[j].hi);
            mul_split_acc(k_low, correction_factor, acc[j]);
            acc[j].hi += (uint64_t)(k_high) * (uint64_t)(correction_factor);

            // apply modular reduction for 96 bit number
            auto i = mod.load_modulus(j);
            uint64_t partial = full_reduce_start<ModReduction::MODBITS>(acc[j], i.modulus);
            out_residues.set(0, j, mod.mod_reduce(partial, i));
        }
    }
};

template<int M1, int M2, int MPT1, int MPT2, int TPI, int BPT, template <template <int, int, int> class, int, int> class ModReduceType>
class MontRNSReduceVec {
    using M1RedType = ModReduceType<RegisterFragment, M1, TPI>;
    using M2RedType = ModReduceType<RegisterFragment, M2, TPI>;

    using E2Type = EpilogueExtraM<0, M1>;

    using ReductionType = RNSVecReduce<M1, M2, MPT1, MPT2, TPI, BPT, M2RedType>;
    using ExpansionType = RNSVecReduce<M2, M1, MPT2, MPT1, TPI, BPT, M1RedType>;

    const M1RedType elementwise_reduction1;
    const M2RedType elementwise_reduction2;

    const ReductionType reduction;
    const ExpansionType expansion;

    const E2Type ee;

    public:
    using M1Type = RegisterFragmentInternal<1, M1, 1, MPT1, TPI>;
    using M2Type = RegisterFragmentInternal<1, M2, 1, MPT2, TPI>;

    __device__ MontRNSReduceVec(const uint32_t* pre_ptr1, const uint32_t* moduli1, const uint32_t* pre_ptr2, const uint32_t* moduli2) : elementwise_reduction1(moduli1), elementwise_reduction2(moduli2), reduction(pre_ptr1), expansion(pre_ptr2), ee() {}

    __device__ void mulRed(const M1Type &a1, const M2Type &a2, const M1Type &b1, const M2Type &b2, M1Type &c1, M2Type &c2) const {
        M1Type s1;
        using E1Type = EpilogueExtraAcc<M2Type, 0, M2>;
        mod_multiply_rns(a1, b1, s1, elementwise_reduction1);
        E1Type eeacc(a2, b2);
        reduction.reduce(s1, c2, eeacc, elementwise_reduction2);
        expansion.reduce(c2, c1, ee, elementwise_reduction1);
    }
};