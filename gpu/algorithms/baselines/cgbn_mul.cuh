#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
// #include "CGBN/include/cgbn/cgbn.h"
#include "../gpu/algorithms/baselines/CGBN/include/cgbn/cgbn.h" // hacky for now

// Compute how many iterations this block should do
__device__ __forceinline__ int block_offsets(int &batches, int total_batches) {
    // The number of multiprocessors might not divide evenly
    int remaining_batches = total_batches - batches * gridDim.x;
    // Each smaller index does batches
    int offset = blockIdx.x * batches;
    if (blockIdx.x < remaining_batches) {
        // I need to do an extra
        batches += 1;
        // Adjust for each smaller index doing batches + 1
        offset += blockIdx.x;
    } else {
        // only the first remaining_batches do an extra
        offset += remaining_batches;
    }
    return offset;
}

// copied from powm_odd.cu
class params {
    public:
        // parameters used by the CGBN context
        static const uint32_t TPB=0;                     // get TPB from blockDim.x  
        static const uint32_t MAX_ROTATION=4;            // good default value
        static const uint32_t SHM_LIMIT=0;               // no shared mem available
        static const bool     CONSTANT_TIME=false;       // constant time implementations aren't available yet

        static const uint32_t TPI=PYREPLACE_TPI;                   // threads per instance
        static const uint32_t BITS=PYREPLACE_MODBITS;                 // instance size
};

typedef cgbn_context_t<params::TPI, params>         context_t;
typedef cgbn_env_t<context_t, params::BITS> env_t;
typedef typename env_t::cgbn_t      bn_t;

// the actual kernel
__device__ __forceinline__ void kernel_mul( env_t bn_env,  cgbn_mem_t<params::BITS>* a_ptr,  cgbn_mem_t<params::BITS>* b_ptr, cgbn_mem_t<params::BITS>* c_ptr,  bn_t &modulus,  int instance) {
    
    uint32_t np0 = PYREPLACE_MONT_FACTOR;                // for fairness, assume this is precomputed
    
    bn_t a, b, c;
    cgbn_load(bn_env, a, a_ptr + instance);  // load my instance's a value
    cgbn_load(bn_env, b, b_ptr + instance);  // load my instance's b value

    // transform one number to montgomery space, get factor, and check correctness so it can be fixed in python
    // uint32_t np0 = cgbn_bn2mont(bn_env, a_mont, a, modulus);
    // uint32_t np0 = cgbn_bn2mont(bn_env, b_mont, b, modulus);
    // printf("np0 %u %u\n",np0, np1);
    // printf("%d %d a %u b %u m %u\n", threadIdx.x, blockIdx.x, a._limbs[0], b._limbs[0], modulus._limbs[0]);
    // cgbn_mont_mul(bn_env, c_mont, a_mont, b_mont, modulus, np0);                       // c=a*b
    // cgbn_mont2bn(bn_env, c, c_mont, modulus, np0);

    int pow = PYREPLACE_POW;
    for (int i = 0; i < pow; i++) {
        cgbn_mont_mul(bn_env, c, a, b, modulus, np0);
        a = c;
    }
    cgbn_store(bn_env, c_ptr + instance, c);   // store product
}

template <uint32_t limbs>
__device__ void mult_cgbn( cgbn_mem_t<params::BITS>* a_ptr,  cgbn_mem_t<params::BITS>* b_ptr,  cgbn_mem_t<params::BITS>* modulus_ptr, cgbn_mem_t<params::BITS>* c_ptr, int total_chunks, int chunks_per_block) {
    // verify alignment
    // static_assert(sizeof(bn_t) == 4 * limbs);
    
    int my_iterations = chunks_per_block;
    int my_offset = block_offsets(my_iterations, total_chunks);
    int block_per = blockDim.x / params::TPI;
    int block_offset = block_per * my_offset;

    context_t      bn_context = context_t();              // ruct a context, ignoring errors
    env_t          bn_env(bn_context.env<env_t>());                     // ruct an environment
    bn_t modulus;
    cgbn_load(bn_env, modulus, modulus_ptr);
    for (int i = 0; i < my_iterations; i++) {
        int instance = block_offset + ((threadIdx.x / params::TPI) * my_iterations) + i; // doesn't handle fractional iterations
        // printf("%d %d %d %d\n", blockIdx.x, threadIdx.x, instance, total_chunks * chunks_per_block);
        kernel_mul(bn_env, a_ptr, b_ptr, c_ptr, modulus, instance);
    }
}

extern "C"
__launch_bounds__(PYREPLACE_BLOCKSIZE, PYREPLACE_BLOCKSPER) __global__
void MULT_CGBN( cgbn_mem_t<params::BITS>* a,  cgbn_mem_t<params::BITS>* b,  cgbn_mem_t<params::BITS>* modulus_ptr, cgbn_mem_t<params::BITS>* c, int total_chunks, int chunks_per_block) {
    mult_cgbn<PYREPLACE_LIMBS>(a, b, modulus_ptr, c, total_chunks, chunks_per_block);
}