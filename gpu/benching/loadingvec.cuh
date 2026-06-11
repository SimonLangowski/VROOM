#define DEBUG PYREPLACE_DEBUG
#define DEBUG_TPI PYREPLACE_TPI
// Some directory issues with how things are setup
#include "../gpu/benching/include.cuh"

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

__device__ void load_and_run(const uint32_t* a_ptr1, const uint32_t* b_ptr1, const uint32_t* pre_ptr1, const uint32_t* moduli1, const uint32_t* a_ptr2, const uint32_t* b_ptr2, const uint32_t* pre_ptr2, const uint32_t* moduli2, uint32_t* c_ptr2, int total_chunks, int chunks_per_block) {
    using MCLASS = PYREPLACE_MULT_CLASS;
    const auto MClass = MCLASS(pre_ptr1, moduli1, pre_ptr2, moduli2);
    // Number of elements process per thread
    using M1Type = MCLASS::M1Type;
    using M2Type = MCLASS::M2Type;
    const int M = 1;

    const int thread_limbs1 = M1Type::COL_LIMBS * M;
    const int thread_limbs2 = M2Type::COL_LIMBS * M;
    const int block_limbs1 = blockDim.x * thread_limbs1;
    const int block_limbs2 = blockDim.x * thread_limbs2;
    // const int num_mults = PYREPLACE_CHUNKS;
    // const int base_iterations = num_mults / (gridDim.x * M);
    int my_iterations = chunks_per_block;
    int my_offset = block_offsets(my_iterations, total_chunks);
    const int block_offset1 = block_limbs1 * my_offset;
    const int block_offset2 = block_limbs2 * my_offset;
    // printf("%d %d %d %d\n", blockIdx.x, threadIdx.x, my_offset, my_iterations);
    for (int iter = 0; iter < my_iterations; iter++) {
        M1Type a1;
        M1Type b1;
        M2Type a2;
        M2Type b2;
        int z1 = 0;
        int z2 = 0;
        #pragma unroll
        for (int i = 0; i < M1Type::ROW_LIMBS; i++) {
            #pragma unroll
            for (int k = 0; k < M1Type::COL_LIMBS; k++) {
                int idx = block_offset1 + (block_limbs1 * iter) + (z1 * blockDim.x + threadIdx.x);
                a1.set(i, k, a_ptr1[idx]);
                b1.set(i, k, b_ptr1[idx]);
                z1+=1;
            }
            #pragma unroll
            for (int k = 0; k < M2Type::COL_LIMBS; k++) {
                int idx = block_offset2 + (block_limbs2 * iter) + (z2 * blockDim.x + threadIdx.x);
                a2.set(i, k, a_ptr2[idx]);
                b2.set(i, k, b_ptr2[idx]);
                z2+=1;
            }
        }
        M1Type c1;
        M2Type c2;
        PYREPLACE_APP_NAME(a1, a2, b1, b2, c1, c2, MClass, PYREPLACE_POW);
        // debug_printvec<4>("c2", c2);
        int z = 0;
        #pragma unroll
        for (int i = 0; i < M2Type::ROW_LIMBS; i++) {
            #pragma unroll
            for (int k = 0; k < M2Type::COL_LIMBS; k++) {
                int idx = block_offset2 + (block_limbs2 * iter) + (z* blockDim.x + threadIdx.x);
                c_ptr2[idx] = c2.get(i, k);
                z+=1;
            }
        }
    }
}

extern "C"
__launch_bounds__(PYREPLACE_BLOCKSIZE, PYREPLACE_BLOCKSPER) __global__
void MULT_APP(const uint32_t* a_ptr1, const uint32_t* b_ptr1, const uint32_t* pre_ptr1, const uint32_t* moduli1, const uint32_t* a_ptr2, const uint32_t* b_ptr2, const uint32_t* pre_ptr2, const uint32_t* moduli2, uint32_t* c_ptr2, int total_chunks, int chunks_per_block) {
    load_and_run(a_ptr1, b_ptr1, pre_ptr1, moduli1, a_ptr2, b_ptr2, pre_ptr2, moduli2, c_ptr2, total_chunks, chunks_per_block);
}