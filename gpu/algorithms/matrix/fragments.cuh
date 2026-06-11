#pragma once
#include <stdint.h>

// TODO: Change to M, K as inputs?
// K_LIMBS = (K + 3) / 4;
// M_LIMBS = (M + 7) / 8;

__device__ __forceinline__ constexpr int ceil_div(const int a, const int b) {
    return (a + b - 1) / b;
}

__device__ __forceinline__ constexpr int padded(const int a, const int b) {
    return b * ceil_div(a, b);
}

// cpp does not allow zero length array, but we would like logical zero length for consistency with variations
// The logical length (used for loop iteration) will be zero so any allocation will be unused, and may be discarded by the compiler
__device__ __forceinline__ constexpr int nonzero(int n) {
    return (n == 0) ? 1 : n;
}

template <class Fragment, class Fragment2, int dst_slice_offset, int src_slice_offset, int slice_length>
__device__ __forceinline__ void copy_slice(Fragment &dst, const Fragment2 &src) {
    static_assert(dst.ROW_LIMBS == src.ROW_LIMBS);
    static_assert(src.COL_LIMBS >= src_slice_offset + slice_length);
    static_assert(dst.COL_LIMBS >= dst_slice_offset + slice_length);
    #pragma unroll
    for (int i = 0; i < src.ROW_LIMBS; i++) {
        #pragma unroll
        for (int j = 0; j < slice_length; j++) {
            dst.set(i, j + dst_slice_offset, src.get(i, j + src_slice_offset));
        }
    }
}

template <class Fragment, class Fragment2>
__device__ __forceinline__ void copy(Fragment &dst, const Fragment2 &src) {
    static_assert(src.COL_LIMBS == dst.COL_LIMBS);
    copy_slice<Fragment, Fragment2, 0, 0, Fragment2::COL_LIMBS>(dst, src);
}

// A class holding a fragment for a M by K matrix in registers
template <int M, int K, int RL, int CL, int STRIDE=32>
class RegisterFragmentInternal {
    public:
        static constexpr int ROWS = M;
        static constexpr int COLS = K;
        static constexpr int ROW_LIMBS = RL;
        static constexpr int COL_LIMBS = CL;
    private:
        uint32_t r[nonzero(ROW_LIMBS)][nonzero(COL_LIMBS)];
        __device__ __forceinline__ int get_index(int row, int col) const {
            return (((row * COL_LIMBS) + col) * STRIDE) + (threadIdx.x % STRIDE);
        }
    public:
        __device__ __forceinline__ RegisterFragmentInternal() {};
        // for abstractness; this constructor is ignored for register fragments
        __device__ __forceinline__ RegisterFragmentInternal(const int offset) {};

        __device__ RegisterFragmentInternal(const uint32_t* pre_ptr, const int offset) {
            #pragma unroll
            for (int i = 0; i < ROW_LIMBS; i++) {
                for (int j = 0; j < COL_LIMBS; j++) {
                    set(i, j, pre_ptr[offset + get_index(i, j)]);
                }
            }
        }
        
        __device__ __forceinline__ uint32_t get(int row, int col) const {
            return r[row][col];
        }
        __device__ __forceinline__ void set(int row, int col, uint32_t val) {
            r[row][col] = val;
        }

        template <class Other>
        __device__ __forceinline__ void copy_from(Other &src) {
            copy(*this, src);
        }

        template <int dst_offset, int src_offset, int length, class Other>
        __device__ __forceinline__ void copy_from_slice(Other &src) {
            static_assert(dst_offset % 4 == 0);
            static_assert(src_offset % 4 == 0);
            copy_slice<RegisterFragment, Other, dst_offset / 4, src_offset / 4, ceil_div(length, 4)>(*this, src);
        }
};

template<int M, int K, int STRIDE=32>
using RegisterFragment = RegisterFragmentInternal<M, K, ceil_div(M, 8), ceil_div(K, 4), STRIDE>;

// A class holding a fragment for a M by K matrix in shared memory 
template<int M, int K, int STRIDE=32>
class SharedFragment {
    public:
        static constexpr int ROWS = M;
        static constexpr int COLS = K;
        static constexpr int ROW_LIMBS = ceil_div(M, 8);
        static constexpr int COL_LIMBS = ceil_div(K, 4);
    private:
        const int sharedOffset;
        __device__ __forceinline__ int get_index(int row, int col) const {
            return sharedOffset + (((row * COL_LIMBS) + col) * STRIDE) + (threadIdx.x % STRIDE);
        }
    public:
        __device__ __forceinline__ SharedFragment(const int offset) : sharedOffset(offset) {}
                
        __device__ SharedFragment(const uint32_t* pre_ptr, const int offset) : sharedOffset(offset) {
            #pragma unroll
            for (int i = 0; i < ROW_LIMBS; i++) {
                for (int j = 0; j < COL_LIMBS; j++) {
                    set(i, j, pre_ptr[get_index(i, j)]); // Careful: get index adds offset so need absolute index
                }
            }
        }

        __device__ __forceinline__ uint32_t get(int row, int col) const {
            extern __shared__ uint32_t shared[];
            return shared[get_index(row, col)];
        }
        __device__ __forceinline__ void set(int row, int col, uint32_t val) {
            extern __shared__ uint32_t shared[];
            shared[get_index(row, col)] = val;
        }

        template <class Other>
        __device__ __forceinline__ void copy_from(Other &src) {
            copy(*this, src);
        }
};


// Do I still want rows and cols to be available?  Do I want a inheritance?
// I want to have logical rows < physical rows for the purpose of small batching

// using RegisterFragmentMatrix = RegisterFragmentInternal<ceil_div(M, 8), ceil_div(K, 4), STRIDE>;

// template<int M, int K, int STRIDE=32>
// class RegisterFragmentMatrix : public RegisterFragment<ceil_div(M, 8), ceil_div(K, 4), STRIDE> {
//     public:
//         constexpr int ROWS = M;
//         constexpr int COLS = N;
// };