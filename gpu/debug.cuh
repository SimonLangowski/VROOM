#pragma once

template<int TPI, class Fragment>
__device__ void debug_printvec(const char* name, const Fragment &f) {
  #if DEBUG
    if ((threadIdx.x < TPI) && (blockIdx.x == 0)) {
      for (int j = 0; j < Fragment::ROW_LIMBS; j++) {
      for (int i = 0; i < Fragment::COL_LIMBS; i++) {
        printf("%s %d %d %d %u\n", name, threadIdx.x, j, i, f.get(j, i));
      }
    }
    }
  #endif
}

template<int tid=0>
__device__ void debug_print1(const char* name, uint64_t v) {
  #if DEBUG
  if ((threadIdx.x == tid) && (blockIdx.x == 0)) {
    printf("%s %lu\n", name, v);
  }
  #endif
}

template<int TPI, int len>
__device__ void debug_printarr(const char* name, const uint32_t (&f)[len]) {
  #if DEBUG
    if ((threadIdx.x < TPI) && (blockIdx.x == 0)) {
      for (int i = 0; i < len; i++) {
        printf("%s %d %d: %u\n", name, threadIdx.x, i, f[i]);
      }
    }
  #endif
}