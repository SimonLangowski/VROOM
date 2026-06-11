#pragma once

template <typename DType, typename DType2, typename MType>
__device__ void bench_modmult2(DType &a1, DType2 &a2, const DType &b1, const DType2 &b2, DType &out1, DType2 &out2, MType &m, int num_mults) {
    for (int i = 0; i < num_mults; i++) {
        m.mulRed(a1, a2, b1, b2, out1, out2);
        a1.copy_from(out1);
        a2.copy_from(out2);
    }
}