#pragma once

#include "../cpu/vector/vector_impl.hpp"
#include <array>
#include <cstdint>

namespace rns_perm {

constexpr int MULBITS = 52;

/** Pre-built perm rows + correction; advance_perm derived from N (mt_c or cyclic mt). */
template <int N>
struct FixedPermTables {
    static_assert(N > 0 && N <= 8, "RNS fixed perm requires 0 < dim <= 8");
    std::array<std::array<uint64_t, 8>, N> perm_mat{};
    std::array<uint64_t, 8> correction{};
};

template <int N>
constexpr FixedPermTables<N> make_fixed_perm_tables(
    const std::array<std::array<uint64_t, 8>, N> &perm_mat,
    const std::array<uint64_t, 8> &correction) {
    FixedPermTables<N> t{};
    t.perm_mat = perm_mat;
    t.correction = correction;
    return t;
}

/** advance_perm_mt_c on lanes 0..N, identity on padding lanes (matches Python AVX8_WIDTH). */
template <int N>
constexpr std::array<uint64_t, 8> make_shift_perm_indices() {
    std::array<uint64_t, 8> p{0, 1, 2, 3, 4, 5, 6, 7};
    p[0] = 1;
    for (int i = 0; i < N; ++i) {
        p[i + 1] = static_cast<uint64_t>(((i - 1 + N) % N) + 1);
    }
    return p;
}

template <int N>
inline AVXVector<8> shift_perm_vector() {
    return AVXVector<8>(make_shift_perm_indices<N>());
}

/** advance_perm_mt: cyclic shift on lanes 0..N-1, identity on padding (8×8 BLS12-381 r1). */
template <int N>
constexpr std::array<uint64_t, 8> make_shift_perm_indices_cyclic() {
    static_assert(N > 0 && N <= 8, "RNS cyclic perm requires 0 < dim <= 8");
    std::array<uint64_t, 8> p{};
    for (int j = 0; j < N; ++j) {
        p[j] = static_cast<uint64_t>((j - 1 + N) % N);
    }
    for (int j = N; j < 8; ++j) {
        p[j] = static_cast<uint64_t>(j);
    }
    return p;
}

template <int N>
inline AVXVector<8> shift_perm_vector_cyclic() {
    return AVXVector<8>(make_shift_perm_indices_cyclic<N>());
}

/** Lane-0 k: hi[0] + (lo[0] >> mulbits), broadcast to all lanes. */
inline AVXVector<8> extract_k_broadcast_m128_lane0(
    const AVXVector<8> &out_hi,
    const AVXVector<8> &out_lo,
    int mulbits = MULBITS) {
    ShortVector k_raw_lo = ShortVector::from_avx_lo(out_lo);
    ShortVector k_raw_hi = ShortVector::from_avx_lo(out_hi);
    k_raw_lo = k_raw_lo.srli(mulbits);
    const ShortVector k_short = k_raw_hi.add(k_raw_lo);
    return k_short.broadcast_to_avx<8>();
}

/** Full-width madd52 correction (mulhi/mullo accumulate correction * k). */
inline void add_correction_perm(
    AVXVector<8> &out_hi,
    AVXVector<8> &out_lo,
    const AVXVector<8> &correction) {
    const AVXVector<8> k_vec = extract_k_broadcast_m128_lane0(out_hi, out_lo);
    out_hi = out_hi.mulhi(correction, k_vec);
    out_lo = out_lo.mullo(correction, k_vec);
}

/** Cumulative perm matmul: one shift_perm, vector mulhi/mullo, residues updated in place. */
template <int N, int batch, bool ApplyCorrection = true, bool CyclicAdvance = false, bool SkipFirstAdvance = false>
inline __attribute__((always_inline)) void rns_reduce_perm_batch_wide(
    const FixedPermTables<N> &tables,
    std::array<AVXVector<8>, batch> &residues_in,
    std::array<AVXVector<8>, batch> &out_hi,
    std::array<AVXVector<8>, batch> &out_lo) {
    static_assert(N > 0 && N <= 8, "RNS fixed perm requires 0 < dim <= 8");
    const AVXVector<8> shift_perm =
        CyclicAdvance ? shift_perm_vector_cyclic<N>() : shift_perm_vector<N>();
    #pragma clang loop unroll(enable)
    for (int d = 0; d < N; ++d) {
        const AVXVector<8> row(tables.perm_mat[d]);
        if (!SkipFirstAdvance || d > 0) {
            #pragma clang loop unroll(enable)
            for (int j = 0; j < batch; ++j) {
                residues_in[j] = residues_in[j].permute(shift_perm);
            }
        }
        #pragma clang loop unroll(enable)
        for (int j = 0; j < batch; ++j) {
            out_hi[j] = out_hi[j].mulhi(row, residues_in[j]);
            out_lo[j] = out_lo[j].mullo(row, residues_in[j]);
        }
    }
    if constexpr (ApplyCorrection) {
        const AVXVector<8> correction_v(tables.correction);
        #pragma clang loop unroll(enable)
        for (int j = 0; j < batch; ++j) {
            add_correction_perm(out_hi[j], out_lo[j], correction_v);
        }
    }
}

} // namespace rns_perm
