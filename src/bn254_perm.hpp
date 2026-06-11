#pragma once

#include "bn254_perm_params.hpp"

namespace bn254_r1_perm {

/** IntRNS4 perm wide: AVXVector<8> in/out, advance_perm_mt_c matmul, optional k correction. */
template <const rns_perm::FixedPermTables<BN254_PERM_DIM> &fixed_tables, bool ApplyCorrection>
struct FixedPermWideReducer {
    static constexpr bool SKIP_K_CORRECTION = !ApplyCorrection;

    template <int batch>
    inline __attribute__((always_inline)) void rns_reduce_raw_batch(
        std::array<AVXVector<BN254_PADDED_LIMBS>, batch> *residues1,
        std::array<AVXVector<BN254_PADDED_LIMBS>, batch> &out_hi,
        std::array<AVXVector<BN254_PADDED_LIMBS>, batch> &out_lo) const {
        static_assert(BN254_PADDED_LIMBS == 8, "perm reducer requires 8-wide AVX storage");
        rns_perm::rns_reduce_perm_batch_wide<BN254_PERM_DIM, batch, ApplyCorrection>(
            fixed_tables, *residues1, out_hi, out_lo);
    }
};

using R1FixedPermWideReducer = FixedPermWideReducer<bn254_r1_fixed_perm_tables, false>;
using R2FixedPermWideReducer = FixedPermWideReducer<bn254_r2_fixed_perm_tables, true>;

inline constexpr R1FixedPermWideReducer r1_fixed_perm_wide_reducer{};
inline constexpr R2FixedPermWideReducer r2_fixed_perm_wide_reducer{};

} // namespace bn254_r1_perm
