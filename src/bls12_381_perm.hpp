#pragma once

#include "bls12_381_perm_params.hpp"

namespace bls12_381_r1_perm {

/** BLS12-381 r1: 8×8 cyclic perm, skip k correction (matches Python skip_correction=True). */
template <const rns_perm::FixedPermTables<BLS12_381_PERM_DIM> &fixed_tables>
struct FixedPermWideReducer {
    static constexpr bool SKIP_K_CORRECTION = true;

    template <int batch>
    inline __attribute__((always_inline)) void rns_reduce_raw_batch(
        std::array<AVXVector<BLS12_381_PADDED_LIMBS>, batch> *residues1,
        std::array<AVXVector<BLS12_381_PADDED_LIMBS>, batch> &out_hi,
        std::array<AVXVector<BLS12_381_PADDED_LIMBS>, batch> &out_lo) const {
        static_assert(BLS12_381_PADDED_LIMBS == 8, "BLS12-381 perm requires 8-wide AVX storage");
        rns_perm::rns_reduce_perm_batch_wide<BLS12_381_PERM_DIM, batch, false, true, true>(
            fixed_tables, *residues1, out_hi, out_lo);
    }
};

inline constexpr FixedPermWideReducer<bls12_381_r1_fixed_perm_tables> r1_fixed_perm_wide_reducer{};

} // namespace bls12_381_r1_perm
