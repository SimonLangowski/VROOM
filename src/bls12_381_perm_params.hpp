#pragma once

#include "rns_reduce_perm.hpp"
#include "../scripts/bls12_381_intrns4_params.hpp"

// BLS12-381 IntRNS4: r1 uses 8×8 pre_permute_mt (cyclic advance), no k correction.

constexpr int BLS12_381_TRUE_LIMBS = 8;
constexpr int BLS12_381_PADDED_LIMBS = bls12_381_intrns4::BLS12_381_PADDED_DIM;
constexpr int BLS12_381_PERM_DIM = BLS12_381_TRUE_LIMBS;
static_assert(bls12_381_intrns4::BLS12_381_RNS_DIM == BLS12_381_PERM_DIM,
    "export BLS12_381_RNS_DIM must match true limb count");
static_assert(BLS12_381_PERM_DIM == BLS12_381_PADDED_LIMBS,
    "BLS12-381 uses full 8-wide residues (no mod-1 padding)");

constexpr rns_perm::FixedPermTables<BLS12_381_PERM_DIM> bls12_381_r1_fixed_perm_tables =
    rns_perm::make_fixed_perm_tables<BLS12_381_PERM_DIM>(
        bls12_381_intrns4::r1_perm_mat_fixed, bls12_381_intrns4::r1_correction);
