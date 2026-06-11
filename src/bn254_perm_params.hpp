#pragma once

#include "rns_reduce_perm.hpp"
#include "../scripts/bn254_intrns4_params.hpp"

// IntRNS4 perm: 6 pre_permute_mt_c rows × 8 AVX lanes; shift_perm derived at compile time.

constexpr int BN254_TRUE_LIMBS = 6;
constexpr int BN254_PADDED_LIMBS = bn254_intrns4::BN254_PADDED_DIM;
constexpr int BN254_PERM_DIM = BN254_TRUE_LIMBS;
static_assert(bn254_intrns4::BN254_RNS_DIM == BN254_PERM_DIM,
    "export BN254_RNS_DIM must match true limb count");
static_assert(BN254_PERM_DIM < BN254_PADDED_LIMBS,
    "bn254 perm rows must fit in AVX width");

constexpr rns_perm::FixedPermTables<BN254_PERM_DIM> bn254_r1_fixed_perm_tables =
    rns_perm::make_fixed_perm_tables<BN254_PERM_DIM>(
        bn254_intrns4::r1_perm_mat_fixed, bn254_intrns4::r1_correction);

constexpr rns_perm::FixedPermTables<BN254_PERM_DIM> bn254_r2_fixed_perm_tables =
    rns_perm::make_fixed_perm_tables<BN254_PERM_DIM>(
        bn254_intrns4::r2_perm_mat_fixed, bn254_intrns4::r2_correction);
