#pragma once

#include "bounded_ring.hpp"

// BN254 G1 BoundedRing parameters (IntRNS4, 6×46-bit residues, b=3).
// Same layout and bounds as secp256k1; only curve constant b differs.

namespace bn254_params {

static const char *BN254_P =
    "21888242871839275222246405745257275088696311157297823662689037894645226208583";

constexpr int BN254_BITS = 256;
constexpr int BN254_TRUE_LIMBS = 6;
constexpr int BN254_PADDED_LIMBS = 8;
constexpr int BN254_LIMBS = BN254_TRUE_LIMBS;
constexpr int BN254_ELEMENT_BITS = 46;
constexpr int MIN_RNS = -1386;
constexpr int MAX_RNS = 190;
constexpr int MAX_ADD = 6305;
constexpr int LOG_MULTIPLES = 13;
constexpr int BN254_B = 3;
constexpr int BN254_3B = 3 * BN254_B;
constexpr int BN254_6B = 6 * BN254_B;
constexpr int RNS_HALF_MAX_NEG = 18;

constexpr int BASE_REDUNDANT_MULTIPLE = 4;
constexpr int BASE_REDUNDANT_MULTIPLE_NO_K = BN254_TRUE_LIMBS * 4;

using Bn254Ring = BoundedRing<
    BN254_BITS, BN254_TRUE_LIMBS, BN254_ELEMENT_BITS,
    MIN_RNS, MAX_RNS, LOG_MULTIPLES, BN254_3B, MAX_ADD, BN254_6B,
    BASE_REDUNDANT_MULTIPLE>;

using Bn254RingFixedPerm = BoundedRing<
    BN254_BITS, BN254_TRUE_LIMBS, BN254_ELEMENT_BITS,
    MIN_RNS, MAX_RNS, LOG_MULTIPLES, BN254_3B, MAX_ADD, BN254_6B,
    BASE_REDUNDANT_MULTIPLE_NO_K, BN254_PADDED_LIMBS, ChangeBaseVariation::FixedPerm>;

using Bn254RingProduction = Bn254RingFixedPerm;

} // namespace bn254_params
