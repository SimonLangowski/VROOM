#pragma once

#include "bounded_ring.hpp"

namespace bls12_381_params {

static const char *BLS12_381_P =
    "4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787";

constexpr int BLS12_381_BITS = 381;
constexpr int BLS12_381_TRUE_LIMBS = 8;
constexpr int BLS12_381_PADDED_LIMBS = 8;
constexpr int BLS12_381_ELEMENT_BITS = 50;
constexpr int MIN_RNS = -1932;
constexpr int MAX_RNS = 2377;
constexpr int MAX_ADD = 800;
constexpr int LOG_MULTIPLES = 14;
constexpr int BLS12_381_B = 4;
constexpr int BLS12_381_3B = 3 * BLS12_381_B;
constexpr int BLS12_381_6B = 6 * BLS12_381_B;
constexpr int RNS_HALF_MAX_NEG = 48;

constexpr int BASE_REDUNDANT_MULTIPLE_MATRIX = 4;

// r1 skips k correction; quo-est error is 2×TRUE_LIMBS (16). Start at 40; raise if bounds tooling reports overflow.
constexpr int BASE_REDUNDANT_MULTIPLE_NO_K = 40;

/** IntRNS2 scalar r1/r2 (generated moduli, no perm). Baseline for benches. */
using BlsRing = BoundedRing<
    BLS12_381_BITS, BLS12_381_TRUE_LIMBS, BLS12_381_ELEMENT_BITS,
    MIN_RNS, MAX_RNS, LOG_MULTIPLES, BLS12_381_3B, MAX_ADD, BLS12_381_6B,
    BASE_REDUNDANT_MULTIPLE_MATRIX, BLS12_381_PADDED_LIMBS, ChangeBaseVariation::Matrix>;

/** IntRNS4 export: r1 8×8 cyclic perm (no k), r2 scalar matrix, qr rotation. */
using BlsRingMatrixNoK = BoundedRing<
    BLS12_381_BITS, BLS12_381_TRUE_LIMBS, BLS12_381_ELEMENT_BITS,
    MIN_RNS, MAX_RNS, LOG_MULTIPLES, BLS12_381_3B, MAX_ADD, BLS12_381_6B,
    BASE_REDUNDANT_MULTIPLE_NO_K, BLS12_381_PADDED_LIMBS, ChangeBaseVariation::MatrixNoK>;

using BlsRingFixedPerm = BlsRingMatrixNoK;
using BlsRingProduction = BlsRingMatrixNoK;

// EC point ops need deeper reduction tables (matches bench_pairing_50bit G1/G2).
constexpr int LOG_MULTIPLES_EC = 12;

using BlsRingEC = BoundedRing<
    BLS12_381_BITS, BLS12_381_TRUE_LIMBS, BLS12_381_ELEMENT_BITS,
    MIN_RNS, MAX_RNS, LOG_MULTIPLES_EC, BLS12_381_3B, MAX_ADD, BLS12_381_6B,
    BASE_REDUNDANT_MULTIPLE_MATRIX, BLS12_381_PADDED_LIMBS, ChangeBaseVariation::Matrix>;

using BlsRingMatrixNoKEC = BoundedRing<
    BLS12_381_BITS, BLS12_381_TRUE_LIMBS, BLS12_381_ELEMENT_BITS,
    MIN_RNS, MAX_RNS, LOG_MULTIPLES_EC, BLS12_381_3B, MAX_ADD, BLS12_381_6B,
    BASE_REDUNDANT_MULTIPLE_NO_K, BLS12_381_PADDED_LIMBS, ChangeBaseVariation::MatrixNoK>;

} // namespace bls12_381_params
