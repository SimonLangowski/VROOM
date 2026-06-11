#include "../cpu/precompute/gmp_wrapper.hpp"
extern "C" {
#include "../blst/vect.h"
#include "../blst/fields.h"
#include "../blst/consts.h"
#include "../blst/point.h"
}
#include "pairing.hpp"
#include "bls12_381_ring_params.hpp"
#include "fp2.hpp"
#include <iostream>
#include <random>
#include <cstring>

using BlsRing = bls12_381_params::BlsRingMatrixNoK;

const char *bls12_381_modulus_hex =
    "1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab";

extern "C" {
extern const POINTonE1 BLS12_381_G1;
extern const POINTonE2 BLS12_381_G2;
void blst_p1_mult(POINTonE1 *out, const POINTonE1 *a, const byte *scalar, size_t nbits);
void blst_p1_to_affine(POINTonE1_affine *out, const POINTonE1 *a);
void blst_p2_mult(POINTonE2 *out, const POINTonE2 *a, const byte *scalar, size_t nbits);
void blst_p2_to_affine(POINTonE2_affine *out, const POINTonE2 *a);
void blst_miller_loop(vec384fp12 ret, const POINTonE2_affine *Q, const POINTonE1_affine *P);
void blst_final_exp(vec384fp12 ret, const vec384fp12 f);
}

static BigInt vec384_montgomery_to_bigint_local(const vec384 a) {
    vec384 normal;
    from_fp(normal, a);
    BigInt result(0);
    BigInt two_to_64 = BigInt(1) << 64;
    for (int i = 5; i >= 0; i--) {
        result = result * two_to_64 + BigInt(static_cast<unsigned long>(normal[i]));
    }
    return result;
}

static std::pair<BigInt, BigInt> vec384x_montgomery_to_bigint(const vec384x a) {
    return {vec384_montgomery_to_bigint_local(a[0]), vec384_montgomery_to_bigint_local(a[1])};
}

static std::array<BigInt, 12> vec384fp12_montgomery_to_bigint(const vec384fp12 a) {
    std::array<BigInt, 12> result;
    for (int i = 0; i < 12; i++) {
        result[i] = vec384_montgomery_to_bigint_local(a[i / 6][(i % 6) / 2][i % 2]);
    }
    return result;
}

static AffinePoint<typename BlsRing::StandardElement> blst_g1_to_affine_point(
    const POINTonE1_affine &p, const BlsRing &ring) {
    return {ring.from_bigint(vec384_montgomery_to_bigint_local(p.X)),
            ring.from_bigint(vec384_montgomery_to_bigint_local(p.Y))};
}

static AffinePoint<typename FP2Ring<BlsRing>::StandardElement> blst_g2_to_affine_point(
    const POINTonE2_affine &p, const BlsRing &ring) {
    auto x_pair = vec384x_montgomery_to_bigint(p.X);
    auto y_pair = vec384x_montgomery_to_bigint(p.Y);
    return {typename FP2Ring<BlsRing>::StandardElement(
                ring.from_bigint(x_pair.first), ring.from_bigint(x_pair.second)),
            typename FP2Ring<BlsRing>::StandardElement(
                ring.from_bigint(y_pair.first), ring.from_bigint(y_pair.second))};
}

static bool compare_fp12_results(const std::array<BigInt, 12> &blst_result,
    const typename FP12Ring<BlsRing>::StandardElement &our_result,
    const BlsRing &ring, const std::string &test_name) {
    bool all_match = true;
    BigInt q(bls12_381_modulus_hex, 16);
    for (int i = 0; i < 12; i++) {
        BigInt our_mod = ring.to_bigint(our_result[i]) % q;
        if (our_mod < 0) our_mod = our_mod + q;
        BigInt blst_mod = blst_result[i] % q;
        if (blst_mod < 0) blst_mod = blst_mod + q;
        if (our_mod != blst_mod) {
            all_match = false;
            std::cout << "Mismatch at index " << i << " in " << test_name << "\n"
                      << "  BLST: " << blst_mod.to_string(16) << "\n"
                      << "  Ours: " << our_mod.to_string(16) << "\n";
        }
    }
    return all_match;
}

int main() {
    try {
        std::cout << "=== BLS12-381 Pairing (MatrixNoK, 8x50-bit): BLST vs Our Implementation ===\n";

        BigInt q(bls12_381_modulus_hex, 16);
        BlsRing ring(q);
        FP2Ring<BlsRing> fp2ring(ring);
        FP12Ring<BlsRing> fp12ring(ring);

        int passed = 0, failed = 0;
        const int num_tests = 5;
        std::random_device rd;
        std::mt19937_64 gen(rd());

        for (int t = 0; t < num_tests; t++) {
            byte scalar1[32], scalar2[32];
            for (int i = 0; i < 32; i++) {
                scalar1[i] = static_cast<byte>(gen() & 0xff);
                scalar2[i] = static_cast<byte>(gen() & 0xff);
            }

            POINTonE1 p1_proj;
            blst_p1_mult(&p1_proj, &BLS12_381_G1, scalar1, 256);
            POINTonE1_affine p1_affine;
            blst_p1_to_affine(&p1_affine, &p1_proj);

            POINTonE2 p2_proj;
            blst_p2_mult(&p2_proj, &BLS12_381_G2, scalar2, 256);
            POINTonE2_affine p2_affine;
            blst_p2_to_affine(&p2_affine, &p2_proj);

            vec384fp12 blst_miller_result, blst_final_result;
            blst_miller_loop(blst_miller_result, &p2_affine, &p1_affine);
            blst_final_exp(blst_final_result, blst_miller_result);
            auto blst_result = vec384fp12_montgomery_to_bigint(blst_final_result);

            auto P = blst_g1_to_affine_point(p1_affine, ring);
            auto Q = blst_g2_to_affine_point(p2_affine, ring);
            auto our_result = pairing_with_blst_inverter(P, Q, fp2ring, fp12ring, ring);

            if (compare_fp12_results(blst_result, our_result, ring, "Test " + std::to_string(t + 1))) {
                passed++;
                std::cout << "Test " << (t + 1) << " passed\n";
            } else {
                failed++;
                std::cout << "Test " << (t + 1) << " FAILED\n";
            }
        }

        std::cout << "\n=== Summary (MatrixNoK) ===\nPassed: " << passed << "  Failed: " << failed << "\n";
        if (failed == 0) {
            std::cout << "All BLS12-381 pairing tests passed!\n";
            return 0;
        }
        std::cout << "Some pairing tests failed!\n";
        return 1;
    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
