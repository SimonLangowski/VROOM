// Singular modular multiplication: integer ↔ RNS conversion and a*b mod p.
//
// Build: make -C examples 02_singular_modmul
// Run:   ./examples/02_singular_modmul

#include "../src/bls12_381_ring_params.hpp"
#include "../cpu/precompute/gmp_wrapper.hpp"
#include <array>
#include <iostream>

namespace test_values_bls12_381 {
#include "../test_data/test_values_bls12_381.hpp"
}

using Ring = bls12_381_params::BlsRingProduction;

static void print_residues(const char *label, const std::array<int64_t, Ring::LIMBS> &m1,
                           const std::array<int64_t, Ring::LIMBS> &m2, int limbs) {
    std::cout << label << " m1=[";
    for (int i = 0; i < limbs; ++i) {
        if (i) std::cout << ", ";
        std::cout << static_cast<uint64_t>(m1[i]);
    }
    std::cout << "]\n       m2=[";
    for (int i = 0; i < limbs; ++i) {
        if (i) std::cout << ", ";
        std::cout << static_cast<uint64_t>(m2[i]);
    }
    std::cout << "]\n";
}

int main() {
    BigInt p(test_values_bls12_381::modulus);
    BigInt a_int(test_values_bls12_381::a);
    BigInt b_int(test_values_bls12_381::b);
    BigInt c_expected(test_values_bls12_381::c_correct);

    Ring ring(p);

    std::cout << "Singular modular multiplication (BLS12-381 field, MatrixNoK ring)\n\n";
    std::cout << "p = " << p.to_string() << "\n";
    std::cout << "a = " << a_int.to_string() << "\n";
    std::cout << "b = " << b_int.to_string() << "\n\n";

    // --- Integer → RNS (convert_to path) ---
    auto a_rns = ring.from_bigint(a_int);
    auto b_rns = ring.from_bigint(b_int);
    print_residues("a in RNS", a_rns.m1.to_array(), a_rns.m2.to_array(), Ring::TRUE_LIMBS);
    print_residues("b in RNS", b_rns.m1.to_array(), b_rns.m2.to_array(), Ring::TRUE_LIMBS);

    // --- Montgomery modmul in RNS ---
    auto c_rns = ring.modmul(a_rns, b_rns);
    print_residues("a*b in RNS (before decode)", c_rns.m1.to_array(), c_rns.m2.to_array(), Ring::TRUE_LIMBS);

    // --- RNS → Integer (convert_from / from_mont path) ---
    BigInt c_got = ring.to_bigint(c_rns);
    std::cout << "\na*b mod p (decoded) = " << c_got.to_string() << "\n";
    std::cout << "expected            = " << c_expected.to_string() << "\n";

    if (c_got != c_expected) {
        std::cerr << "FAIL: modular multiplication mismatch\n";
        return 1;
    }

    // Cross-check with exact Montgomery encode/decode (bypasses fast convert path)
    auto a_exact = ring.from_bigint_exact(a_int);
    auto b_exact = ring.from_bigint_exact(b_int);
    BigInt c_exact = ring.to_bigint_exact(ring.modmul(a_exact, b_exact));
    if (c_exact != c_expected) {
        std::cerr << "FAIL: exact Montgomery path mismatch\n";
        return 1;
    }

    std::cout << "\nOK — integer→RNS→modmul→integer matches GMP reference value.\n";
    return 0;
}
