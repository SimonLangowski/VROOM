// Singular modular multiplication: integer ↔ RNS conversion and a*b mod p.
//
// Build: make -C examples 02_singular_modmul
// Run:   ./examples/02_singular_modmul

#include "../src/bls12_381_ring_params.hpp"
#include "../cpu/precompute/gmp_wrapper.hpp"
#include <iostream>

namespace test_values_bls12_381 {
#include "../test_data/test_values_bls12_381.hpp"
}

using Ring = bls12_381_params::BlsRingProduction;

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

    auto a_rns = ring.from_bigint(a_int);
    auto b_rns = ring.from_bigint(b_int);
    auto c_rns = ring.modmul(a_rns, b_rns);
    BigInt c_got = ring.to_bigint(c_rns);

    std::cout << "a*b mod p (decoded) = " << c_got.to_string() << "\n";
    std::cout << "expected            = " << c_expected.to_string() << "\n";

    if (c_got != c_expected) {
        std::cerr << "FAIL: modular multiplication mismatch\n";
        return 1;
    }

    std::cout << "\nOK — singular modular multiplication matches GMP reference value.\n";
    return 0;
}
