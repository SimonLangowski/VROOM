// Product of sums: (a+b)*(c+d) with prep before delayed multiply.
//
// After addition each sum can exceed 52 bits per limb; ring.prep() elementwise-
// reduces so IFMA multiply is valid. ring.prep_left() is the same on the base
// field; in field extensions it prepares only the left operand for the
// multiplication layout.
//
// Build: make -C examples 05_product_of_sums
// Run:   ./examples/05_product_of_sums

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
    BigInt c_int = b_int;
    BigInt d_int = a_int;
    BigInt expected = (((a_int + b_int) % p) * ((c_int + d_int) % p)) % p;

    Ring ring(p);

    std::cout << "Product of sums (BLS12-381 field, MatrixNoK ring)\n\n";
    std::cout << "p = " << p.to_string() << "\n";
    std::cout << "a = " << a_int.to_string() << "\n";
    std::cout << "b = " << b_int.to_string() << "\n";
    std::cout << "c = " << c_int.to_string() << "\n";
    std::cout << "d = " << d_int.to_string() << "\n\n";

    auto a_rns = ring.from_bigint(a_int);
    auto b_rns = ring.from_bigint(b_int);
    auto c_rns = ring.from_bigint(c_int);
    auto d_rns = ring.from_bigint(d_int);

    // Sums can be wider than 52 bits; prep before multiply (prep_left == prep here).
    auto product_wide = ring.prep_left(a_rns + b_rns) * ring.prep(c_rns + d_rns);
    auto [result_rns] = ring.batch_reduce_expand(product_wide);

    BigInt got = ring.to_bigint(result_rns);
    std::cout << "(a+b)*(c+d) mod p (decoded) = " << got.to_string() << "\n";
    std::cout << "expected                   = " << expected.to_string() << "\n";

    if (got != expected) {
        std::cerr << "FAIL: product of sums mismatch\n";
        return 1;
    }

    std::cout << "\nOK — prep-before-multiply product of sums matches GMP reference value.\n";
    return 0;
}
