// Difference of products with delayed reduction: (a*b) - (c*d) before a single reduce step.
//
// Subtraction of delayed products (e.g. ab - cd) is not supported. Negate one
// multiplicand before multiply — here ring.negate(c) — then accumulate with the
// fused add path used for sums of products.
//
// Build: make -C examples 04_difference_of_products
// Run:   ./examples/04_difference_of_products

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
    BigInt expected = (a_int * b_int - c_int * d_int) % p;

    Ring ring(p);

    std::cout << "Difference of products with delayed reduction (BLS12-381 field, MatrixNoK ring)\n\n";
    std::cout << "p = " << p.to_string() << "\n";
    std::cout << "a = " << a_int.to_string() << "\n";
    std::cout << "b = " << b_int.to_string() << "\n";
    std::cout << "c = " << c_int.to_string() << "\n";
    std::cout << "d = " << d_int.to_string() << "\n\n";

    auto a_rns = ring.from_bigint(a_int);
    auto b_rns = ring.from_bigint(b_int);
    auto c_rns = ring.from_bigint(c_int);
    auto d_rns = ring.from_bigint(d_int);

    auto ab = a_rns * b_rns;
    // Not supported: auto diff_wide = ab - (c_rns * d_rns);
    // Rewrite (a*b) - (c*d) as (a*b) + (-c)*d by negating c before multiply.
    auto neg_c_d = ring.negate(c_rns) * d_rns;

    auto diff_wide = ab + neg_c_d;
    auto [result_rns] = ring.batch_reduce_expand(diff_wide);

    BigInt got = ring.to_bigint(result_rns);
    std::cout << "(a*b - c*d) mod p (decoded) = " << got.to_string() << "\n";
    std::cout << "expected                  = " << expected.to_string() << "\n";

    if (got != expected) {
        std::cerr << "FAIL: difference of products mismatch\n";
        return 1;
    }

    std::cout << "\nOK — negate-before-multiply difference matches GMP reference value.\n";
    return 0;
}
