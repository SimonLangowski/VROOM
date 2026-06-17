// Sum of products with delayed reduction: (a*b) + (c*d) before a single reduce step.
//
// operator* defers multiply; complete() accumulates wide limbs;
// batch_reduce_expand() runs ready() + reduction in one call.
//
// Build: make -C examples 03_sum_of_products
// Run:   ./examples/03_sum_of_products

#include "../src/bls12_381_ring_params.hpp"
#include "../cpu/precompute/gmp_wrapper.hpp"
#include <array>
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
    BigInt expected = (a_int * b_int + c_int * d_int) % p;

    Ring ring(p);

    std::cout << "Sum of products with delayed reduction (BLS12-381 field, MatrixNoK ring)\n\n";
    std::cout << "p = " << p.to_string() << "\n";
    std::cout << "a = " << a_int.to_string() << "\n";
    std::cout << "b = " << b_int.to_string() << "\n";
    std::cout << "c = " << c_int.to_string() << "\n";
    std::cout << "d = " << d_int.to_string() << "\n\n";

    auto a_rns = ring.from_bigint(a_int);
    auto b_rns = ring.from_bigint(b_int);
    auto c_rns = ring.from_bigint(c_int);
    auto d_rns = ring.from_bigint(d_int);

    // Delayed products (Karatsuba-style: prep_left defers reduction on the left factor)
    auto ab = a_rns * b_rns;
    auto cd = c_rns * d_rns;

    // Materialise each product, add wide results, single reduction pass
    auto sum_wide = ab + cd;
    auto [result_rns] = ring.batch_reduce_expand(sum_wide);

    BigInt got = ring.to_bigint(result_rns);
    std::cout << "(a*b + c*d) mod p (decoded) = " << got.to_string() << "\n";
    std::cout << "expected                 = " << expected.to_string() << "\n";

    if (got != expected) {
        std::cerr << "FAIL: sum of products mismatch\n";
        return 1;
    }
}
