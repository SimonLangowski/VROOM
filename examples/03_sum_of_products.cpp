// Sum of products with delayed reduction: (a*b) + (c*d) before a single reduce step.
//
// Multiplication via operator* returns a DelayedProduct; .complete() materialises
// each product; addition accumulates wide limbs; reduce_auto finishes Montgomery
// reduction per residue lane.
//
// Build: make -C examples 03_sum_of_products
// Run:   ./examples/03_sum_of_products

#include "../src/elementwise_reduce.hpp"
#include "../src/ring_element.hpp"
#include "../src/bounded_element.hpp"
#include "../cpu/vector/vector_impl.hpp"
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>

using int128_t = __int128_t;

template <int limbs, int element_bits>
void demo_sum_of_products() {
    std::cout << "\n--- limbs=" << limbs << ", element_bits=" << element_bits << " ---\n";

    // Small primes t_i; moduli p_i = 2^element_bits - t_i (same pattern as test_elementwise_reduce)
    std::array<uint64_t, limbs> t_values = {1, 3, 5, 7, 11, 13, 17, 19};
    std::array<uint64_t, limbs> moduli{};
    for (int i = 0; i < limbs; ++i) {
        moduli[i] = (1ull << element_bits) - t_values[i];
    }

    MontgomeryReduce2<limbs, element_bits> reducer(moduli);

    // Fixed test values (deterministic, no RNG)
    std::array<uint64_t, limbs> arr_a = {100, 200, 300, 400, 500, 600, 700, 800};
    std::array<uint64_t, limbs> arr_b = {  7,  11,  13,  17,  19,  23,  29,  31};
    std::array<uint64_t, limbs> arr_c = {  3,   5,   7,  11,  13,  17,  19,  23};
    std::array<uint64_t, limbs> arr_d = { 41,  43,  47,  53,  59,  61,  67,  71};

    // Lift to Montgomery form per lane
    auto to_mont_array = [&](const std::array<uint64_t, limbs> &vals) {
        std::array<uint64_t, limbs> mont{};
        for (int i = 0; i < limbs; ++i) {
            std::array<uint64_t, limbs> single{};
            single[i] = vals[i];
            AVXVector<limbs> v;
            v.load(single.data());
            mont[i] = reducer.to_mont(v).data.to_array()[i];
        }
        return mont;
    };

    auto mont_a = to_mont_array(arr_a);
    auto mont_b = to_mont_array(arr_b);
    auto mont_c = to_mont_array(arr_c);
    auto mont_d = to_mont_array(arr_d);

    auto a = from_constant<limbs, element_bits>(mont_a, mont_a);
    auto b = from_constant<limbs, element_bits>(mont_b, mont_b);
    auto c = from_constant<limbs, element_bits>(mont_c, mont_c);
    auto d = from_constant<limbs, element_bits>(mont_d, mont_d);

    // Delayed products: a*b and c*d are not reduced yet
    auto ab = a * b;
    auto cd = c * d;

    // Materialise each product, then add wide results
    auto sum_products = ab.complete() + cd.complete();

    // Single reduction pass over the accumulated sum
    auto reduced = reducer.template reduce_auto<1>(sum_products.m1);
    auto normalized = reducer.min_sub_reduce(reduced);
    auto recovered = reducer.from_mont(normalized).to_array();

    for (int i = 0; i < limbs; ++i) {
        uint64_t p = moduli[i];
        uint64_t expected = (((uint128_t)arr_a[i] * arr_b[i]) + ((uint128_t)arr_c[i] * arr_d[i])) % p;
        if (recovered[i] != expected) {
            std::cerr << "FAIL lane " << i << ": got " << recovered[i] << ", expected " << expected << "\n";
            assert(false);
        }
    }

    std::cout << "  (a*b + c*d) mod p per lane — OK\n";
    std::cout << "  Example lane 0: (" << arr_a[0] << "*" << arr_b[0] << " + " << arr_c[0] << "*" << arr_d[0]
              << ") mod " << moduli[0] << " = " << recovered[0] << "\n";
}

int main() {
    std::cout << "Sum of products with delayed reduction\n";
    std::cout << "operator* defers multiply; complete() + reduce_auto() amortises reduction.\n";

    demo_sum_of_products<8, 52>();
    demo_sum_of_products<8, 50>();

    std::cout << "\nAll delayed-reduction examples passed.\n";
    return 0;
}
