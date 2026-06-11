#include "bls12_381_ring_params.hpp"
#include "../cpu/precompute/gmp_wrapper.hpp"
#include <array>
#include <iostream>

namespace test_values_bls12_381 {
#include "test_data/test_values_bls12_381.hpp"
}

using BlsRing = bls12_381_params::BlsRingProduction;

int main() {
    BigInt p(test_values_bls12_381::modulus);
    BigInt a(test_values_bls12_381::a);
    BigInt b(test_values_bls12_381::b);
    BigInt c_exp(test_values_bls12_381::c_correct);

    BlsRing ring(p);
    auto ra = ring.from_bigint(a);
    auto rb = ring.from_bigint(b);
    std::array<BlsRing::StandardElement, 1> as = {ra};
    std::array<BlsRing::StandardElement, 1> bs = {rb};
    BigInt c_got = ring.to_bigint(ring.batch_modmul<1>(as, bs)[0]);

    if (c_got != c_exp) {
        std::cerr << "BLS12-381 FixedPerm field mul FAIL\n  got " << c_got.to_string()
                  << "\n  exp " << c_exp.to_string() << "\n";
        return 1;
    }
    std::cout << "BLS12-381 FixedPerm field mul OK (a*b mod p)" << std::endl;

    // Stress: repeated modmul (bounds / skip-k accumulation)
    auto acc = ra;
    BigInt acc_int = a;
    constexpr int kStress = 10;
    for (int i = 0; i < kStress; i++) {
        std::array<BlsRing::StandardElement, 1> acc_arr = {acc};
        std::array<BlsRing::StandardElement, 1> b_arr = {rb};
        acc = ring.batch_modmul<1>(acc_arr, b_arr)[0];
        acc_int = (acc_int * b) % p;
        if (ring.to_bigint(acc) != acc_int) {
            std::cerr << "BLS12-381 stress modmul FAIL at i=" << i << "\n";
            return 1;
        }
    }
    std::cout << "BLS12-381 stress modmul OK (" << kStress << " muls)" << std::endl;
    return 0;
}
