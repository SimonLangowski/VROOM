// Demonstrates RNS parameter choices: modulus bit width, residue count,
// limb bit width, and change-base strategy (Matrix vs MatrixNoK vs FixedPerm).
//
// Build: make -C examples 01_parameter_setup
// Run:   ./examples/01_parameter_setup

#include "../src/bls12_381_ring_params.hpp"
#include "../src/bn254_ring_params.hpp"
#include "../cpu/precompute/gmp_wrapper.hpp"
#include <cstdint>
#include <iostream>

template <typename Ring>
void print_moduli(const char *label, const std::array<uint64_t, Ring::LIMBS> &moduli, int count) {
    std::cout << "  " << label << " (" << count << " active): ";
    for (int i = 0; i < count; ++i) {
        if (i) std::cout << ", ";
        std::cout << moduli[i];
    }
    std::cout << "\n";
}

template <typename Ring>
void describe_ring(const char *name, const char *modulus_decimal, int modulus_bits, int element_bits,
                   const char *change_base) {
    BigInt p(modulus_decimal);
    Ring ring(p);

    std::cout << "\n=== " << name << " ===\n";
    std::cout << "  Target modulus bits: " << modulus_bits << "\n";
    std::cout << "  True RNS limbs:      " << Ring::TRUE_LIMBS << "\n";
    std::cout << "  Padded AVX limbs:    " << Ring::LIMBS << "\n";
    std::cout << "  Element (residue) bits: " << element_bits << "\n";
    std::cout << "  RNS capacity bits:   " << Ring::TRUE_LIMBS * element_bits - 1 << "\n";
    std::cout << "  Change-base variant: " << change_base << "\n";

    print_moduli<Ring>("moduli1 (m1)", ring.sys.moduli1, Ring::TRUE_LIMBS);
    print_moduli<Ring>("moduli2 (m2)", ring.sys.moduli2, Ring::TRUE_LIMBS);
}

int main() {
    std::cout << "VROOM RNS parameter configurations\n";
    std::cout << "Each ring uses two residue sets (m1, m2) for IntRNS2/IntRNS4 Montgomery multiplication.\n";

    describe_ring<bls12_381_params::BlsRing>(
        "BLS12-381 / Matrix (baseline IntRNS2)",
        bls12_381_params::BLS12_381_P,
        bls12_381_params::BLS12_381_BITS,
        bls12_381_params::BLS12_381_ELEMENT_BITS,
        "Matrix — scalar r1/r2 change-base matrices");

    describe_ring<bls12_381_params::BlsRingMatrixNoK>(
        "BLS12-381 / MatrixNoK (production, 50-bit)",
        bls12_381_params::BLS12_381_P,
        bls12_381_params::BLS12_381_BITS,
        bls12_381_params::BLS12_381_ELEMENT_BITS,
        "MatrixNoK — r1 cyclic permutation, r2 scalar matrix, qr rotation");

    describe_ring<bn254_params::Bn254Ring>(
        "BN254 / Matrix (6×46-bit baseline)",
        bn254_params::BN254_P,
        bn254_params::BN254_BITS,
        bn254_params::BN254_ELEMENT_BITS,
        "Matrix — 6 true limbs, 8 padded for AVX-512");

    describe_ring<bn254_params::Bn254RingFixedPerm>(
        "BN254 / FixedPerm (production)",
        bn254_params::BN254_P,
        bn254_params::BN254_BITS,
        bn254_params::BN254_ELEMENT_BITS,
        "FixedPerm — precomputed r1 permutation from Python export");

    std::cout << "\nModuli/residue combinations are chosen in scripts/rns_secp256k1.py\n";
    std::cout << "and exported to scripts/*_intrns4_params.hpp for FixedPerm/MatrixNoK rings.\n";
    return 0;
}
