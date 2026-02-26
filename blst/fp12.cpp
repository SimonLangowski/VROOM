#include <iostream>
#include <array>
#include <random>
#include <cstdint>
#include <cstring>
#include <cassert>
#include <string>
extern "C" {
#include "fields.h"
}
#include "fp12.hpp"

// Global offset of type FP with value 0 for flattened_code.hpp (additive identity)
FP offset;

#include "flattened_code.hpp"


// Generate random 384-bit field element in Montgomery form
void random_fp(vec384 out) {
    std::random_device rd;
    std::mt19937_64 gen(rd());
    
    // Generate random value less than P
    vec384 temp;
    do {
        for (int i = 0; i < 6; i++) {
            temp[i] = gen();
        }
        // Make sure it's less than P by checking high bits
        // P[5] = 0x1a0111ea397fe69a, so generate values with high bits that are likely < P
        if (temp[5] > BLS12_381_P[5]) {
            temp[5] = gen() % (BLS12_381_P[5] + 1);
        }
    } while (temp[5] == BLS12_381_P[5] && temp[4] >= BLS12_381_P[4]); // Simple check
    
    // Reduce modulo P to ensure it's valid
    vec384 reduced;
    // Check if temp >= P and reduce if needed
    bool ge = false;
    for (int i = 5; i >= 0; i--) {
        if (temp[i] > BLS12_381_P[i]) {
            ge = true;
            break;
        } else if (temp[i] < BLS12_381_P[i]) {
            break;
        }
    }
    if (ge) {
        sub_mod_384(reduced, temp, BLS12_381_P, BLS12_381_P);
    } else {
        vec_copy(reduced, temp, sizeof(vec384));
    }
    // Convert to Montgomery form
    mul_fp(out, reduced, BLS12_381_RR);
}

// Helper function to compare two vec384fp12 results
bool compare_results(const vec384fp12 &standard, const std::array<FP, 12> &flat, const std::string &test_name) {
    bool all_match = true;
    int match_count = 0;
    int mismatch_count = 0;
    
    for (int i = 0; i < 12; i++) {
        vec384 standard_result;
        vec_copy(standard_result, standard[i/6][(i%6)/2][i%2], sizeof(vec384));
        
        vec384 flat_result;
        flat[i].to_vec384(flat_result);
        
        bool match = true;
        for (int j = 0; j < 6; j++) {
            if (standard_result[j] != flat_result[j]) {
                match = false;
                break;
            }
        }
        
        if (match) {
            match_count++;
            std::cout << "Element " << i << " (index " << std::dec << i << std::hex << "): MATCH" << std::endl;
        } else {
            mismatch_count++;
            std::cout << "Mismatch at output element " << std::dec << i << " (index " << i << "):" << std::endl;
            std::cout << "  Standard: ";
            for (int j = 5; j >= 0; j--) {
                std::cout << std::hex << standard_result[j] << " ";
            }
            std::cout << std::endl;
            std::cout << "  Flat:     ";
            for (int j = 5; j >= 0; j--) {
                std::cout << std::hex << flat_result[j] << " ";
            }
            std::cout << std::endl;
            all_match = false;
        }
    }
    std::cout << std::dec << test_name << " Summary: " << match_count << " match, " << mismatch_count << " mismatch" << std::endl;
    
    if (all_match) {
        std::cout << "SUCCESS: All 12 output elements match!" << std::endl;
    } else {
        std::cout << "FAILURE: Results do not match!" << std::endl;
    }
    std::cout << std::endl;
    
    return all_match;
}

bool test_mul_fp12() {
    std::cout << "=== Testing mul_fp12 vs flat_fp12_mul with 24 random field elements ===" << std::endl;
    
    // Generate 24 random field elements (12 for x, 12 for y)
    std::array<vec384, 12> x_mont, y_mont;
    for (int i = 0; i < 12; i++) {
        random_fp(x_mont[i]);
        random_fp(y_mont[i]);
    }
    
    // Convert to FP format for flat_fp12_mul
    std::array<FP, 12> x_fp, y_fp;
    for (int i = 0; i < 12; i++) {
        x_fp[i] = FP(x_mont[i]);
        y_fp[i] = FP(y_mont[i]);
    }
    
    // Convert to vec384fp12 format for mul_fp12_export
    vec384fp12 a, b, c_standard;
    for (int i = 0; i < 12; i++) {
        vec_copy(a[i/6][(i%6)/2][i%2], x_mont[i], sizeof(vec384));
        vec_copy(b[i/6][(i%6)/2][i%2], y_mont[i], sizeof(vec384));
    }
    
    // Compute using standard mul_fp12
    mul_fp12_export(c_standard, a, b);
    
    // Compute using flat_fp12_mul
    RingFP ring;
    std::array<FP, 12> c_flat = fp12_flat_mul(x_fp, y_fp, ring);
    
    return compare_results(c_standard, c_flat, "mul_fp12");
}

bool test_sqr_fp12() {
    std::cout << "=== Testing sqr_fp12 vs flat_fp12_sqr with 12 random field elements ===" << std::endl;
    
    // Generate 12 random field elements for x
    std::array<vec384, 12> x_mont;
    for (int i = 0; i < 12; i++) {
        random_fp(x_mont[i]);
    }
    
    // Convert to FP format for flat_fp12_sqr
    std::array<FP, 12> x_fp;
    for (int i = 0; i < 12; i++) {
        x_fp[i] = FP(x_mont[i]);
    }
    
    // Convert to vec384fp12 format for sqr_fp12_export
    vec384fp12 a, c_standard;
    for (int i = 0; i < 12; i++) {
        vec_copy(a[i/6][(i%6)/2][i%2], x_mont[i], sizeof(vec384));
    }
    
    // Compute using standard sqr_fp12
    sqr_fp12_export(c_standard, a);
    
    // Compute using flat_fp12sqr
    RingFP ring;
    std::array<FP, 12> c_flat = fp12_flat_sqr(x_fp, ring);
    
    return compare_results(c_standard, c_flat, "sqr_fp12");
}

bool test_mul_by_xy00z0() {
    std::cout << "=== Testing mul_by_xy00z0_fp12 vs flat_mulx0y0z00 with 12+6 random field elements ===" << std::endl;
    
    // Generate 12 random field elements for x
    std::array<vec384, 12> x_mont;
    for (int i = 0; i < 12; i++) {
        random_fp(x_mont[i]);
    }
    
    // Generate 6 random field elements for y (xy00z0 pattern)
    std::array<vec384, 6> y_mont;
    for (int i = 0; i < 6; i++) {
        random_fp(y_mont[i]);
    }
    
    // Convert to FP format for flat_mulx0y0z00
    std::array<FP, 12> x_fp;
    for (int i = 0; i < 12; i++) {
        x_fp[i] = FP(x_mont[i]);
    }
    std::array<FP, 6> y_fp;
    for (int i = 0; i < 6; i++) {
        y_fp[i] = FP(y_mont[i]);
    }
    
    // Convert to vec384fp12 format for mul_by_xy00z0_fp12_export
    vec384fp12 a, c_standard;
    for (int i = 0; i < 12; i++) {
        vec_copy(a[i/6][(i%6)/2][i%2], x_mont[i], sizeof(vec384));
    }
    
    // Convert to vec384fp6 format for mul_by_xy00z0_fp12_export
    vec384fp6 xy00z0;
    for (int i = 0; i < 6; i++) {
        vec_copy(xy00z0[i/2][i%2], y_mont[i], sizeof(vec384));
    }
    
    // Compute using standard mul_by_xy00z0_fp12
    mul_by_xy00z0_fp12_export(c_standard, a, xy00z0);
    
    // Compute using flat_mulxy0z00
    RingFP ring;
    std::array<FP, 12> c_flat = flat_mulxy0z00(x_fp, y_fp, ring);
    
    return compare_results(c_standard, c_flat, "mul_by_xy00z0_fp12");
}

bool test_cyclotomic_sqr_fp12() {
    std::cout << "=== Testing cyclotomic_sqr_fp12 vs flat_cyclotomic_sqr_fp12 with 12 random field elements ===" << std::endl;
    
    // Generate 12 random field elements for x
    std::array<vec384, 12> x_mont;
    for (int i = 0; i < 12; i++) {
        random_fp(x_mont[i]);
    }
    
    // Convert to FP format for flat_cyclotomic_sqr_fp12
    std::array<FP, 12> x_fp;
    for (int i = 0; i < 12; i++) {
        x_fp[i] = FP(x_mont[i]);
    }
    
    // Convert to vec384fp12 format for cyclotomic_sqr_fp12_export
    vec384fp12 a, c_standard;
    for (int i = 0; i < 12; i++) {
        vec_copy(a[i/6][(i%6)/2][i%2], x_mont[i], sizeof(vec384));
    }
    
    // Compute using standard cyclotomic_sqr_fp12
    cyclotomic_sqr_fp12_export(c_standard, a);
    
    // Compute using flat_cyclotomic_sqr_fp12
    RingFP ring;
    std::array<FP, 12> c_flat = fp12_flat_cyclotomic_sqr(x_fp, ring);
    
    return compare_results(c_standard, c_flat, "cyclotomic_sqr_fp12");
}

bool test_frobenius_map_fp12(size_t n) {
    std::cout << "=== Testing frobenius_map_fp12 (n=" << n << ") vs flatish_fp12_frobenius with 12 random field elements ===" << std::endl;
    
    // Generate 12 random field elements for x
    std::array<vec384, 12> x_mont;
    for (int i = 0; i < 12; i++) {
        random_fp(x_mont[i]);
    }
    
    // Convert to FP format for flatish_fp12_frobenius
    std::array<FP, 12> x_fp;
    for (int i = 0; i < 12; i++) {
        x_fp[i] = FP(x_mont[i]);
    }
    
    // Convert to vec384fp12 format for frobenius_map_fp12_export
    vec384fp12 a, c_standard;
    for (int i = 0; i < 12; i++) {
        vec_copy(a[i/6][(i%6)/2][i%2], x_mont[i], sizeof(vec384));
    }
    
    // Compute using standard frobenius_map_fp12
    frobenius_map_fp12_export(c_standard, a, n);
    
    // Compute using flat_frob1, flat_frob2, or flat_frob3
    RingFP ring;
    std::array<FP, 12> c_flat;
    if (n == 1) {
        c_flat = fp12_flat_frobenius_map1(ring, x_fp);
    } else if (n == 2) {
        c_flat = fp12_flat_frobenius_map2(ring, x_fp);
    } else if (n == 3) {
        c_flat = fp12_flat_frobenius_map3(ring, x_fp);
    } else {
        std::cout << "ERROR: Invalid n value " << n << std::endl;
        return false;
    }
    
    std::string test_name = "frobenius_map_fp12_n" + std::to_string(n);
    return compare_results(c_standard, c_flat, test_name);
}

int main() {
    bool all_passed = true;
    
    all_passed &= test_mul_fp12();
    all_passed &= test_sqr_fp12();
    all_passed &= test_mul_by_xy00z0();
    all_passed &= test_cyclotomic_sqr_fp12();
    all_passed &= test_frobenius_map_fp12(1);
    all_passed &= test_frobenius_map_fp12(2);
    all_passed &= test_frobenius_map_fp12(3);
    
    if (all_passed) {
        std::cout << "=== ALL TESTS PASSED ===" << std::endl;
        return 0;
    } else {
        std::cout << "=== SOME TESTS FAILED ===" << std::endl;
        return 1;
    }
}

