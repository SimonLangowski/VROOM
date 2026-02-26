#include <iostream>
#include <array>
#include <random>
#include <cstdint>
#include <cstring>
#include <cassert>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <cctype>
extern "C" {
#include "fields.h"
// Forward declarations for blst functions
void blst_final_exp(vec384fp12 ret, const vec384fp12 f);
void blst_fp12_inverse(vec384fp12 ret, const vec384fp12 a);
void blst_fp12_conjugate(vec384fp12 a);
void blst_fp12_mul(vec384fp12 ret, const vec384fp12 a, const vec384fp12 b);
void blst_fp12_frobenius_map(vec384fp12 ret, const vec384fp12 a, size_t n);
void blst_fp12_cyclotomic_sqr(vec384fp12 ret, const vec384fp12 a);
}
#include "fp12.hpp"

// Only need wrapper for fp12_flat_inverse (not used by this test; kept for API consistency)
template<class Ring>
inline auto fp12_flat_inverse(const std::array<typename Ring::StandardElement, 12> &x, Ring &ring) {
    // Convert to vec384fp12
    vec384fp12 a, result;
    for (int i = 0; i < 12; i++) {
        x[i].to_vec384(a[i/6][(i%6)/2][i%2]);
    }
    
    // Compute inverse using blst
    blst_fp12_inverse(result, a);
    
    // Convert back to array - result[i/6][(i%6)/2][i%2] is a vec384&
    std::array<typename Ring::StandardElement, 12> ret;
    for (int i = 0; i < 12; i++) {
        ret[i] = typename Ring::StandardElement(result[i/6][(i%6)/2][i%2]);
    }
    return ret;
}

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

// Helper function to check if an fp12 element is zero
bool is_zero_fp12(const std::array<FP, 12> &x) {
    for (int i = 0; i < 12; i++) {
        for (int j = 0; j < 6; j++) {
            if (x[i].d[j] != 0) return false;
        }
    }
    return true;
}

// Helper function to print fp12 element (first element only for brevity)
void print_fp12_element(const std::string &name, const std::array<FP, 12> &x) {
    std::cout << name << "[0] = ";
    for (int j = 5; j >= 0; j--) {
        std::cout << std::hex << x[0].d[j] << " ";
    }
    std::cout << std::endl;
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

// Helper to convert vec384fp12 to array<FP, 12>
std::array<FP, 12> vec384fp12_to_array(const vec384fp12 &x) {
    std::array<FP, 12> result;
    for (int i = 0; i < 12; i++) {
        result[i] = FP(x[i/6][(i%6)/2][i%2]);
    }
    return result;
}

// Helper to print first element of fp12 for comparison
void print_fp12_first(const std::string &name, const std::array<FP, 12> &x) {
    std::cout << name << "[0] = ";
    for (int j = 5; j >= 0; j--) {
        std::cout << std::hex << x[0].d[j] << " ";
    }
    std::cout << std::endl;
}

void print_fp12_first(const std::string &name, const vec384fp12 &x) {
    std::cout << name << "[0] = ";
    vec384 first;
    vec_copy(first, x[0][0][0], sizeof(vec384));
    for (int j = 5; j >= 0; j--) {
        std::cout << std::hex << first[j] << " ";
    }
    std::cout << std::endl;
}

// Helper function to convert vec384 from Montgomery to standard form and format as hex string
std::string vec384_to_hex_string(const vec384 &mont_value) {
    vec384 standard;
    from_fp(standard, mont_value);
    
    std::ostringstream oss;
    oss << "0x";
    // Print from most significant to least significant
    for (int j = 5; j >= 0; j--) {
        oss << std::hex << std::setfill('0') << std::setw(16) << standard[j];
    }
    return oss.str();
}

// Helper to convert FP to hex string via vec384
std::string fp_to_hex_string(const FP &fp) {
    vec384 v;
    fp.to_vec384(v);
    return vec384_to_hex_string(v);
}

bool test_final_exp() {
    std::cout << "=== Generating pairing values using blst functions ===" << std::endl;
    
    // Generate random fp12 element
    std::array<vec384, 12> f_mont;
    for (int i = 0; i < 12; i++) {
        random_fp(f_mont[i]);
    }
    
    // Convert to vec384fp12 format for blst functions
    vec384fp12 f_standard;
    for (int i = 0; i < 12; i++) {
        vec_copy(f_standard[i/6][(i%6)/2][i%2], f_mont[i], sizeof(vec384));
    }
    
    // Compute inverse using blst_fp12_inverse
    vec384fp12 f_inverse_standard;
    blst_fp12_inverse(f_inverse_standard, f_standard);
    
    // Compute final exponentiation using blst_final_exp
    vec384fp12 result_standard;
    blst_final_exp(result_standard, f_standard);
    
    // Write values to pairing_values.hpp using blst function results
    std::ofstream outfile("pairing_values.hpp");
    outfile << "#pragma once\n";
    outfile << "// Generated pairing test values in standard (non-Montgomery) form\n";
    outfile << "// Using blst_fp12_inverse and blst_final_exp\n\n";
    
    // Write f coordinates (from blst input)
    outfile << "// Input f coordinates\n";
    for (int i = 0; i < 12; i++) {
        vec384 f_coord;
        vec_copy(f_coord, f_standard[i/6][(i%6)/2][i%2], sizeof(vec384));
        outfile << "const char* f_" << i << " = \"" << vec384_to_hex_string(f_coord) << "\";\n";
    }
    outfile << "\n";
    
    // Write f_inverse coordinates (from blst_fp12_inverse)
    outfile << "// Inverse f coordinates (from blst_fp12_inverse)\n";
    for (int i = 0; i < 12; i++) {
        vec384 f_inv_coord;
        vec_copy(f_inv_coord, f_inverse_standard[i/6][(i%6)/2][i%2], sizeof(vec384));
        outfile << "const char* f_inv_" << i << " = \"" << vec384_to_hex_string(f_inv_coord) << "\";\n";
    }
    outfile << "\n";
    
    // Write pairing result coordinates (from blst_final_exp)
    outfile << "// Pairing result coordinates (from blst_final_exp)\n";
    for (int i = 0; i < 12; i++) {
        vec384 result_coord;
        vec_copy(result_coord, result_standard[i/6][(i%6)/2][i%2], sizeof(vec384));
        outfile << "const char* pairing_result" << i << " = \"" << vec384_to_hex_string(result_coord) << "\";\n";
    }
    
    outfile.close();
    std::cout << "Values written to pairing_values.hpp using blst functions\n" << std::endl;
    
    return true;
}

int main() {
    bool all_passed = true;
    
    all_passed &= test_final_exp();
    
    if (all_passed) {
        std::cout << "=== ALL TESTS PASSED ===" << std::endl;
        return 0;
    } else {
        std::cout << "=== SOME TESTS FAILED ===" << std::endl;
        return 1;
    }
}

