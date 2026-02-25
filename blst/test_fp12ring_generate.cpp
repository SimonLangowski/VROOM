#include <iostream>
#include <array>
#include <cstdint>
#include <cstring>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <cctype>
#include <random>
extern "C" {
#include "fields.h"
void blst_fp12_mul(vec384fp12 ret, const vec384fp12 a, const vec384fp12 b);
void blst_fp12_cyclotomic_sqr(vec384fp12 ret, const vec384fp12 a);
void blst_fp12_frobenius_map(vec384fp12 ret, const vec384fp12 a, size_t n);
void blst_fp12_conjugate(vec384fp12 a);
}
#include "fp12.hpp"
#include "pairing_values.hpp"

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

void generate_test_values() {
    std::cout << "=== Generating FP12Ring test values using blst ===" << std::endl;
    
    // Use seeded values from pairing_values.hpp
    // Convert f values from hex strings to Montgomery form
    std::array<vec384, 12> f_mont;
    const char* f_hex[12] = {f_0, f_1, f_2, f_3, f_4, f_5, f_6, f_7, f_8, f_9, f_10, f_11};
    for (int i = 0; i < 12; i++) {
        std::string hex_str = std::string(f_hex[i]);
        if (hex_str.size() < 2 || (hex_str.substr(0, 2) != "0x" && hex_str.substr(0, 2) != "0X")) {
            hex_str = "0x" + hex_str;
        }
        FP fp(hex_str);
        fp.to_vec384(f_mont[i]);
    }
    
    // Convert to vec384fp12 format
    vec384fp12 f_standard;
    for (int i = 0; i < 12; i++) {
        vec_copy(f_standard[i/6][(i%6)/2][i%2], f_mont[i], sizeof(vec384));
    }
    
    // Convert f_inverse
    std::array<vec384, 12> f_inv_mont;
    const char* f_inv_hex[12] = {f_inv_0, f_inv_1, f_inv_2, f_inv_3, f_inv_4, f_inv_5, 
                                 f_inv_6, f_inv_7, f_inv_8, f_inv_9, f_inv_10, f_inv_11};
    for (int i = 0; i < 12; i++) {
        std::string hex_str = std::string(f_inv_hex[i]);
        if (hex_str.size() < 2 || (hex_str.substr(0, 2) != "0x" && hex_str.substr(0, 2) != "0X")) {
            hex_str = "0x" + hex_str;
        }
        FP fp(hex_str);
        fp.to_vec384(f_inv_mont[i]);
    }
    
    vec384fp12 f_inv_standard;
    for (int i = 0; i < 12; i++) {
        vec_copy(f_inv_standard[i/6][(i%6)/2][i%2], f_inv_mont[i], sizeof(vec384));
    }
    
    // Write test values to header file
    std::ofstream outfile("fp12ring_test_values.hpp");
    outfile << "#pragma once\n";
    outfile << "// Generated FP12Ring test values using blst (standard form)\n\n";
    
    // Test 1: fp12_flat_mul (f * f_inv = 1)
    std::cout << "Generating fp12_flat_mul test (f * f_inv)..." << std::endl;
    vec384fp12 mul_result_standard;
    blst_fp12_mul(mul_result_standard, f_standard, f_inv_standard);
    outfile << "// fp12_flat_mul(f, f_inv) result\n";
    for (int i = 0; i < 12; i++) {
        vec384 coord;
        vec_copy(coord, mul_result_standard[i/6][(i%6)/2][i%2], sizeof(vec384));
        outfile << "const char* fp12_mul_result_" << i << " = \"" << vec384_to_hex_string(coord) << "\";\n";
    }
    outfile << "\n";
    
    // Test 1b: fp12_flat_mul with different inputs (non-trivial result)
    // Use f * f (squaring) instead of f * f_inv to get a non-trivial result
    std::cout << "Generating fp12_flat_mul test (f * f)..." << std::endl;
    vec384fp12 mul_ff_result_standard;
    blst_fp12_mul(mul_ff_result_standard, f_standard, f_standard);
    outfile << "// fp12_flat_mul(f, f) result (non-trivial multiplication)\n";
    for (int i = 0; i < 12; i++) {
        vec384 coord;
        vec_copy(coord, mul_ff_result_standard[i/6][(i%6)/2][i%2], sizeof(vec384));
        outfile << "const char* fp12_mul_ff_result_" << i << " = \"" << vec384_to_hex_string(coord) << "\";\n";
    }
    outfile << "\n";
    
    // Test 2: fp12_flat_cyclotomic_sqr
    std::cout << "Generating fp12_flat_cyclotomic_sqr test..." << std::endl;
    vec384fp12 sqr_result_standard;
    blst_fp12_cyclotomic_sqr(sqr_result_standard, f_standard);
    outfile << "// fp12_flat_cyclotomic_sqr(f) result\n";
    for (int i = 0; i < 12; i++) {
        vec384 coord;
        vec_copy(coord, sqr_result_standard[i/6][(i%6)/2][i%2], sizeof(vec384));
        outfile << "const char* fp12_cyclotomic_sqr_result_" << i << " = \"" << vec384_to_hex_string(coord) << "\";\n";
    }
    outfile << "\n";
    
    // Test 3: fp12_flat_conjugate
    std::cout << "Generating fp12_flat_conjugate test..." << std::endl;
    vec384fp12 conj_result_standard;
    vec_copy(conj_result_standard, f_standard, sizeof(vec384fp12));
    blst_fp12_conjugate(conj_result_standard);
    outfile << "// fp12_flat_conjugate(f) result\n";
    for (int i = 0; i < 12; i++) {
        vec384 coord;
        vec_copy(coord, conj_result_standard[i/6][(i%6)/2][i%2], sizeof(vec384));
        outfile << "const char* fp12_conjugate_result_" << i << " = \"" << vec384_to_hex_string(coord) << "\";\n";
    }
    outfile << "\n";
    
    // Test 4: fp12_flat_frobenius_map1
    std::cout << "Generating fp12_flat_frobenius_map1 test..." << std::endl;
    vec384fp12 frob1_result_standard;
    blst_fp12_frobenius_map(frob1_result_standard, f_standard, 1);
    outfile << "// fp12_flat_frobenius_map1(f) result\n";
    for (int i = 0; i < 12; i++) {
        vec384 coord;
        vec_copy(coord, frob1_result_standard[i/6][(i%6)/2][i%2], sizeof(vec384));
        outfile << "const char* fp12_frobenius_map1_result_" << i << " = \"" << vec384_to_hex_string(coord) << "\";\n";
    }
    outfile << "\n";
    
    // Test 5: fp12_flat_frobenius_map2
    std::cout << "Generating fp12_flat_frobenius_map2 test..." << std::endl;
    vec384fp12 frob2_result_standard;
    blst_fp12_frobenius_map(frob2_result_standard, f_standard, 2);
    outfile << "// fp12_flat_frobenius_map2(f) result\n";
    for (int i = 0; i < 12; i++) {
        vec384 coord;
        vec_copy(coord, frob2_result_standard[i/6][(i%6)/2][i%2], sizeof(vec384));
        outfile << "const char* fp12_frobenius_map2_result_" << i << " = \"" << vec384_to_hex_string(coord) << "\";\n";
    }
    outfile << "\n";
    
    // Test 6: fp12_flat_frobenius_map3
    std::cout << "Generating fp12_flat_frobenius_map3 test..." << std::endl;
    vec384fp12 frob3_result_standard;
    blst_fp12_frobenius_map(frob3_result_standard, f_standard, 3);
    outfile << "// fp12_flat_frobenius_map3(f) result\n";
    for (int i = 0; i < 12; i++) {
        vec384 coord;
        vec_copy(coord, frob3_result_standard[i/6][(i%6)/2][i%2], sizeof(vec384));
        outfile << "const char* fp12_frobenius_map3_result_" << i << " = \"" << vec384_to_hex_string(coord) << "\";\n";
    }
    outfile << "\n";
    
    outfile.close();
    std::cout << "Test values written to fp12ring_test_values.hpp\n" << std::endl;
}

int main() {
    generate_test_values();
    return 0;
}

