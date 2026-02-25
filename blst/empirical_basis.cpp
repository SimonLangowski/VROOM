#include <iostream>
#include <vector>
#include <array>
#include <functional>
#include <cstdint>
#include <iomanip>
#include <sstream>
#include <cassert>
extern "C" {
#include "fields.h"
#include "point.h"
}
#include "inversion_parts.hpp"

extern "C" {
void psi_wrapper_export(POINTonE2 *out, const POINTonE2 *in);
}

// Forward declarations
bool is_negative_384(const vec384 value);
bool is_zero_384(const vec384 value);

class Coefficient {
    public:
    bool is_small;
    int64_t value;  // For small values
    vec384 large_value;  // For 384-bit values
    int var0;
    int var1;
    
    Coefficient(int64_t v, int v0, int v1) : is_small(true), value(v), var0(v0), var1(v1) {
        // Initialize large_value to zero
        vec_copy(large_value, ZERO_384, sizeof(vec384));
    }
    
    Coefficient(const vec384 &v, int v0, int v1) : is_small(false), value(0), var0(v0), var1(v1) {
        vec_copy(large_value, v, sizeof(vec384));
    }
    
    void print() const {
        if (is_small) {
            std::cout << value;
        } else {
            // Print as hex for large values
            std::cout << "0x";
            bool leading = true;
            for (int i = 5; i >= 0; i--) {
                if (large_value[i] != 0 || !leading || i == 0) {
                    if (!leading) {
                        std::cout << std::setfill('0') << std::setw(16);
                    }
                    std::cout << std::hex << large_value[i];
                    leading = false;
                }
            }
            std::cout << std::dec;
        }
    }
    
    bool is_zero() const {
        if (is_small) {
            return value == 0;
        } else {
            // Check if all limbs are zero
            for (int i = 0; i < 6; i++) {
                if (large_value[i] != 0) return false;
            }
            return true;
        }
    }
};

class Polynomial {
    public:
    std::vector<Coefficient> coefficients;
    Polynomial() {
        coefficients = std::vector<Coefficient>();
    }

    void simplify() {
        for (size_t i = 0; i < coefficients.size(); i++) {
            for (size_t j = i + 1; j < coefficients.size(); j++) {
                if (coefficients[i].var0 == coefficients[j].var0 && coefficients[i].var1 == coefficients[j].var1) {
                    // Combine like terms
                    if (coefficients[i].is_small && coefficients[j].is_small) {
                        // Both small: add int64_t values
                        coefficients[i].value += coefficients[j].value;
                    } else if (!coefficients[i].is_small && !coefficients[j].is_small) {
                        // Both large: add vec384 values
                        vec384 sum;
                        add_mod_384(sum, coefficients[i].large_value, coefficients[j].large_value, BLS12_381_P);
                        vec_copy(coefficients[i].large_value, sum, sizeof(vec384));
                        coefficients[i].is_small = false;
                        coefficients[i].value = 0;
                    } else {
                        // Mixed: convert small to large and add
                        if (coefficients[i].is_small) {
                            // Convert i to large
                            vec384 i_large = {(limb_t)coefficients[i].value, 0, 0, 0, 0, 0};
                            if (coefficients[i].value < 0) {
                                // Handle negative: P + value
                                vec384 p_val;
                                sub_mod_384(p_val, BLS12_381_P, i_large, BLS12_381_P);
                                vec_copy(i_large, p_val, sizeof(vec384));
                            }
                            vec384 sum;
                            add_mod_384(sum, i_large, coefficients[j].large_value, BLS12_381_P);
                            vec_copy(coefficients[i].large_value, sum, sizeof(vec384));
                            coefficients[i].is_small = false;
                            coefficients[i].value = 0;
                        } else {
                            // Convert j to large
                            vec384 j_large = {(limb_t)coefficients[j].value, 0, 0, 0, 0, 0};
                            if (coefficients[j].value < 0) {
                                vec384 p_val;
                                sub_mod_384(p_val, BLS12_381_P, j_large, BLS12_381_P);
                                vec_copy(j_large, p_val, sizeof(vec384));
                            }
                            vec384 sum;
                            add_mod_384(sum, coefficients[i].large_value, j_large, BLS12_381_P);
                            vec_copy(coefficients[i].large_value, sum, sizeof(vec384));
                        }
                    }
                    coefficients.erase(coefficients.begin() + j);
                    j--;
                }
            }
            // Remove zero coefficients
            bool is_zero = false;
            if (coefficients[i].is_small) {
                is_zero = (coefficients[i].value == 0);
            } else {
                is_zero = is_zero_384(coefficients[i].large_value);
            }
            if (is_zero) {
                coefficients.erase(coefficients.begin() + i);
                i--;
            }
        }
    }

    void print(int var, bool all_x = false) {
        bool first = true;
        for (const auto &coefficient : coefficients) {
            if (coefficient.is_zero()) {
                continue;
            }
            if (!first) {
                std::cout << " ";
            }
            // For 384-bit coefficients, we don't handle negatives - values are in [0, P)
            bool is_positive = coefficient.is_small ? (coefficient.value > 0) : true;
            if (is_positive && !first) {
                std::cout << "+";
            }
            first = false;
            coefficient.print();
            if (coefficient.var0 != -1) {
                if (all_x) {
                    std::cout << " x" << coefficient.var0;
                } else {
                    std::cout << ((coefficient.var0 < var) ? " x" : " y") << coefficient.var0 % var;
                }
            }
            if (coefficient.var1 != -1) {
                if (all_x) {
                    std::cout << "x" << coefficient.var1;
                } else {
                    std::cout << ((coefficient.var1 < var) ? "x" : "y") << coefficient.var1 % var;
                }
            }
        }
        std::cout << std::endl;
    }

    size_t size() const {
        return coefficients.size();
    }
};

template<size_t length>
std::array<int64_t, length> array_subtract(const std::array<int64_t, length> &a, const std::array<int64_t, length> &b) {
    std::array<int64_t, length> result;
    for (size_t i = 0; i < length; i++) {
        result[i] = a[i] - b[i];
    }
    return result;
}

template<int basis_in, int basis_out>
void generate_basis(std::function<std::array<int64_t, basis_out>(const std::array<int64_t, basis_in> &)> basis_conversion, bool all_x = false) {
    std::array<Polynomial, basis_out> polynomials;
    // Get constant term
    std::array<int64_t, basis_in> input = {0};
    std::array<int64_t, basis_out> constant = basis_conversion(input);
    for (int k = 0; k < basis_out; k++) {
        polynomials[k].coefficients.push_back(Coefficient(constant[k], -1, -1));
    }

    
    std::array<std::array<int64_t, basis_out>, basis_in> offsets = {0};
    for (int i = 0; i < basis_in; i++) {
        std::array<int64_t, basis_in> input = {0};
        input[i] = 1;
        // gives s_i x_i^2 + l_i x_i  = s_i + l_i (- c to cancel out the constant term) terms
        offsets[i] = array_subtract(basis_conversion(input), constant);
        // Input of minus 1 gives s_i x_i^2 - l_i x_i = s_i - l_i
        input[i] = -1;
        std::array<int64_t, basis_out> minus_one = array_subtract(basis_conversion(input), constant);
        for (int k = 0; k < basis_out; k++) {
            int64_t sq_term = offsets[i][k] + minus_one[k];
            assert(sq_term % 2 == 0);
            sq_term /= 2; // for vec384 this will have to be modular division (multiplication by 2^-1 mod p)
            polynomials[k].coefficients.push_back(Coefficient(sq_term, i, i));
            int64_t linear_term = offsets[i][k] - minus_one[k];
            assert(linear_term % 2 == 0);
            linear_term /= 2;
            polynomials[k].coefficients.push_back(Coefficient(linear_term, i, -1));
        }
    }

    // x_i = 1, x_j = 1 gives x_i^2 + x_i + x_j^2 + x_j + x_i x_j + c terms
    // Need to subtract off x_i^2 + x_i and x_j^2 + x_j to get x_i x_j term
    for (int i = 0; i < basis_in; i++) {
        for (int j = i + 1; j < basis_in; j++) {
            std::array<int64_t, basis_in> input = {0};
            input[i] = 1;
            input[j] = 1;
            std::array<int64_t, basis_out> output = basis_conversion(input);
            for (int k = 0; k < basis_out; k++) {
                int64_t val = output[k] - offsets[i][k] - offsets[j][k] - constant[k];
                if (val != 0) {
                    polynomials[k].coefficients.push_back(Coefficient(val, i, j));
                }
            }
        }
    }
    for (int i = 0; i < basis_out; i++) {
        polynomials[i].simplify();
        std::cout << "Polynomial " << i << "(" << polynomials[i].size() << ") : ";
        polynomials[i].print(basis_out, all_x);
        std::cout << std::endl;
    }
}

// Compute 2^-1 mod P in regular form (not Montgomery)
// 2^-1 mod P = (P+1)/2 mod P
void compute_inv2(vec384 ret) {
    // Compute (P+1)/2 mod P
    // Since P is odd, (P+1) is even, so (P+1)/2 is an integer
    // P + 1 will not overflow the low limb, so we can add 1 directly
    vec_copy(ret, BLS12_381_P, sizeof(vec384));
    ret[0] += 1;
    
    // Divide by 2: right shift the entire number
    // Start from the most significant limb and propagate the carry
    // The LSB of each limb becomes the MSB carry for the next (lower) limb
    limb_t carry = 0;
    for (int i = 5; i >= 0; i--) {
        limb_t val = ret[i];
        // Right shift first: LSB becomes carry for next (lower) limb
        limb_t next_carry = val & 1; // Extract LSB as carry for next limb
        ret[i] = val >> 1; // Right shift by 1 (divide by 2)
        ret[i] += (carry << 63);
        carry = next_carry;
    }
    // Note: carry from limb 0 is discarded (it's 0 since P+1 is even)
}

// Version for 384-bit coefficients (for frobenius_map)
template<int basis_in, int basis_out>
void generate_basis_384(std::function<std::array<vec384, basis_out>(const std::array<int64_t, basis_in> &, size_t)> basis_conversion, size_t n) {
    std::array<Polynomial, basis_out> polynomials;
    
    // Precompute 2^-1 mod P in regular form (will convert to Montgomery when needed)
    vec384 inv2;
    compute_inv2(inv2);
    
    // Get constant term
    std::array<int64_t, basis_in> input = {0};
    std::array<vec384, basis_out> constant = basis_conversion(input, n);
    for (int k = 0; k < basis_out; k++) {
        vec384 const_val;
        from_fp(const_val, constant[k]);
        if (!is_zero_384(const_val)) {
            polynomials[k].coefficients.push_back(Coefficient(const_val, -1, -1));
        }
    }
    
    // Extract linear and square terms for each input variable
    std::array<std::array<vec384, basis_out>, basis_in> linear_terms;
    std::array<std::array<vec384, basis_out>, basis_in> square_terms;
    
    for (int i = 0; i < basis_in; i++) {
        std::array<int64_t, basis_in> input = {0};
        input[i] = 1;
        // gives x_i^2 + x_i + c
        std::array<vec384, basis_out> output_1 = basis_conversion(input, n);
        
        // Input of minus 1 gives x_i^2 - x_i + c
        input[i] = -1;
        std::array<vec384, basis_out> output_minus1 = basis_conversion(input, n);
        
        for (int k = 0; k < basis_out; k++) {
            vec384 val1, val_minus1, const_val, diff1, diff_minus1;
            from_fp(val1, output_1[k]);
            from_fp(val_minus1, output_minus1[k]);
            from_fp(const_val, constant[k]);
            sub_mod_384(diff1, val1, const_val, BLS12_381_P);
            sub_mod_384(diff_minus1, val_minus1, const_val, BLS12_381_P);
            
            // Square term: (f(1) + f(-1))/2 = (diff1 + diff_minus1) / 2
            vec384 sq_term, sum;
            add_mod_384(sum, diff1, diff_minus1, BLS12_381_P);
            // Convert to Montgomery form, multiply by 2^-1 mod P, then convert back
            vec384 sum_mont, inv2_mont, sq_term_mont;
            mul_fp(sum_mont, sum, BLS12_381_RR);
            mul_fp(inv2_mont, inv2, BLS12_381_RR);
            mul_fp(sq_term_mont, sum_mont, inv2_mont);
            from_fp(sq_term, sq_term_mont);
            vec_copy(square_terms[i][k], sq_term, sizeof(vec384));
            if (!is_zero_384(sq_term)) {
                polynomials[k].coefficients.push_back(Coefficient(sq_term, i, i));
            }
            
            // Linear term: (f(1) - f(-1))/2 = (diff1 - diff_minus1) / 2
            vec384 linear_term, diff_linear;
            sub_mod_384(diff_linear, diff1, diff_minus1, BLS12_381_P);
            // Convert to Montgomery form, multiply by 2^-1 mod P, then convert back
            vec384 diff_linear_mont, linear_term_mont;
            mul_fp(diff_linear_mont, diff_linear, BLS12_381_RR);
            mul_fp(linear_term_mont, diff_linear_mont, inv2_mont);
            from_fp(linear_term, linear_term_mont);
            vec_copy(linear_terms[i][k], linear_term, sizeof(vec384));
            if (!is_zero_384(linear_term)) {
                polynomials[k].coefficients.push_back(Coefficient(linear_term, i, -1));
            }
        }
    }
    
    // Extract cross terms: x_i * x_j
    // For linear functions, cross terms should only exist between components of the same Fp2 element
    for (int i = 0; i < basis_in; i++) {
        for (int j = i + 1; j < basis_in; j++) {
            std::array<int64_t, basis_in> input = {0};
            input[i] = 1;
            input[j] = 1;
            std::array<vec384, basis_out> output = basis_conversion(input, n);
            for (int k = 0; k < basis_out; k++) {
                vec384 val, const_val, linear_i, linear_j, sq_i, sq_j, cross_term;
                from_fp(val, output[k]);
                from_fp(const_val, constant[k]);
                
                // linear_terms and square_terms are already in regular form (not Montgomery)
                // so we can use them directly without from_fp
                vec_copy(linear_i, linear_terms[i][k], sizeof(vec384));
                vec_copy(linear_j, linear_terms[j][k], sizeof(vec384));
                vec_copy(sq_i, square_terms[i][k], sizeof(vec384));
                vec_copy(sq_j, square_terms[j][k], sizeof(vec384));
                
                // cross_term = f(1,1) - constant - linear_i - linear_j - sq_i - sq_j
                // This extracts the true cross term x_i * x_j
                // For a linear function, this should be zero
                sub_mod_384(cross_term, val, const_val, BLS12_381_P);
                sub_mod_384(cross_term, cross_term, linear_i, BLS12_381_P);
                sub_mod_384(cross_term, cross_term, linear_j, BLS12_381_P);
                sub_mod_384(cross_term, cross_term, sq_i, BLS12_381_P);
                sub_mod_384(cross_term, cross_term, sq_j, BLS12_381_P);
                
                // Only add cross term if it's non-zero
                // For linear functions, cross terms should be exactly zero
                if (!is_zero_384(cross_term)) {
                    polynomials[k].coefficients.push_back(Coefficient(cross_term, i, j));
                }
            }
        }
    }
    
    for (int i = 0; i < basis_out; i++) {
        polynomials[i].simplify();
        std::cout << "Polynomial " << i << "(" << polynomials[i].size() << ") : ";
        polynomials[i].print(basis_out, false);  // generate_basis_384 doesn't support all_x yet
        std::cout << std::endl;
    }
}

// Check if a 384-bit value is greater than P/2 (interpret as negative)
// Compare limb by limb from high to low
bool is_negative_384(const vec384 value) {
    vec384 ret;
    from_fp(ret, value);
    
    // Compare with P/2 by comparing 2*ret with P
    // If 2*ret > P, then ret > P/2
    // Compare ret with P/2 directly by comparing limbs
    // P/2 ≈ 0x0d0088f51cbfff34d... (half of each limb, with carries)
    
    // Compare from most significant to least significant limb
    for (int i = 5; i >= 0; i--) {
        limb_t p_half_i = BLS12_381_P[i] / 2;
        if (i < 5 && (BLS12_381_P[i] % 2 == 1)) {
            // If lower limb was odd, add carry to this limb's half
            p_half_i += (1ULL << 63);
        }
        if (ret[i] > p_half_i) {
            return true;
        }
        if (ret[i] < p_half_i) {
            return false;
        }
    }
    return false; // Equal to P/2, not negative
}

bool is_zero_384(const vec384 value) {
    // Check if value is zero (works for both Montgomery and regular form)
    for (int i = 0; i < 6; i++) {
        if (value[i] != 0) return false;
    }
    return true;
}

int64_t from_mont(vec384 value) {
    vec384 ret;
    from_fp(ret, value);
    int64_t result = (int64_t)ret[0];
    // If value > P/2, interpret as negative
    if (is_negative_384(value)) {
        result -= (int64_t)BLS12_381_P[0];
    }
    return result;
}


void to_mont(vec384 ret, int64_t value) {
    if (value == 0) {
        // Zero in Montgomery form
        vec_copy(ret, ZERO_384, sizeof(vec384));
    } else {
        // Convert value to Montgomery form: value * R mod P
        // where R = 2^384 mod P (BLS12_381_Rx.p)
        if (value > 0) {
            vec384 value_vec = {(limb_t)value, 0, 0, 0, 0, 0};
            mul_fp(ret, value_vec, BLS12_381_RR);
        } else {
            vec384 value_vec = {(limb_t)(-value), 0, 0, 0, 0, 0};
            mul_fp(ret, value_vec, BLS12_381_RR);
            neg_fp(ret, ret);
        }
    }
}

std::array<int64_t, 2> mul_fp2_wrapper(const std::array<int64_t, 4> &ab) {
    vec384x a, b, c;
    to_mont(a[0], ab[0]);
    to_mont(a[1], ab[1]);
    to_mont(b[0], ab[2]);
    to_mont(b[1], ab[3]);
    mul_fp2(c, a, b);
    return {from_mont(c[0]), from_mont(c[1])};
}

std::array<int64_t, 2> sqr_fp2_wrapper(const std::array<int64_t, 2> &ab) {
    vec384x a, c;
    to_mont(a[0], ab[0]);
    to_mont(a[1], ab[1]);
    sqr_fp2(c, a);
    return {from_mont(c[0]), from_mont(c[1])};
}

std::array<int64_t, 12> mul_fp12_wrapper(const std::array<int64_t, 24> &ab) {
    vec384fp12 a, b, c;
    // vec384fp12 is vec384fp6[2], where vec384fp6 is vec384x[3], and vec384x is vec384[2]
    // So to access element i (0-11): a[i/6][(i%6)/2][i%2]
    for (int i = 0; i < 12; i++) {
        to_mont(a[i/6][(i%6)/2][i%2], ab[i]);
        to_mont(b[i/6][(i%6)/2][i%2], ab[i + 12]);
    }
    mul_fp12_export(c, a, b);
    std::array<int64_t, 12> result;
    for (int i = 0; i < 12; i++) {
        result[i] = from_mont(c[i/6][(i%6)/2][i%2]);
    }
    return result;
}

std::array<int64_t, 12> sqr_fp12_wrapper(const std::array<int64_t, 12> &ab) {
    vec384fp12 a, c;
    for (int i = 0; i < 12; i++) {
        to_mont(a[i/6][(i%6)/2][i%2], ab[i]);
    }
    sqr_fp12_export(c, a);
    std::array<int64_t, 12> result;
    for (int i = 0; i < 12; i++) {
        result[i] = from_mont(c[i/6][(i%6)/2][i%2]);
    }
    return result;
}

std::array<int64_t, 12> cyclotomic_sqr_fp12_wrapper(const std::array<int64_t, 12> &ab) {
    vec384fp12 a, c;
    for (int i = 0; i < 12; i++) {
        to_mont(a[i/6][(i%6)/2][i%2], ab[i]);
    }
    cyclotomic_sqr_fp12_export(c, a);
    std::array<int64_t, 12> result;
    for (int i = 0; i < 12; i++) {
        result[i] = from_mont(c[i/6][(i%6)/2][i%2]);
    }
    return result;
}

std::array<int64_t, 12> conjugate_fp12_wrapper(const std::array<int64_t, 12> &ab) {
    vec384fp12 a;
    for (int i = 0; i < 12; i++) {
        to_mont(a[i/6][(i%6)/2][i%2], ab[i]);
    }
    conjugate_fp12_export(a);
    std::array<int64_t, 12> result;
    for (int i = 0; i < 12; i++) {
        result[i] = from_mont(a[i/6][(i%6)/2][i%2]);
    }
    return result;
}

// For frobenius_map, coefficients are 384-bit, so we return vec384 (in Montgomery form, not converted)
// The generate_basis_384 function will convert from Montgomery when needed
std::array<vec384, 12> frobenius_map_fp12_wrapper(const std::array<int64_t, 12> &ab, size_t n) {
    vec384fp12 a, c;
    for (int i = 0; i < 12; i++) {
        to_mont(a[i/6][(i%6)/2][i%2], ab[i]);
    }
    frobenius_map_fp12_export(c, a, n);
    std::array<vec384, 12> result;
    for (int i = 0; i < 12; i++) {
        // Keep in Montgomery form for generate_basis_384
        vec_copy(result[i], c[i/6][(i%6)/2][i%2], sizeof(vec384));
    }
    return result;
}

// mul_fp6: multiplies two fp6 elements
// Input: 12 elements (6 for a, 6 for b)
// vec384fp6 is vec384x[3], where vec384x is vec384[2]
// So to access element i (0-5): a[i/2][i%2]
std::array<int64_t, 6> mul_fp6_wrapper(const std::array<int64_t, 12> &ab) {
    vec384fp6 a, b, c;
    for (int i = 0; i < 6; i++) {
        to_mont(a[i/2][i%2], ab[i]);
        to_mont(b[i/2][i%2], ab[i + 6]);
    }
    mul_fp6_export(c, a, b);
    std::array<int64_t, 6> result;
    for (int i = 0; i < 6; i++) {
        result[i] = from_mont(c[i/2][i%2]);
    }
    return result;
}

// sqr_fp4: squares an fp4 element
// Input: 4 elements (2 for a0, 2 for a1)
// vec384fp4 is vec384x[2], where vec384x is vec384[2]
// Input: [a0[0], a0[1], a1[0], a1[1]]
// Output: vec384fp4 = vec384x[2] = vec384[2][2]
std::array<int64_t, 4> sqr_fp4_wrapper(const std::array<int64_t, 4> &ab) {
    vec384x a0, a1;
    vec384fp4 c;
    to_mont(a0[0], ab[0]);
    to_mont(a0[1], ab[1]);
    to_mont(a1[0], ab[2]);
    to_mont(a1[1], ab[3]);
    sqr_fp4_export(c, a0, a1);
    std::array<int64_t, 4> result;
    // c is vec384x[2], so c[0] and c[1] are vec384x (each is vec384[2])
    result[0] = from_mont(c[0][0]);
    result[1] = from_mont(c[0][1]);
    result[2] = from_mont(c[1][0]);
    result[3] = from_mont(c[1][1]);
    return result;
}

// mul_by_xy00z0_fp12: multiplies fp12 by fp6
// Input: 12 elements for fp12, 6 elements for fp6 = 18 total
// vec384fp6 is vec384x[3], where vec384x is vec384[2]
// So to access element i (0-5) in fp6: xy00z0[i/2][i%2]
std::array<int64_t, 12> mul_by_xy00z0_fp12_wrapper(const std::array<int64_t, 18> &ab) {
    vec384fp12 a, c;
    vec384fp6 xy00z0;
    // First 12 elements are the fp12
    for (int i = 0; i < 12; i++) {
        to_mont(a[i/6][(i%6)/2][i%2], ab[i]);
    }
    // Next 6 elements are the fp6
    for (int i = 0; i < 6; i++) {
        to_mont(xy00z0[i/2][i%2], ab[i + 12]);
    }
    mul_by_xy00z0_fp12_export(c, a, xy00z0);
    std::array<int64_t, 12> result;
    for (int i = 0; i < 12; i++) {
        result[i] = from_mont(c[i/6][(i%6)/2][i%2]);
    }
    return result;
}

// inverse_fp12_to_fp6: takes fp12 (12 elements) -> outputs fp6 (6 elements)
std::array<int64_t, 6> inverse_fp12_to_fp6_wrapper(const std::array<int64_t, 12> &a) {
    vec384fp12 a_fp12;
    vec384fp6 ret;
    // Convert input: fp12 is vec384fp6[2], where vec384fp6 is vec384x[3], and vec384x is vec384[2]
    // So to access element i (0-11): a_fp12[i/6][(i%6)/2][i%2]
    for (int i = 0; i < 12; i++) {
        to_mont(a_fp12[i/6][(i%6)/2][i%2], a[i]);
    }
    inverse_fp12_to_fp6_wrapper_export(ret, a_fp12);
    std::array<int64_t, 6> result;
    // ret is vec384fp6 = vec384x[3] = vec384[2][3]
    // So to access element i (0-5): ret[i/2][i%2]
    for (int i = 0; i < 6; i++) {
        result[i] = from_mont(ret[i/2][i%2]);
    }
    return result;
}

// inverse_fp6_to_fp12: takes fp12 (12 elements) and fp6 (6 elements) = 18 inputs -> outputs fp12 (12 elements)
// Input order: [a (fp12), t1 (fp6)] so that a is x (12 elements) and t1 is y (6 elements)
std::array<int64_t, 12> inverse_fp6_to_fp12_wrapper(const std::array<int64_t, 18> &ab) {
    vec384fp6 t1;
    vec384fp12 a, ret;
    // First 12 elements are a (fp12) - these become x variables
    for (int i = 0; i < 12; i++) {
        to_mont(a[i/6][(i%6)/2][i%2], ab[i]);
    }
    // Next 6 elements are t1 (fp6) - these become y variables
    for (int i = 0; i < 6; i++) {
        to_mont(t1[i/2][i%2], ab[i + 12]);
    }
    inverse_fp6_to_fp12_wrapper_export(ret, t1, a);
    std::array<int64_t, 12> result;
    for (int i = 0; i < 12; i++) {
        result[i] = from_mont(ret[i/6][(i%6)/2][i%2]);
    }
    return result;
}

// inverse_fp6_to_fp2 part 1: fp6 (6 elements) -> c0, c1, c2 (6 elements). Quadratic in a.
std::array<int64_t, 6> inverse_fp6_to_fp2_c0_c1_c2_wrapper(const std::array<int64_t, 6> &a) {
    vec384fp6 a_fp6;
    vec384fp2 c0_out, c1_out, c2_out;
    for (int i = 0; i < 6; i++) {
        to_mont(a_fp6[i/2][i%2], a[i]);
    }
    inverse_fp6_to_fp2_c0_c1_c2_wrapper_export(c0_out, c1_out, c2_out, a_fp6);
    std::array<int64_t, 6> result;
    result[0] = from_mont(c0_out[0]);
    result[1] = from_mont(c0_out[1]);
    result[2] = from_mont(c1_out[0]);
    result[3] = from_mont(c1_out[1]);
    result[4] = from_mont(c2_out[0]);
    result[5] = from_mont(c2_out[1]);
    return result;
}

// inverse_fp6_to_fp2 part 2: c0, c1, c2 (6 elements) and a (6 elements) = 12 inputs -> ret (2). Bilinear in (c,a), quadratic for empirical basis.
std::array<int64_t, 2> inverse_fp6_to_fp2_ret_wrapper(const std::array<int64_t, 12> &ca) {
    vec384fp2 c0, c1, c2, ret;
    vec384fp6 a_fp6;
    to_mont(c0[0], ca[0]);
    to_mont(c0[1], ca[1]);
    to_mont(c1[0], ca[2]);
    to_mont(c1[1], ca[3]);
    to_mont(c2[0], ca[4]);
    to_mont(c2[1], ca[5]);
    for (int i = 0; i < 6; i++) {
        to_mont(a_fp6[i/2][i%2], ca[6 + i]);
    }
    inverse_fp6_to_fp2_ret_wrapper_export(ret, c0, c1, c2, a_fp6);
    std::array<int64_t, 2> result;
    result[0] = from_mont(ret[0]);
    result[1] = from_mont(ret[1]);
    return result;
}

// inverse_fp2_to_fp6: takes fp2 c0, c1, c2 (6 elements) and fp2 t1 (2 elements) = 8 inputs -> outputs fp6 (6 elements)
// Input order: [c0, c1, c2 (6 elements), t1 (2 elements)] so that c0, c1, c2 are x (6 elements) and t1 is y (2 elements)
std::array<int64_t, 6> inverse_fp2_to_fp6_wrapper(const std::array<int64_t, 8> &inputs) {
    vec384fp2 t1, c0, c1, c2;
    vec384fp6 ret;
    // First 2 elements are c0 - these become x variables
    to_mont(c0[0], inputs[0]);
    to_mont(c0[1], inputs[1]);
    // Next 2 elements are c1 - these become x variables
    to_mont(c1[0], inputs[2]);
    to_mont(c1[1], inputs[3]);
    // Next 2 elements are c2 - these become x variables
    to_mont(c2[0], inputs[4]);
    to_mont(c2[1], inputs[5]);
    // Last 2 elements are t1 - these become y variables
    to_mont(t1[0], inputs[6]);
    to_mont(t1[1], inputs[7]);
    inverse_fp2_to_fp6_wrapper_export(ret, t1, c0, c1, c2);
    std::array<int64_t, 6> result;
    for (int i = 0; i < 6; i++) {
        result[i] = from_mont(ret[i/2][i%2]);
    }
    return result;
}

// inverse_fp2_to_fp: takes fp2 (2 elements) -> outputs norm (1 element)
// Computes a^2 + b^2
std::array<int64_t, 1> inverse_fp2_to_fp_wrapper(const std::array<int64_t, 2> &a) {
    vec384fp2 a_fp2;
    vec384 norm_out;
    to_mont(a_fp2[0], a[0]);
    to_mont(a_fp2[1], a[1]);
    inverse_fp2_to_fp_wrapper_export(norm_out, a_fp2);
    std::array<int64_t, 1> result;
    result[0] = from_mont(norm_out);
    return result;
}

// inverse_fp_to_fp2: takes fp2 (a, b) (2 elements) and fp r (1 element) = 3 inputs -> outputs fp2 (2 elements)
// Input order: [a, b, r] so that r is treated as x and (a, b) as y in empirical basis
// Computes a*r and -b*r
std::array<int64_t, 2> inverse_fp_to_fp2_wrapper(const std::array<int64_t, 3> &inputs) {
    vec384 r;
    vec384fp2 a_fp2, ret;
    // inputs[0] = a, inputs[1] = b (y variables)
    // inputs[2] = r (x variable)
    to_mont(a_fp2[0], inputs[0]);
    to_mont(a_fp2[1], inputs[1]);
    to_mont(r, inputs[2]);
    inverse_fp_to_fp2_wrapper_export(ret, r, a_fp2);
    std::array<int64_t, 2> result;
    result[0] = from_mont(ret[0]);
    result[1] = from_mont(ret[1]);
    return result;
}

// For psi_neg, coefficients are 384-bit, so we return vec384 (in Montgomery form, not converted)
// The generate_basis_384 function will convert from Montgomery when needed
std::array<vec384, 6> psi_neg_wrapper(const std::array<int64_t, 6> &a, size_t n) {
    (void)n; // Unused parameter, kept for compatibility with generate_basis_384
    POINTonE2 a_e2;
    POINTonE2 ret;
    // Input: a[0]=X[0], a[1]=X[1], a[2]=Y[0], a[3]=Y[1], a[4]=Z[0], a[5]=Z[1]
    to_mont(a_e2.X[0], a[0]);
    to_mont(a_e2.X[1], a[1]);
    to_mont(a_e2.Y[0], a[2]);
    to_mont(a_e2.Y[1], a[3]);
    to_mont(a_e2.Z[0], a[4]);
    to_mont(a_e2.Z[1], a[5]);
    psi_wrapper_export(&ret, &a_e2);
    std::array<vec384, 6> result;
    // Keep in Montgomery form for generate_basis_384
    // Output: result[0]=X[0], result[1]=X[1], result[2]=Y[0], result[3]=Y[1], result[4]=Z[0], result[5]=Z[1]
    vec_copy(result[0], ret.X[0], sizeof(vec384));
    vec_copy(result[1], ret.X[1], sizeof(vec384));
    vec_copy(result[2], ret.Y[0], sizeof(vec384));
    vec_copy(result[3], ret.Y[1], sizeof(vec384));
    vec_copy(result[4], ret.Z[0], sizeof(vec384));
    vec_copy(result[5], ret.Z[1], sizeof(vec384));
    return result;
}

int main() {
    std::cout << "Generating basis for Fp2 multiplication..." << std::endl;
    generate_basis<4, 2>(mul_fp2_wrapper);
    std::cout << "Generating basis for Fp2 squaring..." << std::endl;
    generate_basis<2, 2>(sqr_fp2_wrapper);
    std::cout << "Generating basis for Fp12 multiplication..." << std::endl;
    generate_basis<24, 12>(mul_fp12_wrapper);
    std::cout << "Generating basis for Fp12 squaring..." << std::endl;
    generate_basis<12, 12>(sqr_fp12_wrapper);
    std::cout << "Generating basis for Fp12 cyclotomic squaring..." << std::endl;
    generate_basis<12, 12>(cyclotomic_sqr_fp12_wrapper);
    std::cout << "Generating basis for Fp12 conjugation..." << std::endl;
    generate_basis<12, 12>(conjugate_fp12_wrapper);
    std::cout << "Generating basis for Fp12 Frobenius map (n=1)..." << std::endl;
    generate_basis_384<12, 12>(frobenius_map_fp12_wrapper, 1);
    std::cout << "Generating basis for Fp12 Frobenius map (n=2)..." << std::endl;
    generate_basis_384<12, 12>(frobenius_map_fp12_wrapper, 2);
    std::cout << "Generating basis for Fp12 Frobenius map (n=3)..." << std::endl;
    generate_basis_384<12, 12>(frobenius_map_fp12_wrapper, 3);
    std::cout << "Generating basis for Fp12 mul_by_xy00z0..." << std::endl;
    generate_basis<18, 12>(mul_by_xy00z0_fp12_wrapper);
    std::cout << "Generating basis for Fp6 multiplication..." << std::endl;
    generate_basis<12, 6>(mul_fp6_wrapper);
    std::cout << "Generating basis for Fp4 squaring..." << std::endl;
    generate_basis<4, 4>(sqr_fp4_wrapper);
    std::cout << "Generating basis for inverse_fp12_to_fp6..." << std::endl;
    generate_basis<12, 6>(inverse_fp12_to_fp6_wrapper, true);  // All x (fp12 only)
    std::cout << "Generating basis for inverse_fp6_to_fp12..." << std::endl;
    generate_basis<18, 12>(inverse_fp6_to_fp12_wrapper, false);  // x (fp12) and y (fp6)
    std::cout << "Generating basis for inverse_fp6_to_fp2 (part 1: a -> c0,c1,c2)..." << std::endl;
    generate_basis<6, 6>(inverse_fp6_to_fp2_c0_c1_c2_wrapper, true);  // Quadratic in a
    std::cout << "Generating basis for inverse_fp6_to_fp2 (part 2: c0,c1,c2,a -> ret)..." << std::endl;
    generate_basis<12, 2>(inverse_fp6_to_fp2_ret_wrapper, true);  // Bilinear in (c,a), quadratic
    std::cout << "Generating basis for inverse_fp2_to_fp6..." << std::endl;
    generate_basis<8, 6>(inverse_fp2_to_fp6_wrapper, false);  // x (c0,c1,c2) and y (t1)
    std::cout << "Generating basis for inverse_fp2_to_fp..." << std::endl;
    generate_basis<2, 1>(inverse_fp2_to_fp_wrapper, true);  // All x (fp2 only)
    std::cout << "Generating basis for inverse_fp_to_fp2..." << std::endl;
    generate_basis<3, 2>(inverse_fp_to_fp2_wrapper, false);  // x (fp2) and y (reciprocal)
    std::cout << "Generating basis for psi_neg..." << std::endl;
    generate_basis_384<6, 6>(psi_neg_wrapper, 0);  // POINTonE2: X, Y, Z (6 components) -> X, Y, Z (6 components), using vec384
    return 0;
}