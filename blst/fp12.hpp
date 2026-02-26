#pragma once
#include <array>
#include <string>
#include <cstring>
extern "C" {
#include "fields.h"
}

// FP class for flattened formulas
class FP {
public:
    vec384 d;

    FP() {
        vec_copy(d, ZERO_384, sizeof(vec384));
    }
    
    FP(const vec384 &v) {
        vec_copy(d, v, sizeof(vec384));
    }

    FP(const std::string &s) {
        // Parse hex string - must be exactly 96 hex digits (384 bits = 6 * 64 bits = 6 * 16 hex digits)
        // Hex string is big-endian (most significant digit first), but vec384 is little-endian (least significant limb first)
        if (s.size() < 2 || s.substr(0, 2) != "0x") {
            throw std::runtime_error("FP hex string must start with '0x'");
        }
        std::string hex_str = s.substr(2); // Remove "0x" prefix
        
        // Check if too long
        if (hex_str.size() > 96) {
            throw std::runtime_error("FP hex string too long: expected 96 hex digits, got " + std::to_string(hex_str.size()));
        }
        
        // Pad with leading zeros to make it exactly 96 hex digits
        while (hex_str.size() < 96) {
            hex_str = "0" + hex_str;
        }
        
        // Parse 6 limbs of 16 hex digits each
        // Hex string is big-endian, so first 16 digits are most significant -> go to d[5]
        // Last 16 digits are least significant -> go to d[0]
        for (int i = 0; i < 6; i++) {
            int hex_offset = (5 - i) * 16; // Reverse order: limb 0 gets last 16 digits, limb 5 gets first 16 digits
            std::string limb_str = hex_str.substr(hex_offset, 16);
            d[i] = std::stoull(limb_str, nullptr, 16);
        }
        mul_fp(d, d, BLS12_381_RR);
    }

    FP operator+(const FP &other) const {
        FP result;
        add_fp(result.d, d, other.d);
        return result;
    }

    FP operator-(const FP &other) const {
        FP result;
        sub_fp(result.d, d, other.d);
        return result;
    }

    FP operator*(const FP &other) const {
        FP result;
        mul_fp(result.d, d, other.d);
        return result;
    }

    FP operator*(int factor) const {
        if (factor < 0) {
            return -(*this * -factor);
        } else if (factor == 0) {
            return FP();
        } else if (factor == 1) {
            return *this;
        } else if (factor == 2) {
            return *this + *this;
        } else if (factor == 3) {
            return *this + *this + *this;
        } else {
            throw std::runtime_error("Invalid factor");
        }
    }
    
    FP operator-() const {
        FP result;
        neg_fp(result.d, d);
        return result;
    }
    
    bool operator==(const FP &other) const {
        for (int i = 0; i < 6; i++) {
            if (d[i] != other.d[i]) return false;
        }
        return true;
    }

    void to_vec384(vec384 out) const {
        vec_copy(out, d, sizeof(vec384));
    }
};

// RingFP wrapper
class RingFP {
public:
    using StandardElement = FP;

    auto negate(const FP& x) const {
        return -x;
    }

    template<int off>
    auto offsets() const {
        return FP(); // zero
    }

    FP from_hex(const std::string &s) const {
        return FP(s);
    }

    FP wrap(const FP &x) const {
        return x;
    }

    std::array<FP, 10> batch_reduce_expand10(const FP &poly0, const FP &poly1, const FP &poly2, const FP &poly3, const FP &poly4, const FP &poly5, const FP &poly6, const FP &poly7, const FP &poly8, const FP &poly9) const {
        return std::array<FP, 10> {poly0, poly1, poly2, poly3, poly4, poly5, poly6, poly7, poly8, poly9};
    }
    std::array<FP, 6> batch_reduce_expand6(const FP &poly0, const FP &poly1, const FP &poly2, const FP &poly3, const FP &poly4, const FP &poly5) const {
        return std::array<FP, 6> {poly0, poly1, poly2, poly3, poly4, poly5};
    }
    std::array<FP, 8> batch_reduce_expand8(const FP &poly0, const FP &poly1, const FP &poly2, const FP &poly3, const FP &poly4, const FP &poly5, const FP &poly6, const FP &poly7) const {
        return std::array<FP, 8> {poly0, poly1, poly2, poly3, poly4, poly5, poly6, poly7};
    }
    std::array<FP, 12> batch_reduce_expand12(const FP &poly0, const FP &poly1, const FP &poly2, const FP &poly3, const FP &poly4, const FP &poly5, const FP &poly6, const FP &poly7, const FP &poly8, const FP &poly9, const FP &poly10, const FP &poly11) const {
        return std::array<FP, 12> {poly0, poly1, poly2, poly3, poly4, poly5, poly6, poly7, poly8, poly9, poly10, poly11};
    }

    // batch_reduce_expand overloads for flattened_code.hpp (same semantics, variable arity)
    std::array<FP, 1> batch_reduce_expand(const FP &poly0) const {
        return std::array<FP, 1>{poly0};
    }
    std::array<FP, 2> batch_reduce_expand(const FP &poly0, const FP &poly1) const {
        return std::array<FP, 2>{poly0, poly1};
    }
    std::array<FP, 6> batch_reduce_expand(const FP &poly0, const FP &poly1, const FP &poly2, const FP &poly3, const FP &poly4, const FP &poly5) const {
        return std::array<FP, 6>{poly0, poly1, poly2, poly3, poly4, poly5};
    }
    std::array<FP, 8> batch_reduce_expand(const FP &poly0, const FP &poly1, const FP &poly2, const FP &poly3, const FP &poly4, const FP &poly5, const FP &poly6, const FP &poly7) const {
        return std::array<FP, 8>{poly0, poly1, poly2, poly3, poly4, poly5, poly6, poly7};
    }
    std::array<FP, 10> batch_reduce_expand(const FP &poly0, const FP &poly1, const FP &poly2, const FP &poly3, const FP &poly4, const FP &poly5, const FP &poly6, const FP &poly7, const FP &poly8, const FP &poly9) const {
        return std::array<FP, 10>{poly0, poly1, poly2, poly3, poly4, poly5, poly6, poly7, poly8, poly9};
    }
    std::array<FP, 12> batch_reduce_expand(const FP &poly0, const FP &poly1, const FP &poly2, const FP &poly3, const FP &poly4, const FP &poly5, const FP &poly6, const FP &poly7, const FP &poly8, const FP &poly9, const FP &poly10, const FP &poly11) const {
        return std::array<FP, 12>{poly0, poly1, poly2, poly3, poly4, poly5, poly6, poly7, poly8, poly9, poly10, poly11};
    }
};




