/*
 * C++ wrapper for vec384 type
 * Provides std::array-based Vec384 class compatible with C array types
 */

#ifndef __VEC384_HPP__
#define __VEC384_HPP__

#include <array>
#include <cstdint>
#include <cstring>

#ifdef __cplusplus

// C++ wrapper class for vec384 (6 x 64-bit limbs)
// Note: Include fields.h before using the arithmetic member functions
class Vec384 {
public:
    using ArrayType = std::array<uint64_t, 6>;
    
    // Default constructor - zero-initialized
    Vec384() : data_{} {}
    
    // Constructor from C array
    Vec384(const uint64_t arr[6]) {
        std::memcpy(data_.data(), arr, sizeof(data_));
    }
    
    // Constructor from std::array
    Vec384(const ArrayType& arr) : data_(arr) {}
    
    // Implicit conversion to C array pointer (for C function calls)
    operator uint64_t*() { return data_.data(); }
    operator const uint64_t*() const { return data_.data(); }
    
    // Access operators
    uint64_t& operator[](size_t i) { return data_[i]; }
    const uint64_t& operator[](size_t i) const { return data_[i]; }
    
    // Get underlying array
    ArrayType& array() { return data_; }
    const ArrayType& array() const { return data_; }
    
    // Pointer access
    uint64_t* data() { return data_.data(); }
    const uint64_t* data() const { return data_.data(); }
    
    // Size
    static constexpr size_t size() { return 6; }
    
    // Field arithmetic operations (require fields.h to be included)
    // These functions operate on BLS12-381 field Fp
    // All operations are non-modifying and return new Vec384 instances
    // Note: Include fields.h before using these methods
    
    // Addition: result = this + other (non-modifying)
    Vec384 operator+(const Vec384& other) const;
    
    // Subtraction: result = this - other (non-modifying)
    Vec384 operator-(const Vec384& other) const;
    
    // Multiplication: result = this * other (non-modifying)
    Vec384 operator*(const Vec384& other) const;
    
    // Square: result = this^2 (non-modifying)
    Vec384 sqr() const;
    
    // Multiply by 3: result = this * 3 (non-modifying)
    Vec384 mul_by_3() const;
    
    // Multiply by 8: result = this * 8 (non-modifying)
    Vec384 mul_by_8() const;
    
    // Left shift: result = this << count (non-modifying)
    Vec384 lshift(size_t count) const;
    
    // Right shift: result = this >> count (non-modifying)
    Vec384 rshift(size_t count) const;
    
    // Divide by 2: result = this / 2 (non-modifying)
    Vec384 div_by_2() const;
    
    // Conditional negate: result = -this if flag is true (non-modifying)
    Vec384 cneg(uint64_t flag) const;
    
    // Negate: result = -this (non-modifying)
    Vec384 neg() const;
    
    // From Montgomery: result = convert from Montgomery representation (non-modifying)
    Vec384 from() const;
    
    // Reduce: static method to reduce vec768 (12 limbs) to vec384
    static Vec384 redc(const uint64_t a[12]);

private:
    ArrayType data_;
};

// C++ wrapper class for vec384x (Fp2: 2 x 64-bit limbs)
// Provides high-level Fp2 operations that mirror the shortcuts in fields.h
class Vec384x {
public:
    // Underlying storage: 2 * 64-bit limbs (matches C type vec384x)
    using ArrayType = std::array<uint64_t, 12>;

    // Default constructor - zero-initialized
    Vec384x() : data_{} {}

    // Constructor from C-style 2x6 array
    Vec384x(const uint64_t arr[2][6]) {
        for (size_t i = 0; i < 2; i++)
            std::memcpy(&data_[i * 6], arr[i], 6 * sizeof(uint64_t));
    }

    // Constructor from std::array
    Vec384x(const ArrayType& arr) : data_(arr) {}

    // Access individual components as Vec384 views
    Vec384 operator[](size_t i) const {
        Vec384 v;
        std::memcpy(v.data(), &data_[i * 6], 6 * sizeof(uint64_t));
        return v;
    }

    // Direct access to underlying limbs
    uint64_t* data() { return data_.data(); }
    const uint64_t* data() const { return data_.data(); }

    ArrayType& array() { return data_; }
    const ArrayType& array() const { return data_; }

    // Field arithmetic operations (Fp2), mirroring fields.h shortcuts.
    // All operations are non-modifying and return new Vec384x instances.

    // Addition: result = this + other
    Vec384x operator+(const Vec384x& other) const;

    // Subtraction: result = this - other
    Vec384x operator-(const Vec384x& other) const;

    // Multiply by 3: result = this * 3
    Vec384x mul_by_3() const;

    // Multiply by 8: result = this * 8
    Vec384x mul_by_8() const;

    // Left shift: result = this << count
    Vec384x lshift(size_t count) const;

    // Multiplication: result = this * other
    Vec384x mul(const Vec384x& other) const;

    // Square: result = this^2
    Vec384x sqr() const;

    // Conditional negate: result = -this if flag is true
    Vec384x cneg(uint64_t flag) const;

    // Negate: result = -this
    Vec384x neg() const { return cneg(1); }

private:
    ArrayType data_;
};

// Type aliases for higher towers
using Vec384fp6 = std::array<Vec384x, 3>;     // Fp6: 3 x Vec384x
using Vec384fp12 = std::array<Vec384fp6, 2>;  // Fp12: 2 x Vec384fp6

// Helper functions for C array compatibility
inline Vec384 make_vec384(const uint64_t arr[6]) {
    return Vec384(arr);
}

// Conversion helpers for nested types
inline Vec384x make_vec384x(const uint64_t arr[2][6]) {
    return Vec384x(arr);
}

inline Vec384fp6 make_vec384fp6(const uint64_t arr[3][2][6]) {
    Vec384fp6 result;
    for (size_t i = 0; i < 3; i++) {
        result[i] = Vec384x(&arr[i][0]);
    }
    return result;
}

inline Vec384fp12 make_vec384fp12(const uint64_t arr[2][3][2][6]) {
    Vec384fp12 result;
    for (size_t i = 0; i < 2; i++) {
        result[i] = make_vec384fp6(&arr[i][0][0]);
    }
    return result;
}

// Include fields.h to get the field operation functions
// This must be included after the Vec384 class definition
// to avoid circular dependencies
#include "fields.h"

// Inline implementations of Vec384 member functions
// These call the static inline functions from fields.h
// All operations are non-modifying and return new Vec384 instances
inline Vec384 Vec384::operator+(const Vec384& other) const {
    Vec384 result;
    add_fp(result.data_.data(), data_.data(), other.data_.data());
    return result;
}

inline Vec384 Vec384::operator-(const Vec384& other) const {
    Vec384 result;
    sub_fp(result.data_.data(), data_.data(), other.data_.data());
    return result;
}

inline Vec384 Vec384::operator*(const Vec384& other) const {
    Vec384 result;
    mul_fp(result.data_.data(), data_.data(), other.data_.data());
    return result;
}

inline Vec384 Vec384::sqr() const {
    Vec384 result;
    sqr_fp(result.data_.data(), data_.data());
    return result;
}

inline Vec384 Vec384::mul_by_3() const {
    Vec384 result;
    mul_by_3_fp(result.data_.data(), data_.data());
    return result;
}

inline Vec384 Vec384::mul_by_8() const {
    Vec384 result;
    mul_by_8_fp(result.data_.data(), data_.data());
    return result;
}

inline Vec384 Vec384::lshift(size_t count) const {
    Vec384 result;
    lshift_fp(result.data_.data(), data_.data(), count);
    return result;
}

inline Vec384 Vec384::rshift(size_t count) const {
    Vec384 result;
    rshift_fp(result.data_.data(), data_.data(), count);
    return result;
}

inline Vec384 Vec384::div_by_2() const {
    Vec384 result;
    div_by_2_fp(result.data_.data(), data_.data());
    return result;
}

inline Vec384 Vec384::cneg(uint64_t flag) const {
    Vec384 result;
    cneg_fp(result.data_.data(), data_.data(), flag);
    return result;
}

inline Vec384 Vec384::neg() const {
    return cneg(1);
}

inline Vec384 Vec384::from() const {
    Vec384 result;
    from_fp(result.data_.data(), data_.data());
    return result;
}

inline Vec384 Vec384::redc(const uint64_t a[12]) {
    Vec384 result;
    redc_fp(result.data_.data(), a);
    return result;
}

// Inline implementations of Vec384x member functions
// These call the static inline Fp2 functions from fields.h directly on the
// underlying 2x64-bit-limb storage via vec384x views.

inline Vec384x Vec384x::operator+(const Vec384x& other) const {
    Vec384x result;
    auto *r = reinterpret_cast<vec384x*>(result.data());
    const auto *a = reinterpret_cast<const vec384x*>(this->data());
    const auto *b = reinterpret_cast<const vec384x*>(other.data());
    add_fp2(*r, *a, *b);
    return result;
}

inline Vec384x Vec384x::operator-(const Vec384x& other) const {
    Vec384x result;
    auto *r = reinterpret_cast<vec384x*>(result.data());
    const auto *a = reinterpret_cast<const vec384x*>(this->data());
    const auto *b = reinterpret_cast<const vec384x*>(other.data());
    sub_fp2(*r, *a, *b);
    return result;
}

inline Vec384x Vec384x::mul_by_3() const {
    Vec384x result;
    auto *r = reinterpret_cast<vec384x*>(result.data());
    const auto *a = reinterpret_cast<const vec384x*>(this->data());
    mul_by_3_fp2(*r, *a);
    return result;
}

inline Vec384x Vec384x::mul_by_8() const {
    Vec384x result;
    auto *r = reinterpret_cast<vec384x*>(result.data());
    const auto *a = reinterpret_cast<const vec384x*>(this->data());
    mul_by_8_fp2(*r, *a);
    return result;
}

inline Vec384x Vec384x::lshift(size_t count) const {
    Vec384x result;
    auto *r = reinterpret_cast<vec384x*>(result.data());
    const auto *a = reinterpret_cast<const vec384x*>(this->data());
    lshift_fp2(*r, *a, count);
    return result;
}

inline Vec384x Vec384x::mul(const Vec384x& other) const {
    Vec384x result;
    auto *r = reinterpret_cast<vec384x*>(result.data());
    const auto *a = reinterpret_cast<const vec384x*>(this->data());
    const auto *b = reinterpret_cast<const vec384x*>(other.data());
    mul_fp2(*r, *a, *b);
    return result;
}

inline Vec384x Vec384x::sqr() const {
    Vec384x result;
    auto *r = reinterpret_cast<vec384x*>(result.data());
    const auto *a = reinterpret_cast<const vec384x*>(this->data());
    sqr_fp2(*r, *a);
    return result;
}

inline Vec384x Vec384x::cneg(uint64_t flag) const {
    Vec384x result;
    auto *r = reinterpret_cast<vec384x*>(result.data());
    const auto *a = reinterpret_cast<const vec384x*>(this->data());
    cneg_fp2(*r, *a, flag ? 1 : 0);
    return result;
}

#endif // __cplusplus

#endif // __VEC384_HPP__


