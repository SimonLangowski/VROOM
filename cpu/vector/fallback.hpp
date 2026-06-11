#pragma once
#include <array>
#include <stdint.h>
#include <cstring>
#include <cstdio>
#include <immintrin.h>
#include <cassert>
#include <type_traits>

// Fallback implementation of AVXVector using standard operations
// This provides the same interface as the AVX512 implementation but uses
// loops and standard operations instead of SIMD instructions

// Helper functions
inline constexpr int ceil_div(int n, int d) {
    return (n + d - 1) / d;
}

// Helper functions for 52-bit operations
constexpr uint64_t MASK_52_BITS = (1ULL << 52) - 1;

inline uint64_t apply_52bit_mask(uint64_t value) {
    return value & MASK_52_BITS;
}

// Scalar type for compatibility
typedef uint64_t AVXScalar;

inline AVXScalar Scalar(uint64_t v) {
    return v;
}

constexpr AVXScalar get_scalar_mask(int bits) {
    return ((uint64_t)(1) << bits) - 1;
}

inline AVXScalar scalar_shift(const AVXScalar &s, int bits) {
    return s >> bits;
}

// Fallback AVXVector implementation
template <int limbs>
class AVXVector {
public:
    static constexpr int LIMBS_PER_VEC = 8;  // Keep for compatibility
    static constexpr int VEC_LIMBS = (limbs + LIMBS_PER_VEC - 1) / LIMBS_PER_VEC;

    std::array<uint64_t, LIMBS_PER_VEC * VEC_LIMBS> data;

    AVXVector() = default;
    
    inline AVXVector(uint64_t v) {
        for (int i = 0; i < limbs; i++) {
            data[i] = v;
        }
    }
    
    inline AVXVector(const std::array<uint64_t, limbs> &arr) {
        for (int i = 0; i < limbs; i++) {
            data[i] = arr[i];
        }
    }

    /** Widen to W lanes; lanes limbs..W-1 must be 0 (maskz). */
    template <int W, typename = std::enable_if_t<(W > limbs) && (W <= 8) && (VEC_LIMBS == 1) && (AVXVector<W>::VEC_LIMBS == 1)>>
    AVXVector<W> widening_cast() const {
        AVXVector<W> out;
        __builtin_memcpy(out.data.data(), data.data(), 8 * sizeof(uint64_t));
        return out;
    }

    /** Narrow to W lanes; upper lanes of source are ignored. */
    template <int W, typename = std::enable_if_t<(W < limbs) && (W > 0) && (VEC_LIMBS == 1) && (AVXVector<W>::VEC_LIMBS == 1)>>
    AVXVector<W> narrowing_cast() const {
        AVXVector<W> out;
        for (int i = 0; i < W; ++i) {
            out.data[i] = data[i];
        }
        return out;
    }

    /** Lanes 0..N-1 from array; lanes N..7 zero (maskz). AVXVector<8> only. */
    template <int N, typename = std::enable_if_t<(limbs == 8) && (N > 0) && (N < 8)>>
    static AVXVector<8> widening_cast(const std::array<uint64_t, N> &narrow) {
        AVXVector<8> out;
        for (int i = 0; i < N; ++i) {
            out.data[i] = narrow[i];
        }
        return out;
    }

    inline constexpr AVXVector(const std::array<__m512i, VEC_LIMBS> &array) {
        union {
            __m512i m;
            uint64_t u[8];
        } u;
        for (int i = 0; i < VEC_LIMBS; i++) {
            u.m = array[i];
            for (int j = 0; j < 8; j++) {
                data[i * 8 + j] = u.u[j];
            }
        }
    }

    inline __m512i getm512() const {
        union {
            __m512i m;
            uint64_t u[8];
        } u = {};
        for (int i = 0; i < std::min(limbs, 8); i++) {
            u.u[i] = data[i];
        }
        return u.m;
    }

    // Multiplication operations
    inline AVXVector mulhi(const AVXVector &a, const AVXVector &b) const {
        AVXVector out;
        for (int i = 0; i < limbs; i++) {
            // Python: self.data[i] + (((a.data[i] & mask) * (b.data[i] & mask)) >> mulbits)
            // Mask inputs first, multiply (need 128-bit for 52-bit * 52-bit = 104-bit product), shift RIGHT, then add to data[i]
            uint64_t a_masked = a.data[i] & MASK_52_BITS;
            uint64_t b_masked = b.data[i] & MASK_52_BITS;
            __uint128_t product = (__uint128_t)a_masked * (__uint128_t)b_masked;
            uint64_t hi_part = (uint64_t)(product >> 52);  // Shift BEFORE adding
            out.data[i] = data[i] + hi_part;
        }
        return out;
    }

    inline AVXVector mullo(const AVXVector &a, const AVXVector &b) const {
        AVXVector out;
        for (int i = 0; i < limbs; i++) {
            // Python: self.data[i] + (((a.data[i] & mask) * (b.data[i] & mask)) & mask)
            // Mask inputs first, multiply (need 128-bit for 52-bit * 52-bit = 104-bit product), mask low bits, then add to data[i]
            uint64_t a_masked = a.data[i] & MASK_52_BITS;
            uint64_t b_masked = b.data[i] & MASK_52_BITS;
            __uint128_t product = (__uint128_t)a_masked * (__uint128_t)b_masked;
            uint64_t lo_part = (uint64_t)(product & MASK_52_BITS);  // Mask low 52 bits
            out.data[i] = data[i] + lo_part;
        }
        return out;
    }

    inline AVXVector mulhi_scalar(const AVXVector &a, const AVXScalar &s) const {
        AVXVector out;
        for (int i = 0; i < limbs; i++) {
            // Python: self.data[i] + (((a.data[i] & mask) * (s & mask)) >> mulbits)
            // Mask inputs first, multiply (need 128-bit for 52-bit * 52-bit = 104-bit product), shift RIGHT, then add to data[i]
            uint64_t a_masked = a.data[i] & MASK_52_BITS;
            uint64_t s_masked = s & MASK_52_BITS;
            __uint128_t product = (__uint128_t)a_masked * (__uint128_t)s_masked;
            uint64_t hi_part = (uint64_t)(product >> 52);  // Shift BEFORE adding
            out.data[i] = data[i] + hi_part;
        }
        return out;
    }

    inline AVXVector mullo_scalar(const AVXVector &a, const AVXScalar &s) const {
        AVXVector out;
        for (int i = 0; i < limbs; i++) {
            // Python: self.data[i] + (((a.data[i] & mask) * (s & mask)) & mask)
            // Mask inputs first, multiply (need 128-bit for 52-bit * 52-bit = 104-bit product), mask low bits, then add to data[i]
            uint64_t a_masked = a.data[i] & MASK_52_BITS;
            uint64_t s_masked = s & MASK_52_BITS;
            __uint128_t product = (__uint128_t)a_masked * (__uint128_t)s_masked;
            uint64_t lo_part = (uint64_t)(product & MASK_52_BITS);  // Mask low 52 bits
            out.data[i] = data[i] + lo_part;
        }
        return out;
    }

    inline AVXVector maskz_madd52lo(uint8_t mask, const AVXVector &b, const AVXVector &c) const {
        AVXVector out;
        for (int vi = 0; vi < VEC_LIMBS; vi++) {
            for (int j = 0; j < LIMBS_PER_VEC; j++) {
                const int idx = vi * LIMBS_PER_VEC + j;
                if (((mask >> j) & 1) == 0) {
                    out.data[idx] = 0;
                    continue;
                }
                const uint64_t z = idx < limbs ? data[idx] : 0;
                const uint64_t bv = idx < limbs ? b.data[idx] : 0;
                const uint64_t cv = idx < limbs ? c.data[idx] : 0;
                const uint64_t b_m = bv & MASK_52_BITS;
                const uint64_t c_m = cv & MASK_52_BITS;
                const __uint128_t product = (__uint128_t)b_m * (__uint128_t)c_m;
                const uint64_t lo_part = (uint64_t)(product & MASK_52_BITS);
                out.data[idx] = z + lo_part;
            }
        }
        return out;
    }

    inline AVXVector maskz_madd52hi(uint8_t mask, const AVXVector &b, const AVXVector &c) const {
        AVXVector out;
        for (int vi = 0; vi < VEC_LIMBS; vi++) {
            for (int j = 0; j < LIMBS_PER_VEC; j++) {
                const int idx = vi * LIMBS_PER_VEC + j;
                if (((mask >> j) & 1) == 0) {
                    out.data[idx] = 0;
                    continue;
                }
                const uint64_t z = idx < limbs ? data[idx] : 0;
                const uint64_t bv = idx < limbs ? b.data[idx] : 0;
                const uint64_t cv = idx < limbs ? c.data[idx] : 0;
                const uint64_t b_m = bv & MASK_52_BITS;
                const uint64_t c_m = cv & MASK_52_BITS;
                const __uint128_t product = (__uint128_t)b_m * (__uint128_t)c_m;
                const uint64_t hi_part = (uint64_t)(product >> 52);
                out.data[idx] = z + hi_part;
            }
        }
        return out;
    }

    // Shift operations
    inline AVXVector srli(const int bits) const {
        AVXVector out;
        for (int i = 0; i < limbs; i++) {
            out.data[i] = data[i] >> bits;
        }
        return out;
    }

    inline AVXVector slli(const int bits) const {
        AVXVector out;
        for (int i = 0; i < limbs; i++) {
            out.data[i] = data[i] << bits;
        }
        return out;
    }

    // Arithmetic operations
    inline AVXVector add(const AVXVector &other) const {
        AVXVector out;
        for (int i = 0; i < limbs; i++) {
            out.data[i] = data[i] + other.data[i];
        }
        return out;
    }

    inline AVXVector sub(const AVXVector &other) const {
        AVXVector out;
        for (int i = 0; i < limbs; i++) {
            out.data[i] = data[i] - other.data[i];
        }
        return out;
    }

    inline AVXVector sub_min(const AVXVector &other) const {
        AVXVector out;
        for (int i = 0; i < limbs; i++) {
            out.data[i] = std::min(data[i], data[i] - other.data[i]);
        }
        return out;
    }

    inline AVXVector add_min(const AVXVector &other) const {
        AVXVector out;
        for (int i = 0; i < limbs; i++) {
            out.data[i] = std::min(data[i], data[i] + other.data[i]);
        }
        return out;
    }

    inline AVXVector min(const AVXVector &other) const {
        AVXVector out;
        for (int i = 0; i < limbs; i++) {
            out.data[i] = std::min(data[i], other.data[i]);
        }
        return out;
    }

    // Bitwise operations
    inline AVXVector and_scalar(const AVXScalar &other) const {
        AVXVector out;
        for (int i = 0; i < limbs; i++) {
            out.data[i] = data[i] & other;
        }
        return out;
    }

    inline AVXVector bit_or(const AVXVector &other) const {
        AVXVector out;
        for (int i = 0; i < limbs; i++) {
            out.data[i] = data[i] | other.data[i];
        }
        return out;
    }

    inline AVXVector ctime_select(const bool selection, const AVXVector &b) const {
        // fallback, not worried about constant time, just correctness
        return selection ? *this : b;
    }

    inline AVXVector permute(AVXVector<8> perm) const {
        AVXVector out;
        for (int j = 0; j < limbs; j++) {
            const uint64_t to = perm.data[j];
            assert(to < 8);
            out.data[j] = data[to];
        }
        return out;
    }

    // Access operations
    inline uint64_t extract0() const {
        return data[0];
    }

    /** Return the first lane as AVXScalar (use when all lanes are equal, e.g. after horizontal sum). */
    inline AVXScalar as_scalar() const {
        return data[0];
    }

    inline void store(uint64_t* dst) const {
        for (int i = 0; i < limbs; i++) {
            dst[i] = data[i];
        }
    }

    inline void load(uint64_t* src) {
        for (int i = 0; i < limbs; i++) {
            data[i] = src[i];
        }
    }

    inline void print(const char* label) const {
        printf("%s: [", label);
        for (int i = 0; i < limbs; i++) {
            printf("%lu ", data[i]);
        }
        printf("]\n");
    }

    // Additional helper methods for debugging and compatibility
    inline uint64_t get_limb(int index) const {
        return data[index];
    }

    inline void set_limb(int index, uint64_t value) {
        data[index] = value;
    }

    // Apply 52-bit mask to all limbs
    inline void apply_mask() {
        for (int i = 0; i < limbs; i++) {
            data[i] = apply_52bit_mask(data[i]);
        }
    }

    inline std::array<uint64_t, limbs> to_array() const {
        std::array<uint64_t, limbs> result{};
        for (int i = 0; i < limbs; i++) {
            result[i] = data[i];
        }
        return result;
    }
};

/** Scalar pair matching the low 128 bits of the first ZMM group (fallback for ShortVector / __m128i). */
class ShortVector {
    uint64_t q0_;
    uint64_t q1_;

public:
    ShortVector() : q0_(0), q1_(0) {}
    ShortVector(uint64_t q0, uint64_t q1) : q0_(q0), q1_(q1) {}

    template <int limbs>
    static ShortVector from_avx_lo(const AVXVector<limbs> &a) {
        return ShortVector(a.get_limb(0), a.get_limb(1));
    }

    ShortVector add(const ShortVector &o) const {
        return ShortVector(q0_ + o.q0_, q1_ + o.q1_);
    }

    ShortVector srli(int bits) const {
        return ShortVector(q0_ >> bits, q1_ >> bits);
    }

    template <int limbs>
    AVXVector<limbs> broadcast_to_avx() const {
        AVXVector<limbs> out;
        const int span = AVXVector<limbs>::LIMBS_PER_VEC * AVXVector<limbs>::VEC_LIMBS;
        for (int i = 0; i < span; i++) {
            out.set_limb(i, q0_);
        }
        return out;
    }

    uint64_t q0() const { return q0_; }
    uint64_t q1() const { return q1_; }
};
