#pragma once
#include "../cpu/vector/vector_impl.hpp"
#include "../cpu/precompute/gmp_wrapper.hpp"
#include <array>
#include <stdint.h>
#include <cstring>
#include <string>
#include "constexpr_utils.hpp"
#include "../blst/vect.h"
#include "../blst/fields.h"

// Forward declaration - function is in blst/recip.c
extern "C" {
    void ct_inverse_mod_384_wrapper(vec384 out, const vec384 inp);
    void blst_fp12_inverse(vec384fp12 ret, const vec384fp12 a);
}
#include "../blst/consts.h"
#include "fp12.hpp"

// Pack vec384 (6×64-bit LE) into convert_to digit width (50 on MatrixNoK, 52 otherwise).
// Fully unrolled — no runtime division.
template<int word_bits, int ndigits>
inline void vec384_to_convert_digits(const vec384 v, std::array<uint64_t, ndigits> &digits) {
    static_assert(word_bits == 50 || word_bits == 52);
    static_assert(ndigits == 8);
    if constexpr (word_bits == 52) {
        constexpr uint64_t MASK = (1ULL << 52) - 1;
        digits[0] = v[0] & MASK;
        digits[1] = ((v[0] >> 52) | (v[1] << 12)) & MASK;
        digits[2] = ((v[1] >> 40) | (v[2] << 24)) & MASK;
        digits[3] = ((v[2] >> 28) | (v[3] << 36)) & MASK;
        digits[4] = ((v[3] >> 16) | (v[4] << 48)) & MASK;
        digits[5] = (v[4] >> 4) & MASK;
        digits[6] = ((v[4] >> 56) | (v[5] << 8)) & MASK;
        digits[7] = (v[5] >> 44) & MASK;
    } else {
        constexpr uint64_t MASK = (1ULL << 50) - 1;
        digits[0] = v[0] & MASK;
        digits[1] = ((v[0] >> 50) | (v[1] << 14)) & MASK;
        digits[2] = ((v[1] >> 36) | (v[2] << 28)) & MASK;
        digits[3] = ((v[2] >> 22) | (v[3] << 42)) & MASK;
        digits[4] = (v[3] >> 8) & MASK;
        digits[5] = ((v[3] >> 58) | (v[4] << 6)) & MASK;
        digits[6] = ((v[4] >> 44) | (v[5] << 20)) & MASK;
        digits[7] = (v[5] >> 30) & MASK;
    }
}

template<int limbs>
AVXVector<limbs> to_vec(uint8_t *le_buf) {
    static_assert(limbs == 8, "to_vec currently only supports 8 limbs");
    std::array<uint64_t, 8> digits{};
    vec384_to_convert_digits<52, 8>(*(const vec384 *)le_buf, digits);
    return AVXVector<limbs>(digits);
}

// rns_reduce_raw (lo, hi) -> vec768 integer N before blst_prepare.
// Combine lo/hi as 64-bit lanes, then accumulate at 52-bit offsets LSB-first.
template<int digit_idx, int bits_per_digit=52>
inline void accumulate_digit(std::array<uint64_t, 8> &sum_array, uint64_t digit) {
    constexpr int offset = digit_idx * bits_per_digit;
    constexpr int sum_word = offset / 64;
    constexpr int sum_bit = offset % 64;
    __uint128_t tmp = (__uint128_t)(sum_array[sum_word]) + ((__uint128_t)(digit) << sum_bit);
    sum_array[sum_word] = (uint64_t)tmp;
    uint64_t c = (uint64_t)(tmp >> 64);
    static_assert(sum_word + 1 < 8, "Accumulated digit overflow");
    sum_array[sum_word + 1] = c;
}

inline void lo_hi_to_vec768(const std::array<uint64_t, 8> &digits_lo,
    const std::array<uint64_t, 8> &digits_hi, vec768 &sum) {
    std::array<uint64_t, 9> digits{};
    digits[0] = digits_lo[0];
    for (int i = 1; i < 8; i++) {
        digits[i] = digits_hi[i - 1] + digits_lo[i];
    }
    digits[8] = digits_hi[7];
    std::array<uint64_t, 8> sum_array{};
    sum_array[0] = digits[0];
    accumulate_digit<1>(sum_array, digits[1]);
    accumulate_digit<2>(sum_array, digits[2]);
    accumulate_digit<3>(sum_array, digits[3]);
    accumulate_digit<4>(sum_array, digits[4]);
    accumulate_digit<5>(sum_array, digits[5]);
    accumulate_digit<6>(sum_array, digits[6]);
    accumulate_digit<7>(sum_array, digits[7]);
    accumulate_digit<8>(sum_array, digits[8]);
    for (int i = 0; i < 8; i++) {
        sum[i] = sum_array[i];
    }
    for (int i = 8; i < 12; i++) {
        sum[i] = 0;
    }
}

template<bool montgomery=false>
void blst_prepare(vec768 &sum, vec384 ret) {
    // Always: redc_fp -> multiply by RR -> redc_fp (gives normal form)
    vec384 temp;
    redc_fp(temp, sum);  // temp = sum * R^-1 mod p
    vec768 prod;
    mul_384(prod, temp, BLS12_381_RR);  // prod = sum * R
    redc_fp(ret, prod);  // ret = sum (normal form)

    if (montgomery) {
        // For Montgomery form: multiply by RR again
        vec768 prod2;
        mul_384(prod2, ret, BLS12_381_RR);  // prod2 = sum * R^2
        redc_fp(ret, prod2);  // ret = sum * R (Montgomery form)
    }
}

template<bool montgomery=false>
void from_vec(AVXVector<8> &vec_lo, AVXVector<8> &vec_hi, vec384 ret) {
    std::array<uint64_t, 8> digits_lo;
    std::array<uint64_t, 8> digits_hi;
    vec_hi.store(digits_hi.data());
    vec_lo.store(digits_lo.data());
    vec768 sum;
    lo_hi_to_vec768(digits_lo, digits_hi, sum);
    blst_prepare<montgomery>(sum, ret);
}

template<class Ring, bool montgomery=false>
void convert_to_int(const Ring & ring, typename Ring::StandardElement & a, vec384 ret) {
    AVXVector<8> out_hi = AVXVector<8>(0);
    AVXVector<8> out_lo = AVXVector<8>(0);
    ring.sys.convert_from.conversion_matrix.rns_reduce_raw(a.m2.data, out_hi, out_lo);
    from_vec<montgomery>(out_lo, out_hi, ret);
}

template<class Ring, bool montgomery=false>
auto convert_to_rns(const Ring & ring, const vec384 ret) {
    vec384 temp_ret;
    if constexpr (montgomery) {
        from_fp(temp_ret, ret);
    } else {
        vec_copy(temp_ret, ret, sizeof(vec384));
    }
    // MatrixNoK convert_to uses 50-bit digits (CONVERT_TO_WORD_BITS); to_vec is 52-bit.
    constexpr int wb = Ring::CONVERT_TO_WORD_BITS;
    constexpr int ndigits = 8;
    std::array<uint64_t, ndigits> digits{};
    vec384_to_convert_digits<wb, ndigits>(temp_ret, digits);
    AVXVector<Ring::LIMBS> l;
    l.load(digits.data());
    return ring.template convert_to_rns_batch<1>(std::array<AVXVector<Ring::LIMBS>, 1>{l})[0];
}

template<class Ring>
auto invert_via_blst(const Ring & ring, typename Ring::StandardElement & a) {
    vec384 ret;
    convert_to_int<Ring, false>(ring, a, ret);
    ct_inverse_mod_384_wrapper(ret, ret);
    return convert_to_rns<Ring, false>(ring, ret);
}

template<class Ring>
class BLSTInverter {
    public:
    using StandardElement = typename Ring::StandardElement;

    StandardElement invert(const StandardElement &a, const Ring &ring) const {
        return invert_via_blst(ring, const_cast<StandardElement&>(a));
    }
};

// Convert FP12 (12 ring elements) from ring form to BLST Montgomery form,
// call blst_fp12_inverse, convert back to ring form
template<class Ring>
std::array<typename Ring::StandardElement, 12> invert_fp12_via_blst(
    const Ring& ring,
    const std::array<typename Ring::StandardElement, 12>& fp12_ring
) {
    // Step 1: Convert ring form to BLST Montgomery form using convert_to_int
    vec384fp12 fp12_mont;
    // vec384fp12 is vec384fp6[2], vec384fp6 is vec384x[3], vec384x is vec384[2]
    // For index i (0-11): fp12_mont[i/6][(i%6)/2][i%2]
    for (int i = 0; i < 12; i++) {
        convert_to_int<Ring, true>(ring, const_cast<typename Ring::StandardElement&>(fp12_ring[i]), 
                       fp12_mont[i/6][(i%6)/2][i%2]);
    }
    
    // Step 2: Call blst_fp12_inverse
    vec384fp12 fp12_inv_mont;
    blst_fp12_inverse(fp12_inv_mont, fp12_mont);
    
    // Step 3: Convert BLST Montgomery form back to ring form using convert_to_rns
    std::array<typename Ring::StandardElement, 12> fp12_inv_ring;
    for (int i = 0; i < 12; i++) {
        fp12_inv_ring[i] = convert_to_rns<Ring, true>(ring, fp12_inv_mont[i/6][(i%6)/2][i%2]);
    }
    
    return fp12_inv_ring;
}

// --- BigInt-based FP12 inversion via BLST ---
// This path converts ring -> BigInt (normal form) -> vec384 Montgomery (BLST),
// calls blst_fp12_inverse, then converts back to BigInt (normal form) -> ring.
static inline BigInt vec384_normal_to_bigint(const vec384 a) {
    BigInt result(0);
    BigInt two_to_64 = BigInt(1) << 64;
    for (int i = 5; i >= 0; i--) {
        result = result * two_to_64 + BigInt(static_cast<unsigned long>(a[i]));
    }
    return result;
}

static inline BigInt vec384_montgomery_to_bigint(const vec384 a) {
    vec384 normal;
    from_fp(normal, a);
    return vec384_normal_to_bigint(normal);
}

static inline void bigint_to_vec384_normal(vec384 out, const BigInt& a) {
    // Assumes `a` is already in normal form and non-negative (e.g. output of ring.to_bigint()).
    std::string hex_str = a.to_string(16);
    while (hex_str.length() < 96) {
        hex_str = "0" + hex_str;
    }

    memset(out, 0, sizeof(vec384));
    for (int i = 0; i < 6; i++) {
        std::string limb_hex = hex_str.substr((5 - i) * 16, 16);
        out[i] = std::stoull(limb_hex, nullptr, 16);
    }
}

static inline void bigint_to_vec384_montgomery(vec384 out, const BigInt& a) {
    vec384 normal;
    bigint_to_vec384_normal(normal, a);

    // out = normal * R^2 mod p
    vec768 prod;
    mul_384(prod, normal, BLS12_381_RR);
    redc_fp(out, prod);
}

template<class Ring>
std::array<typename Ring::StandardElement, 12> invert_fp12_via_blst_bigint(
    const Ring& ring,
    const std::array<typename Ring::StandardElement, 12>& fp12_ring
) {
    vec384fp12 fp12_mont;
    for (int i = 0; i < 12; i++) {
        BigInt v = ring.to_bigint(fp12_ring[i]);
        bigint_to_vec384_montgomery(fp12_mont[i/6][(i%6)/2][i%2], v);
    }

    vec384fp12 fp12_inv_mont;
    blst_fp12_inverse(fp12_inv_mont, fp12_mont);

    std::array<BigInt, 12> inv_bigints;
    for (int i = 0; i < 12; i++) {
        inv_bigints[i] = vec384_montgomery_to_bigint(fp12_inv_mont[i/6][(i%6)/2][i%2]);
    }

    return ring.template from_bigint_batch<12>(inv_bigints);
}

// FP12 inversion provider using invert_fp12_via_blst_bigint
template<class Ring>
class FP12BLSTInverter {
public:
    using StandardElement = std::array<typename Ring::StandardElement, 12>;
    using FP12R = FP12Ring<Ring>;
    StandardElement invert(const StandardElement& f, const FP12R& /* fp12ring */, const Ring& ring) const {
        return invert_fp12_via_blst(ring, f);
    }
};

template<class Ring>
class FP12RNS_BLSTInverter {
public:
    using FP12R = FP12Ring<Ring>;
    using StandardElement = typename FP12R::StandardElement;

    using Inverter = BLSTInverter<Ring>;
    Inverter inverter;

    StandardElement invert(const StandardElement& f, const FP12R& fp12ring, const Ring& ring) const {
        return fp12ring.template fp12_flat_inverse<Inverter>(f, inverter, ring);
    }
};