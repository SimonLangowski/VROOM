#pragma once
#include <stdint.h>
// [0, 2**32 * p) -> [0, 2p)
// since 2p * 2p = 4p^2 this limits the moduli to 30 bits.
__device__ __forceinline__ uint32_t montgomery_reduce30(const uint64_t u, const uint32_t p, const uint32_t p_inv) {  
    // Since the input is already u64, this is just 2 multiplication operations
    uint32_t u_0 = u;
    uint32_t u_1 = u >> 32;
    uint32_t q = u_0 * p_inv; // mul
    uint32_t h = __umulhi(q, p); // mulwide
    return u_1 - h + p; // iadd3
}

// [0, 2**32 * p) -> [0, p)
// Effectively the same as friendly reduction but works for any 32 bit modulus
__device__ __forceinline__ uint32_t montgomery_reduce32(const uint64_t u, const uint32_t p, const uint32_t p_inv) {  
    // Since the input is already u64, this is just 2 multiplication operations
    uint32_t u_0 = u;
    uint32_t u_1 = u >> 32;
    uint32_t q = u_0 * p_inv; // mul
    uint32_t h = __umulhi(q, p); // mulwide
    return u_1 - h + (h > u_1) * p; // iadd3, cmp, selp
}


// __device__ __forceinline__ uint32_t montgomery_add32(const uint32_t a, const uint32_t b, const uint32_t modulus) {
//     uint64_t sum = (uint64_t)(a) + b;
//     uint32_t low = sum;
//     // carry = 0 or 1 and we can use as conditional correction bit
//     uint32_t carry = sum >> 32;
//     // if carry is 0, this does nothing, and we are done <- What if it is between modulus and 2^32, and then two of those are added together, and now carry = 1 is not guaranteed to be < 2^32?
//     // if carry is 1, then we have to output (2^32) + low - modulus
//     // We do this as a low operation by knowing that the extra bit must always cancel
//     // ((2^32) + low - modulus) % 2^32 = low - modulus
//     return low - carry * modulus;
// }

// [0, 2p) + [0, 2p) -> [0, 2p)
__device__ __forceinline__ uint32_t montgomery_add30(const uint32_t a, const uint32_t b, const uint32_t two_modulus) {
    // for 30 bits, 2p is 31 bits, and this sum fits in a 32 bit word
    uint32_t sum = a + b;
    // do conditional correction
    sum -= (sum >= two_modulus) * two_modulus;
    return sum;
}

// [0, 2p) + [0, 2p) -> [0, 2p)
__device__ __forceinline__ uint32_t montgomery_sub30(const uint32_t a, const uint32_t b, const uint32_t two_modulus) {
    // two_modulus - b is also in [0, 2p)
    uint32_t sum = a + two_modulus - b; // iadd3
    // do conditional correction
    sum -= (sum >= two_modulus) * two_modulus;
    return sum;
}

// Modular reduction montgmery without corection
template<template <int, int, int> class FragmentType, int ModCount, int TPI>
class MontRedundant30 {
    private:
    using ModType = FragmentType<1, ModCount, TPI>;
    const ModType mod_info;
    const ModType mont_info;

    public:
    class MontRedundant30Info {
        public:
        uint32_t modulus;
        uint32_t mod_factor;
    };
    static constexpr int COLS = ModCount;
    static constexpr int ROW_LIMBS = 2;
    static constexpr int MODBITS = 30;

    __device__ MontRedundant30(const uint32_t* mod_data) : mod_info(mod_data, 0), mont_info(mod_data, padded(ModCount, TPI)) {}

    __device__ __forceinline__  MontRedundant30Info load_modulus(const int j) const {
        MontRedundant30Info m = {mod_info.get(0, j), mont_info.get(0, j)};
        return m;
    }

    __device__ __forceinline__ uint32_t mod_reduce(uint64_t unreduced, const MontRedundant30Info &m) const {
        return montgomery_reduce30(unreduced, m.modulus, m.mod_factor);
    }
};

// Modular reduction for psuedomersenne moduli
template<template <int, int, int> class FragmentType, int ModCount, int TPI>
class Mont32 {
    private:
    using ModType = FragmentType<1, ModCount, TPI>;
    const ModType mod_info;
    const ModType mont_info;

    public:
    class Mont32Info {
        public:
        uint32_t modulus;
        uint32_t mod_factor;
    };
    static constexpr int COLS = ModCount;
    static constexpr int ROW_LIMBS = 2;
    static constexpr int MODBITS = 32;

    __device__ Mont32(const uint32_t* mod_data) : mod_info(mod_data, 0), mont_info(mod_data, padded(ModCount, TPI)) {}

    __device__ __forceinline__  Mont32Info load_modulus(const int j) const {
        Mont32Info m = {mod_info.get(0, j), mont_info.get(0, j)};
        return m;
    }

    __device__ __forceinline__ uint32_t mod_reduce(uint64_t unreduced, const Mont32Info &m) const {
        return montgomery_reduce32(unreduced, m.modulus, m.mod_factor);
    }
};