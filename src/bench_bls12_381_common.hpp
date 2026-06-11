#pragma once

#include <benchmark/benchmark.h>
#include <array>
#include <cstdint>
#include "bls12_381_ring_params.hpp"
#include "fp12.hpp"
#include "fp2.hpp"
#include "final_exponentiation.hpp"
#include "miller.hpp"
#include "inversion.hpp"
#include "conversion_inversion.hpp"
#include "pairing.hpp"
#include "ec.hpp"
#include "scalar_mult.hpp"
#include "../cpu/precompute/gmp_wrapper.hpp"

namespace bench_bls12_381 {

inline BigInt field_modulus() {
    // Avoid BLST global BLS12_381_P (vec384) name clash from conversion_inversion.hpp.
    return BigInt(bls12_381_params::BLS12_381_P);
}

// G1 and G2 generator coordinates (BLS12-381, decimal strings).
inline const char *g1_gen_x =
    "3685416753713387016781088315183077757961620795782546409894578378688607592378376318836054947676345821548104185464507";
inline const char *g1_gen_y =
    "1339506544944476473020471379941921221584933875938349620426543736416511423956333506472724655353366534992391756441569";
inline const char *g2_gen_x_c0 =
    "352701069587466618187139116011060144890029952792775240219908644239793785735715026873347600343865175952761926303160";
inline const char *g2_gen_x_c1 =
    "3059144344244213709971259814753781636986470325476647558659373206291635324768958432433509563104347017837885763365758";
inline const char *g2_gen_y_c0 =
    "1985150602287291935568054521177171638300868978215655730859378665066344726373823718423869104263333984641494340347905";
inline const char *g2_gen_y_c1 =
    "927553665492332455747201965776037880757740193453592970025027978793976877002675564980949289727957565575433344219582";
inline const char *scalar_modulus_hex =
    "73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001";

inline void bigint_to_bytes_le(uint8_t *bytes, const BigInt &value, size_t len) {
    BigInt temp = value;
    for (size_t i = 0; i < len; i++) {
        bytes[i] = static_cast<uint8_t>(temp.to_ulong() & 0xff);
        temp = temp >> 8;
    }
}

template<class Ring>
inline typename FP12Ring<Ring>::StandardElement generate_random_fp12(const Ring &ring, uint64_t seed) {
    typename FP12Ring<Ring>::StandardElement result;
    BigInt p = field_modulus();
    volatile uint64_t s = seed;
    gmp_randclass rng(gmp_randinit_default);
    rng.seed(s);
    for (int i = 0; i < 12; i++) {
        result[i] = ring.from_bigint(BigInt(rng.get_z_range(p.get_mpz())));
    }
    return result;
}

template<class Ring>
inline typename G1<Ring>::ProjPoint generate_random_g1_projective(const Ring &ring, uint64_t seed) {
    BigInt p = field_modulus();
    gmp_randclass rng(gmp_randinit_default);
    rng.seed(seed);
    return typename G1<Ring>::ProjPoint(
        ring.from_bigint(BigInt(rng.get_z_range(p.get_mpz()))),
        ring.from_bigint(BigInt(rng.get_z_range(p.get_mpz()))),
        ring.one());
}

template<class Ring>
inline AffinePoint<typename Ring::StandardElement> generate_random_g1_affine(const Ring &ring, uint64_t seed) {
    BigInt p = field_modulus();
    gmp_randclass rng(gmp_randinit_default);
    rng.seed(seed);
    return AffinePoint<typename Ring::StandardElement>{
        ring.from_bigint(BigInt(rng.get_z_range(p.get_mpz()))),
        ring.from_bigint(BigInt(rng.get_z_range(p.get_mpz())))};
}

template<class Ring>
inline typename G2<Ring>::ProjPoint generate_random_g2_projective(
    const FP2Ring<Ring> &fp2ring, const Ring &ring, uint64_t seed) {
    BigInt p = field_modulus();
    gmp_randclass rng(gmp_randinit_default);
    rng.seed(seed);
    return typename G2<Ring>::ProjPoint(
        typename FP2Ring<Ring>::StandardElement(
            ring.from_bigint(BigInt(rng.get_z_range(p.get_mpz()))),
            ring.from_bigint(BigInt(rng.get_z_range(p.get_mpz())))),
        typename FP2Ring<Ring>::StandardElement(
            ring.from_bigint(BigInt(rng.get_z_range(p.get_mpz()))),
            ring.from_bigint(BigInt(rng.get_z_range(p.get_mpz())))),
        fp2ring.one());
}

template<class Ring>
inline AffinePoint<typename FP2Ring<Ring>::StandardElement> generate_random_g2_affine(const Ring &ring, uint64_t seed) {
    BigInt p = field_modulus();
    gmp_randclass rng(gmp_randinit_default);
    rng.seed(seed);
    return AffinePoint<typename FP2Ring<Ring>::StandardElement>{
        typename FP2Ring<Ring>::StandardElement(
            ring.from_bigint(BigInt(rng.get_z_range(p.get_mpz()))),
            ring.from_bigint(BigInt(rng.get_z_range(p.get_mpz())))),
        typename FP2Ring<Ring>::StandardElement(
            ring.from_bigint(BigInt(rng.get_z_range(p.get_mpz()))),
            ring.from_bigint(BigInt(rng.get_z_range(p.get_mpz()))))};
}

template<class Ring>
inline void setup_g1_g2_generators(
    const Ring &ring,
    AffinePoint<typename Ring::StandardElement> &e1,
    AffinePoint<typename FP2Ring<Ring>::StandardElement> &e2) {
    e1 = {
        ring.from_bigint(BigInt(g1_gen_x, 10)),
        ring.from_bigint(BigInt(g1_gen_y, 10))};
    e2 = {
        typename FP2Ring<Ring>::StandardElement(
            ring.from_bigint(BigInt(g2_gen_x_c0, 10)),
            ring.from_bigint(BigInt(g2_gen_x_c1, 10))),
        typename FP2Ring<Ring>::StandardElement(
            ring.from_bigint(BigInt(g2_gen_y_c0, 10)),
            ring.from_bigint(BigInt(g2_gen_y_c1, 10)))};
}

} // namespace bench_bls12_381
