#pragma once

#include "bench_bls12_381_common.hpp"

namespace bench_bls12_381 {

template<class Ring>
static void BM_ModMul_Impl(benchmark::State& state) {
    Ring ring(field_modulus());
    gmp_randclass rng(gmp_randinit_default);
    rng.seed(42);
    BigInt p = field_modulus();
    auto a = ring.from_bigint(BigInt(rng.get_z_range(p.get_mpz())));
    auto b = ring.from_bigint(BigInt(rng.get_z_range(p.get_mpz())));
    benchmark::DoNotOptimize(a);
    benchmark::DoNotOptimize(b);
    for (auto _ : state) {
        auto result = ring.modmul(a, b);
        benchmark::DoNotOptimize(result);
    }
}

template<int Batch, class Ring>
static void BM_BatchModMul_Impl(benchmark::State& state) {
    Ring ring(field_modulus());
    gmp_randclass rng(gmp_randinit_default);
    rng.seed(42);
    BigInt p = field_modulus();
    std::array<typename Ring::StandardElement, Batch> a_elements;
    std::array<typename Ring::StandardElement, Batch> b_elements;
    for (int i = 0; i < Batch; i++) {
        a_elements[i] = ring.from_bigint(BigInt(rng.get_z_range(p.get_mpz())));
        b_elements[i] = ring.from_bigint(BigInt(rng.get_z_range(p.get_mpz())));
    }
    benchmark::DoNotOptimize(a_elements);
    benchmark::DoNotOptimize(b_elements);
    for (auto _ : state) {
        auto result = ring.template batch_modmul<Batch>(a_elements, b_elements);
        benchmark::DoNotOptimize(result);
    }
}

template<int Batch, class Ring>
static void BM_BatchReduceExpand_Impl(benchmark::State& state) {
    Ring ring(field_modulus());
    gmp_randclass rng(gmp_randinit_default);
    rng.seed(42);
    BigInt p = field_modulus();
    std::array<typename Ring::StandardElement, Batch> a_elements;
    std::array<typename Ring::StandardElement, Batch> b_elements;
    for (int i = 0; i < Batch; i++) {
        a_elements[i] = ring.from_bigint(BigInt(rng.get_z_range(p.get_mpz())));
        b_elements[i] = ring.from_bigint(BigInt(rng.get_z_range(p.get_mpz())));
    }
    std::array<typename Ring::template ReadyElement<Ring::MAX_ADD>, Batch> ready_elements;
    for (int i = 0; i < Batch; i++) {
        ready_elements[i] = ring.template ready<Ring::MAX_ADD>(a_elements[i] * b_elements[i]);
    }
    benchmark::DoNotOptimize(ready_elements);
    for (auto _ : state) {
        auto result = ring.template batch_reduce_expand<Batch>(ready_elements);
        benchmark::DoNotOptimize(result);
    }
}

template<class Ring>
static void BM_FP2_Mul_Impl(benchmark::State& state) {
    Ring ring(field_modulus());
    FP2Ring<Ring> fp2ring(ring);
    gmp_randclass rng(gmp_randinit_default);
    rng.seed(42);
    BigInt p = field_modulus();
    auto a = typename FP2Ring<Ring>::StandardElement(
        ring.from_bigint(BigInt(rng.get_z_range(p.get_mpz()))),
        ring.from_bigint(BigInt(rng.get_z_range(p.get_mpz()))));
    auto b = typename FP2Ring<Ring>::StandardElement(
        ring.from_bigint(BigInt(rng.get_z_range(p.get_mpz()))),
        ring.from_bigint(BigInt(rng.get_z_range(p.get_mpz()))));
    benchmark::DoNotOptimize(a);
    benchmark::DoNotOptimize(b);
    for (auto _ : state) {
        auto result = fp2ring.mul(a, b);
        benchmark::DoNotOptimize(result);
    }
}

template<class Ring>
static void BM_FP2_Sqr_Impl(benchmark::State& state) {
    Ring ring(field_modulus());
    FP2Ring<Ring> fp2ring(ring);
    gmp_randclass rng(gmp_randinit_default);
    rng.seed(42);
    BigInt p = field_modulus();
    auto a = typename FP2Ring<Ring>::StandardElement(
        ring.from_bigint(BigInt(rng.get_z_range(p.get_mpz()))),
        ring.from_bigint(BigInt(rng.get_z_range(p.get_mpz()))));
    benchmark::DoNotOptimize(a);
    for (auto _ : state) {
        auto result = fp2ring.square(a);
        benchmark::DoNotOptimize(result);
    }
}

template<class Ring>
static void BM_FP12_Mul_Impl(benchmark::State& state) {
    Ring ring(field_modulus());
    FP12Ring<Ring> fp12ring(ring);
    volatile uint64_t seed_a = 42;
    volatile uint64_t seed_b = 123;
    auto a = generate_random_fp12<Ring>(ring, seed_a);
    auto b = generate_random_fp12<Ring>(ring, seed_b);
    benchmark::DoNotOptimize(a);
    benchmark::DoNotOptimize(b);
    for (auto _ : state) {
        auto result = fp12ring.fp12_flat_mul(a, b, ring);
        benchmark::DoNotOptimize(result);
    }
}

template<class Ring>
static void BM_FP12_Sqr_Impl(benchmark::State& state) {
    Ring ring(field_modulus());
    FP12Ring<Ring> fp12ring(ring);
    volatile uint64_t seed = 42;
    auto a = generate_random_fp12<Ring>(ring, seed);
    benchmark::DoNotOptimize(a);
    for (auto _ : state) {
        auto result = fp12ring.fp12_flat_sqr(a, ring);
        benchmark::DoNotOptimize(result);
    }
}

template<class Ring>
static void BM_FP12_CyclotomicSqr_Impl(benchmark::State& state) {
    Ring ring(field_modulus());
    FP12Ring<Ring> fp12ring(ring);
    volatile uint64_t seed = 42;
    auto a = generate_random_fp12<Ring>(ring, seed);
    benchmark::DoNotOptimize(a);
    for (auto _ : state) {
        auto result = fp12ring.fp12_flat_cyclotomic_sqr(a, ring);
        benchmark::DoNotOptimize(result);
    }
}

template<class Ring>
static void BM_FP12_Conjugate_Impl(benchmark::State& state) {
    Ring ring(field_modulus());
    FP12Ring<Ring> fp12ring(ring);
    volatile uint64_t seed = 42;
    auto a = generate_random_fp12<Ring>(ring, seed);
    benchmark::DoNotOptimize(a);
    for (auto _ : state) {
        auto result = fp12ring.fp12_flat_conjugate(a, ring);
        benchmark::DoNotOptimize(result);
    }
}

template<class Ring>
static void BM_FP12_FrobeniusMap1_Impl(benchmark::State& state) {
    Ring ring(field_modulus());
    FP12Ring<Ring> fp12ring(ring);
    volatile uint64_t seed = 42;
    auto a = generate_random_fp12<Ring>(ring, seed);
    benchmark::DoNotOptimize(a);
    for (auto _ : state) {
        auto result = fp12ring.fp12_flat_frobenius_map1(a, ring);
        benchmark::DoNotOptimize(result);
    }
}

template<class Ring>
static void BM_FP12_FrobeniusMap2_Impl(benchmark::State& state) {
    Ring ring(field_modulus());
    FP12Ring<Ring> fp12ring(ring);
    volatile uint64_t seed = 42;
    auto a = generate_random_fp12<Ring>(ring, seed);
    benchmark::DoNotOptimize(a);
    for (auto _ : state) {
        auto result = fp12ring.fp12_flat_frobenius_map2(a, ring);
        benchmark::DoNotOptimize(result);
    }
}

template<class Ring>
static void BM_FP12_FrobeniusMap3_Impl(benchmark::State& state) {
    Ring ring(field_modulus());
    FP12Ring<Ring> fp12ring(ring);
    volatile uint64_t seed = 42;
    auto a = generate_random_fp12<Ring>(ring, seed);
    benchmark::DoNotOptimize(a);
    for (auto _ : state) {
        auto result = fp12ring.fp12_flat_frobenius_map3(a, ring);
        benchmark::DoNotOptimize(result);
    }
}

template<class Ring>
static void BM_FP12_FinalExpAfterInverse_Impl(benchmark::State& state) {
    Ring ring(field_modulus());
    FP12Ring<Ring> fp12ring(ring);
    volatile uint64_t seed_f = 42;
    volatile uint64_t seed_f_inv = 123;
    auto f = generate_random_fp12<Ring>(ring, seed_f);
    auto f_inverse = generate_random_fp12<Ring>(ring, seed_f_inv);
    benchmark::DoNotOptimize(f);
    benchmark::DoNotOptimize(f_inverse);
    for (auto _ : state) {
        auto result = final_exp_after_inverse(f, f_inverse, fp12ring, ring);
        benchmark::DoNotOptimize(result);
    }
}

template<class Ring>
static void BM_FP12_MulByXy0z00_Impl(benchmark::State& state) {
    Ring ring(field_modulus());
    FP12Ring<Ring> fp12ring(ring);
    volatile uint64_t seed_x = 42;
    auto x = generate_random_fp12<Ring>(ring, seed_x);
    using RingType = Ring;
    std::array<typename RingType::StandardElement, 6> y;
    gmp_randclass rng(gmp_randinit_default);
    rng.seed(123);
    BigInt p = field_modulus();
    for (int i = 0; i < 6; i++) {
        y[i] = ring.from_bigint(BigInt(rng.get_z_range(p.get_mpz())));
    }
    auto c_term = FP2<typename RingType::StandardElement, typename RingType::StandardElement>(y[0], y[1]);
    auto x_term = FP2<typename RingType::StandardElement, typename RingType::StandardElement>(y[2], y[3]);
    auto y_term = FP2<typename RingType::StandardElement, typename RingType::StandardElement>(y[4], y[5]);
    benchmark::DoNotOptimize(x);
    benchmark::DoNotOptimize(y);
    for (auto _ : state) {
        typename FP12Ring<RingType>::StandardElement result;
        fp12ring.flat_mulxy00z0(result, x, c_term, x_term, y_term, ring);
        benchmark::DoNotOptimize(result);
    }
}

template<class Ring>
static void BM_Inverse_Impl(benchmark::State& state) {
    Ring ring(field_modulus());
    BLS381AddChainInversion<Ring> inverter(ring);
    gmp_randclass rng(gmp_randinit_default);
    rng.seed(42);
    BigInt p = field_modulus();
    BigInt random_val(rng.get_z_range(p.get_mpz()));
    if (random_val == BigInt(0)) random_val = BigInt(1);
    auto element = ring.from_bigint(random_val);
    benchmark::DoNotOptimize(element);
    for (auto _ : state) {
        auto result = inverter.invert(element, ring);
        benchmark::DoNotOptimize(result);
    }
}

template<class Ring>
static void BM_MillerLoop_DoubleStep_Impl(benchmark::State& state) {
    Ring ring(field_modulus());
    FP2Ring<Ring> fp2ring(ring);
    FP12Ring<Ring> fp12ring(ring);
    AffinePoint<typename Ring::StandardElement> e1;
    AffinePoint<typename FP2Ring<Ring>::StandardElement> e2;
    setup_g1_g2_generators(ring, e1, e2);
    Miller<Ring> miller(e1, e2, ring);
    ProjectivePoint<typename FP2Ring<Ring>::StandardElement> P{e2.x, e2.y, FP2Ring<Ring>::one(ring)};
    auto f = fp12ring.one(ring);
    benchmark::DoNotOptimize(P);
    benchmark::DoNotOptimize(f);
    for (auto _ : state) {
        P = miller.double_step(P, f, ring, fp2ring, fp12ring);
        benchmark::DoNotOptimize(P);
        benchmark::DoNotOptimize(f);
    }
}

template<class Ring>
static void BM_MillerLoop_AddStep_Impl(benchmark::State& state) {
    Ring ring(field_modulus());
    FP2Ring<Ring> fp2ring(ring);
    FP12Ring<Ring> fp12ring(ring);
    AffinePoint<typename Ring::StandardElement> e1;
    AffinePoint<typename FP2Ring<Ring>::StandardElement> e2;
    setup_g1_g2_generators(ring, e1, e2);
    Miller<Ring> miller(e1, e2, ring);
    ProjectivePoint<typename FP2Ring<Ring>::StandardElement> P{e2.x, e2.y, FP2Ring<Ring>::one(ring)};
    auto f = fp12ring.one(ring);
    P = miller.double_step(P, f, ring, fp2ring, fp12ring);
    benchmark::DoNotOptimize(P);
    benchmark::DoNotOptimize(f);
    for (auto _ : state) {
        P = miller.add_step(P, f, ring, fp2ring, fp12ring);
        benchmark::DoNotOptimize(P);
        benchmark::DoNotOptimize(f);
    }
}

template<class Ring>
static void BM_MillerLoop_Impl(benchmark::State& state) {
    Ring ring(field_modulus());
    FP2Ring<Ring> fp2ring(ring);
    FP12Ring<Ring> fp12ring(ring);
    AffinePoint<typename Ring::StandardElement> e1;
    AffinePoint<typename FP2Ring<Ring>::StandardElement> e2;
    setup_g1_g2_generators(ring, e1, e2);
    Miller<Ring> miller(e1, e2, ring);
    for (auto _ : state) {
        auto result = miller.miller_loop(ring, fp2ring, fp12ring);
        benchmark::DoNotOptimize(result);
    }
}

template<class Ring>
static void BM_FP12_Inverse_BLST_Impl(benchmark::State& state) {
    Ring ring(field_modulus());
    volatile uint64_t seed = 42;
    auto a = generate_random_fp12<Ring>(ring, seed);
    benchmark::DoNotOptimize(a);
    for (auto _ : state) {
        auto inv = invert_fp12_via_blst(ring, a);
        benchmark::DoNotOptimize(inv);
    }
}

template<class Ring>
static void BM_FP12_Inverse_RNS_BLST_Impl(benchmark::State& state) {
    Ring ring(field_modulus());
    FP12Ring<Ring> fp12ring(ring);
    FP12RNS_BLSTInverter<Ring> inverter;
    volatile uint64_t seed = 42;
    auto a = generate_random_fp12<Ring>(ring, seed);
    benchmark::DoNotOptimize(a);
    for (auto _ : state) {
        auto inv = inverter.invert(a, fp12ring, ring);
        benchmark::DoNotOptimize(inv);
    }
}

template<class Ring>
static void BM_FP12_FinalExp_BLST_Inverter_Impl(benchmark::State& state) {
    Ring ring(field_modulus());
    FP12Ring<Ring> fp12ring(ring);
    FP12BLSTInverter<Ring> inverter;
    volatile uint64_t seed_f = 42;
    auto f = generate_random_fp12<Ring>(ring, seed_f);
    benchmark::DoNotOptimize(f);
    for (auto _ : state) {
        auto result = final_exp(f, inverter, fp12ring, ring);
        benchmark::DoNotOptimize(result);
    }
}

template<class Ring>
static void BM_FP12_FinalExp_RNS_BLST_Inverter_Impl(benchmark::State& state) {
    Ring ring(field_modulus());
    FP12Ring<Ring> fp12ring(ring);
    FP12RNS_BLSTInverter<Ring> inverter;
    volatile uint64_t seed_f = 42;
    auto f = generate_random_fp12<Ring>(ring, seed_f);
    benchmark::DoNotOptimize(f);
    for (auto _ : state) {
        auto result = final_exp(f, inverter, fp12ring, ring);
        benchmark::DoNotOptimize(result);
    }
}

template<class Ring>
static void BM_Pairing_BLST_Inverter_Impl(benchmark::State& state) {
    Ring ring(field_modulus());
    FP2Ring<Ring> fp2ring(ring);
    FP12Ring<Ring> fp12ring(ring);
    AffinePoint<typename Ring::StandardElement> e1;
    AffinePoint<typename FP2Ring<Ring>::StandardElement> e2;
    setup_g1_g2_generators(ring, e1, e2);
    FP12BLSTInverter<Ring> inverter;
    for (auto _ : state) {
        auto result = pairing(e1, e2, inverter, fp2ring, fp12ring, ring);
        benchmark::DoNotOptimize(result);
    }
}

template<class Ring>
static void BM_Pairing_RNS_BLST_Inverter_Impl(benchmark::State& state) {
    Ring ring(field_modulus());
    FP2Ring<Ring> fp2ring(ring);
    FP12Ring<Ring> fp12ring(ring);
    AffinePoint<typename Ring::StandardElement> e1;
    AffinePoint<typename FP2Ring<Ring>::StandardElement> e2;
    setup_g1_g2_generators(ring, e1, e2);
    FP12RNS_BLSTInverter<Ring> inverter;
    for (auto _ : state) {
        auto result = pairing(e1, e2, inverter, fp2ring, fp12ring, ring);
        benchmark::DoNotOptimize(result);
    }
}

template<class Ring>
static void BM_G1_Add_Impl(benchmark::State& state) {
    Ring ring(field_modulus());
    G1<Ring> g1;
    volatile uint64_t seed1 = 42;
    volatile uint64_t seed2 = 123;
    auto P = generate_random_g1_projective<Ring>(ring, seed1);
    auto Q = generate_random_g1_projective<Ring>(ring, seed2);
    benchmark::DoNotOptimize(P);
    benchmark::DoNotOptimize(Q);
    for (auto _ : state) {
        auto result = g1.add_point(P, Q, ring);
        benchmark::DoNotOptimize(result);
    }
}

template<class Ring>
static void BM_G1_MixedAdd_Impl(benchmark::State& state) {
    Ring ring(field_modulus());
    G1<Ring> g1;
    volatile uint64_t seed1 = 42;
    volatile uint64_t seed2 = 123;
    auto P = generate_random_g1_projective<Ring>(ring, seed1);
    auto Q = generate_random_g1_affine<Ring>(ring, seed2);
    benchmark::DoNotOptimize(P);
    benchmark::DoNotOptimize(Q);
    for (auto _ : state) {
        auto result = g1.add_mixed_point(P, Q, ring);
        benchmark::DoNotOptimize(result);
    }
}

template<class Ring>
static void BM_G1_Double_Impl(benchmark::State& state) {
    Ring ring(field_modulus());
    G1<Ring> g1;
    volatile uint64_t seed = 42;
    auto P = generate_random_g1_projective<Ring>(ring, seed);
    benchmark::DoNotOptimize(P);
    for (auto _ : state) {
        auto result = g1.double_point(P, ring);
        benchmark::DoNotOptimize(result);
    }
}

template<class Ring>
static void BM_G2_Add_Impl(benchmark::State& state) {
    Ring ring(field_modulus());
    FP2Ring<Ring> fp2ring(ring);
    G2<Ring> g2;
    volatile uint64_t seed1 = 42;
    volatile uint64_t seed2 = 123;
    auto P = generate_random_g2_projective(fp2ring, ring, seed1);
    auto Q = generate_random_g2_projective(fp2ring, ring, seed2);
    benchmark::DoNotOptimize(P);
    benchmark::DoNotOptimize(Q);
    for (auto _ : state) {
        auto result = g2.add_point(P, Q, fp2ring);
        benchmark::DoNotOptimize(result);
    }
}

template<class Ring>
static void BM_G2_MixedAdd_Impl(benchmark::State& state) {
    Ring ring(field_modulus());
    FP2Ring<Ring> fp2ring(ring);
    G2<Ring> g2;
    volatile uint64_t seed1 = 42;
    volatile uint64_t seed2 = 123;
    auto P = generate_random_g2_projective(fp2ring, ring, seed1);
    auto Q = generate_random_g2_affine<Ring>(ring, seed2);
    benchmark::DoNotOptimize(P);
    benchmark::DoNotOptimize(Q);
    for (auto _ : state) {
        auto result = g2.add_mixed_point(P, Q, fp2ring);
        benchmark::DoNotOptimize(result);
    }
}

template<class Ring>
static void BM_G2_Double_Impl(benchmark::State& state) {
    Ring ring(field_modulus());
    FP2Ring<Ring> fp2ring(ring);
    G2<Ring> g2;
    volatile uint64_t seed = 42;
    auto P = generate_random_g2_projective(fp2ring, ring, seed);
    benchmark::DoNotOptimize(P);
    for (auto _ : state) {
        auto result = g2.double_point(P, fp2ring);
        benchmark::DoNotOptimize(result);
    }
}

template<class Ring>
static void BM_G1_ScalarMult_GLV_Impl(benchmark::State& state) {
    Ring ring(field_modulus());
    G1<Ring> g1_curve;
    GLVMult<Ring> glv_mult(g1_curve, ring);
    volatile uint64_t seed_point = 42;
    auto P = generate_random_g1_projective<Ring>(ring, seed_point);
    BigInt r_scalar(scalar_modulus_hex, 16);
    gmp_randclass rng(gmp_randinit_default);
    rng.seed(123);
    BigInt scalar_bigint(rng.get_z_range(r_scalar.get_mpz()));
    std::array<uint8_t, 32> scalar_bytes = {0};
    bigint_to_bytes_le(scalar_bytes.data(), scalar_bigint, 32);
    benchmark::DoNotOptimize(P);
    benchmark::DoNotOptimize(scalar_bytes);
    for (auto _ : state) {
        auto R = glv_mult.multiply(P, scalar_bytes, ring);
        benchmark::DoNotOptimize(R);
    }
}

template<class Ring>
static void BM_G2_ScalarMult_GLS_Impl(benchmark::State& state) {
    Ring ring(field_modulus());
    FP2Ring<Ring> fp2ring(ring);
    GLSMult<Ring> gls_mult(ring);
    typename GLSMult<Ring>::G2Curve g2_curve;
    volatile uint64_t seed_point = 42;
    auto P = generate_random_g2_projective(fp2ring, ring, seed_point);
    BigInt r_scalar(scalar_modulus_hex, 16);
    gmp_randclass rng(gmp_randinit_default);
    rng.seed(321);
    BigInt scalar_bigint(rng.get_z_range(r_scalar.get_mpz()));
    std::array<uint8_t, 32> scalar_bytes = {0};
    bigint_to_bytes_le(scalar_bytes.data(), scalar_bigint, 32);
    benchmark::DoNotOptimize(P);
    benchmark::DoNotOptimize(scalar_bytes);
    for (auto _ : state) {
        auto R = gls_mult.multiply(P, scalar_bytes, fp2ring, g2_curve, ring);
        benchmark::DoNotOptimize(R);
    }
}

} // namespace bench_bls12_381
