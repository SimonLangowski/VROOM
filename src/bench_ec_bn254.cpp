// BN254 G1 + BatchReduceExpand benchmarks (6×46 RNS, IntRNS4).
#include <benchmark/benchmark.h>
#include <array>
#include <cstdint>
#include "bn254_ring_params.hpp"
#include "ec.hpp"
#include "../cpu/precompute/gmp_wrapper.hpp"

using namespace bn254_params;

template<class Ring>
typename G1<Ring>::ProjPoint generate_random_g1_projective(const Ring& ring, uint64_t seed) {
    BigInt p(BN254_P);
    gmp_randclass rng(gmp_randinit_default);
    rng.seed(seed);
    BigInt x(rng.get_z_range(p.get_mpz()));
    BigInt y(rng.get_z_range(p.get_mpz()));
    return typename G1<Ring>::ProjPoint(ring.from_bigint(x), ring.from_bigint(y), ring.one());
}

template<class Ring>
AffinePoint<typename Ring::StandardElement> generate_random_g1_affine(const Ring& ring, uint64_t seed) {
    BigInt p(BN254_P);
    gmp_randclass rng(gmp_randinit_default);
    rng.seed(seed);
    BigInt x(rng.get_z_range(p.get_mpz()));
    BigInt y(rng.get_z_range(p.get_mpz()));
    return AffinePoint<typename Ring::StandardElement>{ring.from_bigint(x), ring.from_bigint(y)};
}

template<int Batch, class Ring>
static void BM_BatchReduceExpand_Impl(benchmark::State& state, const char *label) {
    BigInt p(BN254_P);
    Ring ring(p);

    gmp_randclass rng(gmp_randinit_default);
    rng.seed(42);

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
    state.SetLabel(label);
}

template<int Batch>
static void BM_BatchReduceExpand(benchmark::State& state) {
    BM_BatchReduceExpand_Impl<Batch, Bn254Ring>(state, ("BatchReduceExpand_Matrix_" + std::to_string(Batch)).c_str());
}

template<int Batch>
static void BM_BatchReduceExpandFixedPerm(benchmark::State& state) {
    BM_BatchReduceExpand_Impl<Batch, Bn254RingFixedPerm>(
        state, ("BatchReduceExpand_FixedPerm_" + std::to_string(Batch)).c_str());
}

template<class Ring>
static void BM_G1_Add_Impl(benchmark::State& state, const char *label) {
    BigInt p(BN254_P);
    Ring ring(p);
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
    state.SetLabel(label);
}

template<class Ring>
static void BM_G1_MixedAdd_Impl(benchmark::State& state, const char *label) {
    BigInt p(BN254_P);
    Ring ring(p);
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
    state.SetLabel(label);
}

template<class Ring>
static void BM_G1_Double_Impl(benchmark::State& state, const char *label) {
    BigInt p(BN254_P);
    Ring ring(p);
    G1<Ring> g1;

    volatile uint64_t seed = 42;
    auto P = generate_random_g1_projective<Ring>(ring, seed);
    benchmark::DoNotOptimize(P);

    for (auto _ : state) {
        auto result = g1.double_point(P, ring);
        benchmark::DoNotOptimize(result);
    }
    state.SetLabel(label);
}

template<class Ring>
static void BM_G1_AddSmall_Impl(benchmark::State& state, const char *label) {
    BigInt p(BN254_P);
    Ring ring(p);
    G1<Ring> g1;

    volatile uint64_t seed1 = 42;
    volatile uint64_t seed2 = 123;
    auto P = generate_random_g1_projective<Ring>(ring, seed1);
    auto Q = generate_random_g1_projective<Ring>(ring, seed2);
    benchmark::DoNotOptimize(P);
    benchmark::DoNotOptimize(Q);

    for (auto _ : state) {
        auto result = g1.add_point_small(P, Q, ring);
        benchmark::DoNotOptimize(result);
    }
    state.SetLabel(label);
}

template<class Ring>
static void BM_ModMul_Impl(benchmark::State& state, const char *label) {
    BigInt p(BN254_P);
    Ring ring(p);

    gmp_randclass rng(gmp_randinit_default);
    rng.seed(42);
    auto a = ring.from_bigint(BigInt(rng.get_z_range(p.get_mpz())));
    auto b = ring.from_bigint(BigInt(rng.get_z_range(p.get_mpz())));
    benchmark::DoNotOptimize(a);
    benchmark::DoNotOptimize(b);

    for (auto _ : state) {
        auto result = ring.modmul(a, b);
        benchmark::DoNotOptimize(result);
    }
    state.SetLabel(label);
}

static void BM_ModMul(benchmark::State& state) { BM_ModMul_Impl<Bn254Ring>(state, "ModMul_Matrix"); }
static void BM_ModMul_FixedPerm(benchmark::State& state) { BM_ModMul_Impl<Bn254RingFixedPerm>(state, "ModMul_FixedPerm"); }

static void BM_G1_Add(benchmark::State& state) { BM_G1_Add_Impl<Bn254Ring>(state, "G1_Add_Matrix"); }
static void BM_G1_Add_FixedPerm(benchmark::State& state) { BM_G1_Add_Impl<Bn254RingFixedPerm>(state, "G1_Add_FixedPerm"); }

static void BM_G1_AddSmall(benchmark::State& state) { BM_G1_AddSmall_Impl<Bn254Ring>(state, "G1_AddSmall_Matrix"); }
static void BM_G1_AddSmall_FixedPerm(benchmark::State& state) {
    BM_G1_AddSmall_Impl<Bn254RingFixedPerm>(state, "G1_AddSmall_FixedPerm");
}

static void BM_G1_MixedAdd(benchmark::State& state) { BM_G1_MixedAdd_Impl<Bn254Ring>(state, "G1_MixedAdd_Matrix"); }
static void BM_G1_MixedAdd_FixedPerm(benchmark::State& state) {
    BM_G1_MixedAdd_Impl<Bn254RingFixedPerm>(state, "G1_MixedAdd_FixedPerm");
}

static void BM_G1_Double(benchmark::State& state) { BM_G1_Double_Impl<Bn254Ring>(state, "G1_Double_Matrix"); }
static void BM_G1_Double_FixedPerm(benchmark::State& state) {
    BM_G1_Double_Impl<Bn254RingFixedPerm>(state, "G1_Double_FixedPerm");
}

BENCHMARK(BM_ModMul);
BENCHMARK(BM_ModMul_FixedPerm);
BENCHMARK(BM_BatchReduceExpand<1>);
BENCHMARK(BM_BatchReduceExpand<2>);
BENCHMARK(BM_BatchReduceExpand<3>);
BENCHMARK(BM_BatchReduceExpand<4>);
BENCHMARK(BM_BatchReduceExpand<5>);
BENCHMARK(BM_BatchReduceExpand<6>);
BENCHMARK(BM_BatchReduceExpandFixedPerm<1>);
BENCHMARK(BM_BatchReduceExpandFixedPerm<2>);
BENCHMARK(BM_BatchReduceExpandFixedPerm<3>);
BENCHMARK(BM_BatchReduceExpandFixedPerm<4>);
BENCHMARK(BM_BatchReduceExpandFixedPerm<5>);
BENCHMARK(BM_BatchReduceExpandFixedPerm<6>);
BENCHMARK(BM_G1_Add);
BENCHMARK(BM_G1_Add_FixedPerm);
BENCHMARK(BM_G1_AddSmall);
BENCHMARK(BM_G1_AddSmall_FixedPerm);
BENCHMARK(BM_G1_MixedAdd);
BENCHMARK(BM_G1_MixedAdd_FixedPerm);
BENCHMARK(BM_G1_Double);
BENCHMARK(BM_G1_Double_FixedPerm);

BENCHMARK_MAIN();
