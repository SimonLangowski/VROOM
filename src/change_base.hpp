#pragma once
#include "../cpu/vector/vector_impl.hpp"
#include "constexpr_utils.hpp"
#include "ring_element.hpp"
#include "bounded_element.hpp"
#include "bounds.hpp"
#include <cstdint>
#include <cstddef>
#include <array>
#include "preprocess.hpp"

enum class ChangeBaseVariation {
    Matrix,
    FixedPerm,
    MatrixNoK
};

/** Logical RNS limb count when storage is AVX-padded (bn254: 8 stored, 6 active). */
template<class ReducedM1>
constexpr int change_base_true_limbs() {
    if constexpr (ReducedM1::LIMBS == 8 && ReducedM1::ELEMENT_BITS == 46) {
        return 6;
    }
    return ReducedM1::LIMBS;
}

/*
template<int limbs>
class RNSMatrix2 {
    constexpr static int PADDED_LIMBS = ceil_div(limbs, AVXVector<limbs>::LIMBS_PER_VEC) * AVXVector<limbs>::LIMBS_PER_VEC;
    alingas(64) const std::array<std::array<uint64_t, PADDED_LIMBS>, limbs> rns_mat;
    alingas(64) const std::array<uint64_t, PADDED_LIMBS> shifted_quotient_estimations;
    alingas(64) const std::array<uint64_t, PADDED_LIMBS> correction;
    alingas(64) const std::array<uint64_t, PADDED_LIMBS> correction_shift;

    template<int batch>
    std::pair<std::array<AVXVector<limbs>, batch>, std::array<AVXVector<limbs>, batch>> rns_reduce_raw_batch(const std::array<AVXVector<limbs>, batch> &residues, const AVXVector<limbs> &hi_acc_c, const AVXVector<limbs> &lo_acc_c) const {
        std::array<AVXVector<limbs>, batch> out_hi;
        std::array<AVXVector<limbs>, batch> out_lo;
        std::array<uint128_t, batch> k_raw;
        alignas(64) std::array<std::array<uint64_t, PADDED_LIMBS>, batch> scalars;
        #pragma clang loop unroll (enable)
        for (int i = 0; i < batch; i++) {
            residues[i].store(scalars[i].data()); // No-op if residues passed by pointer.
            k_raw[i] = 0;
        }
        AVXVector<limbs> hi_acc = hi_acc_c; // Ideally, load from pointer?
        AVXVector<limbs> lo_acc = lo_acc_c; 
        #pragma clang loop unroll (enable)
        for (int i = 0; i < limbs; i++) {
            AVXVector<limbs> rns_mat_i = AVXVector<limbs>(rns_mat[i]);
            std::array<uint64_t, batch> scalars_l;
            std::array<AVXScalar, batch> scalar_v;
            #pragma clang loop unroll (enable)
            for (int j = 0; j < batch; j++) {
                scalars_l[j] = scalars[j][i];
                scalar_v[j] = Scalar(scalars_l[j]);
            }
            // Do integer operations (might hide vector load latency)
            #pragma clang loop unroll (enable)
            for (int j = 0; j < batch; j++) {
                k_raw[j] += (uint128_t)(scalars_l[j]) * (uint128_t)(shifted_quotient_estimations[i]);
            }
            #pragma clang loop unroll (enable)
            for (int j = 0; j < batch; j++) {
                if (i == 0) {
                    out_hi[j] = hi_acc.mulhi_scalar(rns_mat_i, scalar_v[j]);
                    out_lo[j] = lo_acc.mullo_scalar(rns_mat_i, scalar_v[j]);
                } else {
                    out_hi[j] = out_hi[j].mulhi_scalar(rns_mat_i, scalar_v[j]);
                    out_lo[j] = out_lo[j].mullo_scalar(rns_mat_i, scalar_v[j]);
                }
            }
        }
        AVXVector<limbs> correction_v = AVXVector<limbs>(correction);
        AVXVector<limbs> correction_shift_v = AVXVector<limbs>(correction_shift);
        std::array<AVXScalar, batch> k_shifted;
        #pragma clang loop unroll (enable)
        for (int j = 0; j < batch; j++) {
            k_shifted[j] = Scalar((uint64_t)(k_raw[j] >> precision));
        }
        #pragma clang loop unroll (enable)
        for (int j = 0; j < batch; j++) {
            out_hi[j] = out_hi[j].mulhi_scalar(correction_v, k_shifted[j]);
            out_lo[j] = out_lo[j].mullo_scalar(correction_v, k_shifted[j]);
            k_shifted[j] = scalar_shift(k_shifted[j], 52);
            // k = k_hi * 2^52 + k_lo
            // correction_shift_v = (k_hi * 2^52) % modulus
            out_lo[j] = out_lo[j].mullo_scalar(correction_shift_v, k_shifted[j]);
            out_hi[j] = out_hi[j].mulhi_scalar(correction_shift_v, k_shifted[j]);
        }
        // Should be RVO in C++17
        // Since C++ ABI is that all ZMM registers are clobbered, there should be no memory overhead...
        return {out_hi, out_lo};
    }
}
 */

template<class WideM2, class ReducedM1>
class ReadyToReduce {
    public:
    WideM2 acc;
    ReducedM1 input;
    
    // Constructor for aggregate initialization compatibility
    ReadyToReduce() = default;
    ReadyToReduce(const WideM2 &a, const ReducedM1 &i) : acc(a), input(i) {}
};

// TODO: add chunk size parameter and chunk accordingly?
template<std::size_t batch, class RNSMatrix, class WideM2, class ReducedM1>
INLINE auto batch_change_base(const std::array<ReadyToReduce<WideM2, ReducedM1>, batch> &inputs, const RNSMatrix &rns_matrix) {
    constexpr int limbs = ReducedM1::LIMBS;
    constexpr int element_bits = ReducedM1::ELEMENT_BITS;
    ReducedM1::check_mult_bounds();
    std::array<AVXVector<limbs>, batch> out_hi;
    std::array<AVXVector<limbs>, batch> out_lo;
    std::array<AVXVector<limbs>, batch> residues1;
    for (size_t i = 0; i < batch; i++) {
        residues1[i] = inputs[i].input.data;
        out_hi[i] = inputs[i].acc.high.data;
        out_lo[i] = inputs[i].acc.low.data;
    }
    rns_matrix.template rns_reduce_raw_batch<batch, 2>(&residues1, out_hi, out_lo);
    auto change_base_bounds = Bounds<0, 2*change_base_true_limbs<ReducedM1>()*ReducedM1::bounds.upper>{};
    auto total_bounds = inputs[0].acc.bounds + change_base_bounds;
    std::array<WideElement<decltype(total_bounds), limbs, element_bits>, batch> result;
    for (size_t i = 0; i < batch; i++) {
        result[i].low.data = out_lo[i];
        result[i].high.data = out_hi[i];
    }
    return result;
}

template<class PermReducer>
constexpr bool perm_reducer_skips_k() {
    if constexpr (requires { PermReducer::SKIP_K_CORRECTION; }) {
        return PermReducer::SKIP_K_CORRECTION;
    }
    return false;
}

/** Same as batch_change_base but reducer supplies rns_reduce_raw_batch (e.g. bn254 dot-first perm). */
template<std::size_t batch, int QuoEstError = 0, class PermReducer, class WideM2, class ReducedM1>
INLINE auto batch_change_base_perm(
    const std::array<ReadyToReduce<WideM2, ReducedM1>, batch> &inputs,
    const PermReducer &reducer) {
    constexpr int limbs = ReducedM1::LIMBS;
    constexpr int true_limbs = change_base_true_limbs<ReducedM1>();
    constexpr int element_bits = ReducedM1::ELEMENT_BITS;
    ReducedM1::check_mult_bounds();
    if constexpr (perm_reducer_skips_k<PermReducer>()) {
        static_assert(QuoEstError > 0, "skip-k perm reducer requires QuoEstError template argument");
        static_assert(
            ReducedM1::bounds.upper * true_limbs <= QuoEstError,
            "skip-k: m1 upper * true_limbs must be <= QUO_EST_ERROR; increase QUO_EST_ERROR or reduce m1 to STANDARD_BOUND (2)");
    }
    std::array<AVXVector<limbs>, batch> out_hi;
    std::array<AVXVector<limbs>, batch> out_lo;
    std::array<AVXVector<limbs>, batch> residues1;
    for (size_t i = 0; i < batch; i++) {
        residues1[i] = inputs[i].input.data;
        out_hi[i] = inputs[i].acc.high.data;
        out_lo[i] = inputs[i].acc.low.data;
    }
    reducer.template rns_reduce_raw_batch<batch>(&residues1, out_hi, out_lo);
    auto change_base_bounds = [&]() {
        if constexpr (perm_reducer_skips_k<PermReducer>()) {
            return Bounds<0, true_limbs * ReducedM1::bounds.upper>{};
        } else {
            return Bounds<0, 2 * true_limbs * ReducedM1::bounds.upper>{};
        }
    }();
    auto total_bounds = inputs[0].acc.bounds + change_base_bounds;
    std::array<WideElement<decltype(total_bounds), limbs, element_bits>, batch> result;
    for (size_t i = 0; i < batch; i++) {
        result[i].low.data = out_lo[i];
        result[i].high.data = out_hi[i];
    }
    return result;
}

template<std::size_t batch, int QuoEstError, class PermReducer, class WideM2, class ReducedM1>
INLINE auto batch_change_base_perm_qr(
    const std::array<ReadyToReduce<WideM2, ReducedM1>, batch> &inputs,
    const PermReducer &reducer) {
    constexpr int limbs = ReducedM1::LIMBS;
    constexpr int true_limbs = change_base_true_limbs<ReducedM1>();
    constexpr int element_bits = ReducedM1::ELEMENT_BITS;
    ReducedM1::check_mult_bounds();
    static_assert(perm_reducer_skips_k<PermReducer>(), "batch_change_base_perm_qr requires a skip-k perm reducer");
    static_assert(
        ReducedM1::bounds.upper * true_limbs <= QuoEstError,
        "skip-k: m1 upper * true_limbs must be <= QUO_EST_ERROR; increase QUO_EST_ERROR or reduce m1 to STANDARD_BOUND (2)");
    std::array<AVXVector<limbs>, batch> out_hi;
    std::array<AVXVector<limbs>, batch> out_lo;
    std::array<AVXVector<limbs>, batch> residues1;
    for (size_t i = 0; i < batch; i++) {
        residues1[i] = inputs[i].input.data;
        out_hi[i] = inputs[i].acc.high.data;
        out_lo[i] = inputs[i].acc.low.data;
    }
    reducer.template rns_reduce_raw_batch<batch>(&residues1, out_hi, out_lo);
    // No k correction applied, so no increase in bounds
    // Possible Montgomery reduction could fully reduce to bound 2?
    auto change_base_bounds = Bounds<0, true_limbs*ReducedM1::bounds.upper>{};
    auto total_bounds = inputs[0].acc.bounds + change_base_bounds;
    std::array<WideElement<decltype(total_bounds), limbs, element_bits>, batch> result;
    for (size_t i = 0; i < batch; i++) {
        result[i].low.data = out_lo[i];
        result[i].high.data = out_hi[i];
    }
    // No k correction applied, so RNS bounds are larger
    // expanded_element with <0, limbs*ReducedM1::bounds.upper> bounds
    // better be less than N/2
    std::array<ExpandedElement<WideElement<decltype(total_bounds), limbs, element_bits>, Bounds<0, true_limbs*ReducedM1::bounds.upper>>, batch> result_expanded;
    for (size_t i = 0; i < batch; i++) {
        result_expanded[i] = ExpandedElement<WideElement<decltype(total_bounds), limbs, element_bits>, Bounds<0, true_limbs*ReducedM1::bounds.upper>>(result[i], result[i]);
    }
    return result_expanded;
}