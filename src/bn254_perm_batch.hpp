#pragma once

#include "change_base.hpp"
#include "bn254_perm.hpp"

template <class Ring, int batch, int64_t A = Ring::MAX_ADD>
INLINE auto bn254_batch_reduce_perm_bound(
    const Ring &ring,
    const std::array<typename Ring::template ReadyElement<A>, batch> &c_in) {
    static_assert(Ring::TRUE_LIMBS == BN254_TRUE_LIMBS && Ring::LIMBS == BN254_PADDED_LIMBS,
        "bn254 perm batch requires 6×8 ring");
    auto c_m2_acc = batch_change_base_perm<batch, Ring::QUO_EST_ERROR>(c_in, bn254_r1_perm::r1_fixed_perm_wide_reducer);
    using MontHalfType = decltype(ring.montgomery_reduce_m2.reduce_wide(c_m2_acc[0]));
    using OutType = SmallElement<MontHalfType, Bounds<0, 1>>;
    std::array<OutType, batch> result;
    for (int i = 0; i < batch; i++) {
        result[i] = OutType(ring.montgomery_reduce_m2.reduce_wide(c_m2_acc[i]));
    }
    return result;
}

template <class Ring, int batch>
INLINE std::array<typename Ring::StandardElement, batch> bn254_batch_expand_perm(
    const Ring &ring,
    const std::array<typename Ring::SmallStandardElement, batch> &c_m2) {
    static_assert(Ring::TRUE_LIMBS == BN254_TRUE_LIMBS && Ring::LIMBS == BN254_PADDED_LIMBS,
        "bn254 perm batch requires 6×8 ring");
    constexpr int limbs = Ring::StandardBoundedElement::LIMBS;
    constexpr int element_bits = Ring::StandardBoundedElement::ELEMENT_BITS;
    auto wide_zero = WideZero<limbs, element_bits>();
    using ZeroAccElement = ReadyToReduce<decltype(wide_zero), typename Ring::StandardBoundedElement>;
    std::array<ZeroAccElement, batch> zero_acc_elements;
    for (int i = 0; i < batch; i++) {
        zero_acc_elements[i] = ZeroAccElement(wide_zero, c_m2[i].element);
    }
    auto c_m1_acc = batch_change_base_perm<batch>(zero_acc_elements, bn254_r1_perm::r2_fixed_perm_wide_reducer);
    std::array<typename Ring::StandardElement, batch> result;
    for (int i = 0; i < batch; i++) {
        result[i].m1 = ring.montgomery_reduce_m1.template reduce_auto<Ring::STANDARD_BOUND>(c_m1_acc[i]);
        result[i].m2 = c_m2[i].element;
    }
    return result;
}

template <class Ring, int batch>
INLINE std::array<typename Ring::StandardElement, batch> bn254_batch_reduce_expand_perm(
    const Ring &ring,
    const std::array<typename Ring::template ReadyElement<static_cast<int64_t>(Ring::MAX_ADD)>, batch> &elements) {
    static_assert(Ring::TRUE_LIMBS == BN254_TRUE_LIMBS && Ring::LIMBS == BN254_PADDED_LIMBS,
        "bn254 perm batch requires 6×8 ring");
    auto cm2_acc = bn254_batch_reduce_perm_bound<Ring, batch>(ring, elements);
    std::array<typename Ring::SmallStandardElement, batch> c_m2;
    for (int i = 0; i < batch; i++) {
        c_m2[i] = SmallElement<typename Ring::StandardBoundedElement, Bounds<0, 1>>(
            ring.montgomery_reduce_m2.template reduce_auto<Ring::STANDARD_BOUND>(cm2_acc[i].element));
    }
    return bn254_batch_expand_perm<Ring, batch>(ring, c_m2);
}
