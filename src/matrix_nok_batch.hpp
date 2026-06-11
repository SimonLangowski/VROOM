#pragma once

#include "change_base.hpp"
#include "bls12_381_perm.hpp"
#include "bls12_381_perm_params.hpp"

namespace matrix_nok {

/** r1 perm (no k) → m2 reduce_wide; shared by MatrixNoK rings. */
template <class Ring, int batch, class R1PermReducer, int64_t A = Ring::MAX_ADD>
INLINE auto batch_reduce_bound(
    const Ring &ring,
    const std::array<typename Ring::template ReadyElement<A>, batch> &c_in,
    const R1PermReducer &r1_reducer) {
    auto c_m2_acc = batch_change_base_perm<batch, Ring::QUO_EST_ERROR>(c_in, r1_reducer);
    using MontHalfType = decltype(ring.montgomery_reduce_m2.reduce_wide(c_m2_acc[0]));
    using OutType = SmallElement<MontHalfType, Bounds<0, 1>>;
    std::array<OutType, batch> result;
    for (int i = 0; i < batch; i++) {
        result[i] = OutType(ring.montgomery_reduce_m2.reduce_wide(c_m2_acc[i]));
    }
    return result;
}

/** r1 perm (no k) + scalar r2 matrix expand to StandardElement. */
template <class Ring, int batch, class R1PermReducer>
INLINE std::array<typename Ring::StandardElement, batch> batch_reduce_expand(
    const Ring &ring,
    const std::array<typename Ring::template ReadyElement<static_cast<int64_t>(Ring::MAX_ADD)>, batch> &elements,
    const R1PermReducer &r1_reducer) {
    auto cm2_acc = batch_reduce_bound<Ring, batch, R1PermReducer>(ring, elements, r1_reducer);
    std::array<typename Ring::SmallStandardElement, batch> c_m2;
    for (int i = 0; i < batch; i++) {
        c_m2[i] = SmallElement<typename Ring::StandardBoundedElement, Bounds<0, 1>>(
            ring.montgomery_reduce_m2.template reduce_auto<Ring::STANDARD_BOUND>(cm2_acc[i].element));
    }
    constexpr int limbs = Ring::StandardBoundedElement::LIMBS;
    constexpr int element_bits = Ring::StandardBoundedElement::ELEMENT_BITS;
    auto wide_zero = WideZero<limbs, element_bits>();
    using ZeroAccElement = ReadyToReduce<decltype(wide_zero), typename Ring::StandardBoundedElement>;
    std::array<ZeroAccElement, batch> zero_acc_elements;
    for (int i = 0; i < batch; i++) {
        zero_acc_elements[i] = ZeroAccElement(wide_zero, c_m2[i].element);
    }
    auto c_m1_acc = batch_change_base<batch>(zero_acc_elements, ring.sys.intrns2.r2);
    std::array<typename Ring::StandardElement, batch> result;
    for (int i = 0; i < batch; i++) {
        result[i].m1 = ring.montgomery_reduce_m1.template reduce_auto<Ring::STANDARD_BOUND>(c_m1_acc[i]);
        result[i].m2 = c_m2[i].element;
    }
    return result;
}

template <class Ring, int batch, int64_t A = Ring::MAX_ADD>
INLINE auto dispatch_batch_reduce_bound(
    const Ring &ring,
    const std::array<typename Ring::template ReadyElement<A>, batch> &c_in) {
    static_assert(Ring::TRUE_LIMBS == BLS12_381_TRUE_LIMBS && Ring::LIMBS == BLS12_381_PADDED_LIMBS,
        "MatrixNoK requires BLS12-381 8×8 ring");
    return batch_reduce_bound<Ring, batch>(ring, c_in, bls12_381_r1_perm::r1_fixed_perm_wide_reducer);
}

template <class Ring, int batch>
INLINE std::array<typename Ring::StandardElement, batch> dispatch_batch_reduce_expand(
    const Ring &ring,
    const std::array<typename Ring::template ReadyElement<static_cast<int64_t>(Ring::MAX_ADD)>, batch> &elements) {
    static_assert(Ring::TRUE_LIMBS == BLS12_381_TRUE_LIMBS && Ring::LIMBS == BLS12_381_PADDED_LIMBS,
        "MatrixNoK requires BLS12-381 8×8 ring");
    return batch_reduce_expand<Ring, batch>(ring, elements, bls12_381_r1_perm::r1_fixed_perm_wide_reducer);
}

} // namespace matrix_nok
