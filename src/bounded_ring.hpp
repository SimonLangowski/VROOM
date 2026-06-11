#pragma once
#include "ring_element.hpp"
#include "elementwise_reduce.hpp"
#include "change_base.hpp"
#include "../cpu/precompute/gmp_wrapper.hpp"
#include "../cpu/vector/multiplication.hpp"
#include <array>
#include <tuple>
#include <utility>
#include <cstddef>
#include <type_traits>
#include <string>
#include <iostream>
#include "preprocess.hpp"

#include "bn254_perm_batch.hpp"
#include "bn254_intrns4_system.hpp"
#include "matrix_nok_batch.hpp"
#include "bls12_381_intrns4_system.hpp"
#include "bls12_381_perm_params.hpp"
#include "bn254_perm_params.hpp"

constexpr int constexpr_ceil_log2(uint64_t x) {
    if (x == 0) return 0;
    if (x == 1) return 0;
    return 64 - __builtin_clzll(x - 1);
}

// The maximum negative and positive values that are needed.  These occur in the Y coordinate of doubling_step.
// true_limbs: logical RNS width (RNS_BITS capacity). padded_limbs: AVX storage width (often 8 for AVX-512).
template<int modulus_bits=381, int true_limbs=8, int element_bits=52, int RNS_MAX_NEG=-1932, int RNS_MAX_POS=2377, int LOG_MULTIPLES=4, int MUL_3B_COEFF=12, int MAX_ADD_VALUE=800, int RNS_HALF_MAX_NEG=48, int BASE_REDUNDANT_MULTIPLE=3, int padded_limbs=true_limbs, ChangeBaseVariation change_base=ChangeBaseVariation::Matrix>
class BoundedRing {

    // static constexpr int BASE_REDUNDANT_MULTIPLE = 3;
    static constexpr int WIDE_REDUNDANT_MULTIPLE = BASE_REDUNDANT_MULTIPLE * BASE_REDUNDANT_MULTIPLE;

    // -1 to account for each moduli being slighly less than element_bits
    static constexpr int RNS_BITS = element_bits * true_limbs - 1;

    // We store values only in the positive range, so we need a range of [0, MAX_POS - MAX_NEG]
    static constexpr int RNS_MAX_MULTIPLE = (RNS_MAX_POS - RNS_MAX_NEG);
    static_assert(RNS_MAX_NEG <= 0, "RNS_MAX_NEG must be negative");
    static_assert(RNS_MAX_POS > 0, "RNS_MAX_POS must be positive");
    // PointAdd: mul_3b(y3_18) scales RNS tags by MUL_3B_COEFF; y3_18 can be -2 → need 2× coeff offset.
    static_assert(RNS_HALF_MAX_NEG >= 2 * MUL_3B_COEFF, "RNS_HALF_MAX_NEG must cover mul_3b(y3_18) RNS shift");
    
    
    static constexpr int MODULUS_BITS = modulus_bits;

    static constexpr int BASE_BITS = ceil_log2(BASE_REDUNDANT_MULTIPLE);

    static constexpr int COMBINED_MAX = RNS_MAX_MULTIPLE * WIDE_REDUNDANT_MULTIPLE;
    static constexpr int COMBINED_MAX_BITS = ceil_log2(COMBINED_MAX);

    static_assert(COMBINED_MAX_BITS + MODULUS_BITS + MODULUS_BITS <= RNS_BITS + RNS_BITS, "Numbers of RNS max * modulus^2 must fit in RNS of M N");

    // Bits used for computation in mod N without M part - could add verification on small elements for smaller rns bound potentially.
    static constexpr int BITS_FOR_SMALL_COMPUTATION = ceil_log2(RNS_HALF_MAX_NEG);

    static_assert(BASE_BITS + MODULUS_BITS + BITS_FOR_SMALL_COMPUTATION + 1 <= RNS_BITS, "Numbers after reduction must fit in RNS N/2");

    static_assert(BASE_REDUNDANT_MULTIPLE > 0, "Base redundance must be positive");

    // combined max * p^2 / M + 2p < base_redundance * p
    // combined max * p / M + 2 < base_redundance
    // To avoid big integers, we compute combined max / (M / p) where (M / p) is estimated as 1 << (RNS_BITS - modulus_bits)

    BigInt modulus;
    
    public:

    static constexpr int QUO_EST_ERROR =
        (change_base == ChangeBaseVariation::FixedPerm || change_base == ChangeBaseVariation::MatrixNoK)
            ? 2 * true_limbs : 1;

    static_assert(
        (change_base == ChangeBaseVariation::MatrixNoK
            && BASE_REDUNDANT_MULTIPLE >= QUO_EST_ERROR + 4)
        || (change_base != ChangeBaseVariation::MatrixNoK
            && ceil_div(int64_t{COMBINED_MAX}, int64_t{1} << (RNS_BITS - modulus_bits)) + 1 + QUO_EST_ERROR
                <= BASE_REDUNDANT_MULTIPLE),
        "RNS Montgomery reduction must produce base redundance product for whole range");

    static constexpr int TRUE_LIMBS = true_limbs;
    static constexpr int LIMBS = padded_limbs;
    
    // In theory we could pass these reducers to the IntRNS2System constructor (also would fix the hardcoded 52 bit operations for the IntRNS2System)
    static constexpr int CONVERT_TO_WORD_BITS =
        (change_base == ChangeBaseVariation::MatrixNoK && element_bits == 50) ? 50 : 52;
    using SysType = IntRNS2System<modulus_bits, LIMBS, MontgomeryReduce, CONVERT_TO_WORD_BITS>;
    const SysType sys;
    const MontgomeryReduce2<LIMBS, element_bits, LOG_MULTIPLES> montgomery_reduce_m1;
    const MontgomeryReduce2<LIMBS, element_bits, LOG_MULTIPLES> montgomery_reduce_m2;

    // A 52-bit or less element, which can be used directly in 52-bit multiplication operations.
    constexpr static uint64_t MULT_OK_BOUND = 1ull << (52 - element_bits);
    // The output of the change base operation; if possible, we skip the reduction and output with a redundancy of 2
    constexpr static uint64_t STANDARD_BOUND = (MULT_OK_BOUND < 2ull) ? MULT_OK_BOUND : 2ull;
    static constexpr int NumLimbs = LIMBS;
    // To use array types, we recast everything to some maximum number of additions.  You can probably make this bigger if needed.
    static constexpr uint64_t MAX_ADD = MAX_ADD_VALUE;
    using FullyReducedElement = BoundedElement<Bounds<0, 1>, LIMBS, element_bits>;
    using StandardBoundedElement = BoundedElement<Bounds<0, STANDARD_BOUND>, LIMBS, element_bits>;
    using MultOkBoundedElement = BoundedElement<Bounds<0, MULT_OK_BOUND>, LIMBS, element_bits>;
    using ReadyInputElement = std::conditional_t<
        change_base == ChangeBaseVariation::FixedPerm || change_base == ChangeBaseVariation::MatrixNoK,
        StandardBoundedElement, MultOkBoundedElement>;
    // A fully reduced element, the output of RNS reduction.
    using StandardElement = ExpandedElement<StandardBoundedElement, Bounds<0, 1>>;
    template<int64_t A>
    using ReadyElement = ReadyToReduce<WideElement<Bounds<0, A>, LIMBS, element_bits>, ReadyInputElement>;

    template<int64_t rns_bound>
    using OffsetElement = ExpandedElement<FullyReducedElement, Bounds<rns_bound, rns_bound>>;

    template<int64_t ele_bound, int64_t rns_bound>
    using WideOffsetElement = ExpandedElement<WideElement<Bounds<ele_bound, ele_bound+1>, LIMBS, element_bits>, Bounds<rns_bound, rns_bound>>;

    template<int64_t upper_bound>
    using ReducedElement = ExpandedElement<StandardBoundedElement, Bounds<0, upper_bound>>;

    template<int64_t upper_bound>
    using SmallReducedElement = SmallElement<StandardBoundedElement, Bounds<0, upper_bound>>;

    using SmallStandardElement = SmallReducedElement<1>;

    using RawReduceOutput = SmallElement<Bounds<0, MAX_ADD + 2*true_limbs*MULT_OK_BOUND>, Bounds<0, 1>>;

    // k×3p in RNS (Montgomery m2); prep_expand adds this when RNS metadata is negative.
    template<int64_t k>
    using SmallRNSOffsetAt = SmallElement<StandardBoundedElement, Bounds<k, k>>;

    private:

    // One in RNS
    StandardElement one_;
    // Used to make the largest possible negative value positive.
    // WideOffsetElement<0, -RNS_MAX_NEG> max_neg_offset;
    OffsetElement<-RNS_MAX_NEG> max_neg_offset;
    // the modulus in RNS
    OffsetElement<1> modulus_in_rns;

    OffsetElement<6> six_modulus_in_rns;

    // Precomputed k×3p RNS offsets (GMP to_mont_exact at ring init only).
    // k=4 covers small negative tags (e.g. karatsuba / pairwise diffs with lower ≥ -4).
    // k=42 (= RNS_HALF_MAX_NEG) covers mul_3b-scaled tags (e.g. lower = -21, -42).
    SmallRNSOffsetAt<4> small_rns_offset_4;
    SmallRNSOffsetAt<RNS_HALF_MAX_NEG> sub_rns_offset;

    public:

    BoundedRing(const BigInt &modulus)
        : modulus(modulus)
        // 2 * RNS_MAX to allow for negative range
        // 9p * max redundance for elementwise reduction
        , sys([&]() {
            if constexpr (change_base == ChangeBaseVariation::FixedPerm) {
                static_assert(MUL_3B_COEFF == 9, "FixedPerm requires bn254 (3b=9)");
                return make_BN254_IntRNS4SystemFromExport<modulus_bits, LIMBS>(modulus);
            } else if constexpr (change_base == ChangeBaseVariation::MatrixNoK) {
                static_assert(true_limbs == BLS12_381_TRUE_LIMBS && padded_limbs == BLS12_381_PADDED_LIMBS,
                    "MatrixNoK requires BLS12-381 8×8 ring");
                return make_BLS12_381_IntRNS4SystemFromExport<modulus_bits, LIMBS>(modulus);
            } else {
                return make_IntRNS2SystemMontgomery<modulus_bits, LIMBS>(modulus, element_bits, COMBINED_MAX);
            }
        }())
         // There reducers differ in that they return without conditional correction (this should help when element_bits = 50, but has no advantage at 52)
        , montgomery_reduce_m1(sys.moduli1)
        , montgomery_reduce_m2(sys.moduli2)
        // To and from operate on 52 bit chunks independently of element bits, we let the system handle the digit calculation.
    {
        auto one_exact = sys.to_mont_exact(BigInt(1));
        auto one_m1 = StandardBoundedElement(AVXVector<LIMBS>(one_exact.first));
        auto one_m2 = StandardBoundedElement(AVXVector<LIMBS>(one_exact.second));
        one_ = StandardElement(one_m1, one_m2);
        // This also needs the multiplication by (3p)*(3p) = 9p^2 redundance for products.
        auto max_neg_offset_exact = sys.to_mont_exact(BigInt(-RNS_MAX_NEG*WIDE_REDUNDANT_MULTIPLE*modulus*modulus), true);
        auto max_neg_m1 = FullyReducedElement(AVXVector<LIMBS>(max_neg_offset_exact.first));
        auto max_neg_m2 = FullyReducedElement(AVXVector<LIMBS>(max_neg_offset_exact.second));
        max_neg_offset = OffsetElement<-RNS_MAX_NEG>(max_neg_m1, max_neg_m2);

        // 3p = 1 RNS bound.
        auto modulus_in_rns_exact = sys.to_mont_exact(BASE_REDUNDANT_MULTIPLE*modulus);
        auto modulus_in_rns_m1 = FullyReducedElement(AVXVector<LIMBS>(modulus_in_rns_exact.first));
        auto modulus_in_rns_m2 = FullyReducedElement(AVXVector<LIMBS>(modulus_in_rns_exact.second));
        modulus_in_rns = OffsetElement<1>(modulus_in_rns_m1, modulus_in_rns_m2);
        auto six_modulus_in_rns_exact = sys.to_mont_exact(6*BASE_REDUNDANT_MULTIPLE*modulus);
        auto six_modulus_in_rns_m1 = FullyReducedElement(AVXVector<LIMBS>(six_modulus_in_rns_exact.first));
        auto six_modulus_in_rns_m2 = FullyReducedElement(AVXVector<LIMBS>(six_modulus_in_rns_exact.second));
        six_modulus_in_rns = OffsetElement<6>(six_modulus_in_rns_m1, six_modulus_in_rns_m2);

        small_rns_offset_4 = compute_small_rns_offset<4>();
        sub_rns_offset = compute_small_rns_offset<RNS_HALF_MAX_NEG>();
    }

    static INLINE StandardElement zero() {
        return StandardElement{};
    }

    INLINE StandardElement one() const {
        return one_;
    }
    
    static INLINE auto wrap(const std::pair<AVXVector<LIMBS>, AVXVector<LIMBS>> &elements) {
        return StandardElement(
            StandardBoundedElement(elements.first), 
            StandardBoundedElement(elements.second)
        );
    };

    static INLINE auto wrap(const std::pair<std::array<uint64_t, LIMBS>, std::array<uint64_t, LIMBS>> &elements) {
        return wrap(std::make_pair(AVXVector<LIMBS>(elements.first), AVXVector<LIMBS>(elements.second)));
    }

    template<int bytes>
    INLINE auto from_hex(const char(&lit)[bytes]) const {
        // Strip "0x" or "0X" prefix if present, as GMP's mpz_set_str doesn't accept it when base is explicitly set
        std::string hex_str(lit);
        if (hex_str.length() >= 2 && (hex_str.substr(0, 2) == "0x" || hex_str.substr(0, 2) == "0X")) {
            hex_str = hex_str.substr(2);
        }
        BigInt a = BigInt(hex_str.c_str(), 16);
        return wrap(sys.to_mont_exact(a, false));
    }

    INLINE BigInt to_bigint(const StandardElement &element) const {
        return sys.from_mont_avx(element.m2.data);
    }

    template<int lower_bound, int upper_bound, class RNSBounds>
    INLINE auto recast(const ExpandedElement<BoundedElement<Bounds<lower_bound, upper_bound>, LIMBS, element_bits>, RNSBounds> &element) const {
        auto m1_recast = element.m1.template recast<lower_bound, upper_bound>();
        auto m2_recast = element.m2.template recast<lower_bound, upper_bound>();
        return ExpandedElement<decltype(m1_recast), RNSBounds>(m1_recast, m2_recast);
    }

    template<int multiple=1>
    INLINE auto elementwise_moduli() const {
        static_assert(multiple >= 1, "Multiple must be at least 1");
        // We only store log2(multiple) for efficiency.
        constexpr int required_multiple = ceil_log2(multiple);
        auto m1_moduli = montgomery_reduce_m1.template moduli<required_multiple>();
        auto m2_moduli = montgomery_reduce_m2.template moduli<required_multiple>();
        // This is 0 in the RNS system, so we can use the Bounds<0, 0> type.
        return ExpandedElement<decltype(m1_moduli), Bounds<0, 0>>(m1_moduli, m2_moduli);
    }

    template<int64_t lower, int64_t upper, class RNSBoundsType>
    INLINE auto wideify(const ExpandedElement<BoundedElement<Bounds<lower, upper>, LIMBS, element_bits>, RNSBoundsType> &element) const {
        auto recast_zero = BoundedElement<Bounds<0, 0>, LIMBS, element_bits>::zero().template recast<lower, upper>();
        using WideElementType = WideElement<Bounds<lower, upper>, LIMBS, element_bits>;
        auto m1 = WideElementType(recast_zero, element.m1);
        auto m2 = WideElementType(recast_zero, element.m2);
        return ExpandedElement<WideElementType, RNSBoundsType>(m1, m2);
    }
    
    INLINE auto negative_range_offset() const {
        return wideify(max_neg_offset);
    }

    // Prepare an element for a multiplication - elementwise reduce to make 52 bits.
    template<int64_t lower, int64_t upper, class RNSBoundsType, int to=MULT_OK_BOUND>
    INLINE auto prep(const ExpandedElement<BoundedElement<Bounds<lower, upper>, LIMBS, element_bits>, RNSBoundsType> &element) const {
        auto m1_reduced = montgomery_reduce_m1.template reduce_auto<to, upper, lower>(element.m1);
        auto m2_reduced = montgomery_reduce_m2.template reduce_auto<to, upper, lower>(element.m2);
        return ExpandedElement<decltype(m1_reduced), RNSBoundsType>(m1_reduced, m2_reduced);
    }

    template<class ElementType>
    INLINE auto prep_left(const ElementType &element) const {
        return prep(element);
    }

    template<int64_t lower, int64_t upper, class RNSBoundsType, int to=MULT_OK_BOUND>
    INLINE auto prep(const SmallElement<BoundedElement<Bounds<lower, upper>, LIMBS, element_bits>, RNSBoundsType> &element) const {
        auto m2_reduced = montgomery_reduce_m2.template reduce_auto<to, upper, lower>(element.element);
        return SmallElement<decltype(m2_reduced), RNSBoundsType>(m2_reduced);
    }

    template<int64_t lower, int64_t upper, class RNSBoundsType>
    INLINE auto prep_standard(const ExpandedElement<BoundedElement<Bounds<lower, upper>, LIMBS, element_bits>, RNSBoundsType> &element) const {
        return prep<lower, upper, RNSBoundsType, STANDARD_BOUND>(element);
    }

    template<int64_t lower, int64_t upper, class RNSBoundsType>
    INLINE auto prep_standard(const SmallElement<BoundedElement<Bounds<lower, upper>, LIMBS, element_bits>, RNSBoundsType> &element) const {
        return prep<lower, upper, RNSBoundsType, STANDARD_BOUND>(element);
    }
    
    template<int64_t upper_bound, int64_t lower_bound, class RNSBoundsType>
    INLINE auto negate(const ExpandedElement<BoundedElement<Bounds<lower_bound, upper_bound>, LIMBS, element_bits>, RNSBoundsType> &a) const {
        auto negation_offset = elementwise_moduli<upper_bound>();
        return negation_offset - a;
    }

    template<int64_t lower, int64_t upper, class RNSBoundsType>
    INLINE auto negate(const SmallElement<BoundedElement<Bounds<lower, upper>, LIMBS, element_bits>, RNSBoundsType> &a) const {
        auto negation_offset = elementwise_moduli<upper>();
        auto m2_negation_offset = negation_offset.m2;
        auto ns = SmallElement<decltype(m2_negation_offset), decltype(negation_offset.rns_bounds)>(m2_negation_offset);
        return ns - a;
    }

    // not as efficient as delaying negation but useful for making fp12 have a consistent input and output type
    template<int64_t lower, int64_t upper>
    INLINE StandardElement standard_negate(const ExpandedElement<BoundedElement<Bounds<lower, upper>, LIMBS, element_bits>, Bounds<0, 1>> &a) const {
        auto negated = negate(a) + modulus_in_rns;
        return prep<negated.m1.bounds.lower, negated.m1.bounds.upper, Bounds<0, 1>, STANDARD_BOUND>(negated);
    }

    template<int64_t lower, int64_t upper>
    INLINE auto standard_negate(const ExpandedElement<BoundedElement<Bounds<lower, upper>, LIMBS, element_bits>, Bounds<0, 6>> &a) const {
        auto negated = negate(a) + six_modulus_in_rns;
        return prep<negated.m1.bounds.lower, negated.m1.bounds.upper, Bounds<0, 6>, STANDARD_BOUND>(negated);
    }

    template<int64_t lower, int64_t upper>
    INLINE SmallElement<StandardBoundedElement, Bounds<0, 1>> standard_negate(const SmallElement<BoundedElement<Bounds<lower, upper>, LIMBS, element_bits>, Bounds<0, 1>> &a) const {
        auto modulus_in_rns_m2 = modulus_in_rns.m2;
        auto modulus_in_rns_s = SmallElement<decltype(modulus_in_rns_m2), decltype(modulus_in_rns.rns_bounds)>(modulus_in_rns_m2);
        auto negated = negate(a) + modulus_in_rns_s;
        return prep<negated.m2.bounds.lower, negated.m2.bounds.upper, Bounds<0, 1>, STANDARD_BOUND>(negated);
    }

    // broken: TODO: debug
    // template<class BoundedElementType, int64_t rns_lower, int64_t rns_upper>
    // inline BigInt to_bigint(const ExpandedElement<BoundedElementType, Bounds<rns_lower, rns_upper>> &element) const {
    //     if constexpr(rns_lower < 0) {
    //         return to_bigint(convert_neg_offset + element);
    //     } else {
    //         static_assert(rns_lower >= 0, "RNS lower bound must be positive");
    //         static_assert(rns_upper <= RNS_MAX_MULTIPLE, "RNS upper bound must be less than or equal to RNS_MAX_MULTIPLE");
    //         static_assert((rns_upper * 3) * (2 * ceil_log2(limbs)) <= 1ull << (64 - 52), "Maximum to big int conversion bound exceeded.  Increase precision in convert to class.");
    //         auto p = prep(element);
    //         return sys.from_mont_avx(p.m2.data);
    //     }
    // }


    template<class ElementType, int64_t rns_lower, int64_t rns_upper>
    INLINE auto complete(const ExpandedElement<ElementType, Bounds<rns_lower, rns_upper>> &delayed_product) const {
        // automatically add negation offset if needed
        if constexpr(rns_lower < 0) {
            return negative_range_offset() + delayed_product;
        } else {
            // automatically complete delayed products if needed
            return delayed_product.complete();
        }
    }
    
    template<int64_t A, class ElementType, class RNSBoundsType>
    INLINE ReadyElement<A> ready(const ExpandedElement<ElementType, RNSBoundsType> &e) const {
        auto element = complete(e);
        using CompleteElementType = decltype(element);
        using LowType = decltype(element.m2.low);
        using HighType = decltype(element.m2.high);
        constexpr auto rns_lower_bound = CompleteElementType::rns_bounds.lower;
        constexpr auto rns_upper_bound = CompleteElementType::rns_bounds.upper;
        constexpr auto low_lower_bound = LowType::bounds.lower;
        constexpr auto high_lower_bound = HighType::bounds.lower;
        static_assert(rns_lower_bound >= 0, "RNS lower bound must be positive");
        static_assert(rns_upper_bound <= RNS_MAX_MULTIPLE, "RNS upper bound must be less than or equal to RNS_MAX_MULTIPLE");
        auto m2_recast_low = element.m2.low.template recast<0, A>();
        auto m2_recast_high = element.m2.high.template recast<0, A>();
        auto m2_recast = WideElement<Bounds<0, A>, LIMBS, element_bits>(m2_recast_high, m2_recast_low);
        static_assert(low_lower_bound >= 0);
        static_assert(high_lower_bound >= 0);
        if constexpr (change_base == ChangeBaseVariation::FixedPerm || change_base == ChangeBaseVariation::MatrixNoK) {
            auto m1_mont_ele_reduced = montgomery_reduce_m1.template reduce_auto<STANDARD_BOUND>(element.m1);
            auto m1_recast = m1_mont_ele_reduced.template recast<0, STANDARD_BOUND>();
            return ReadyToReduce<WideElement<Bounds<0, A>, LIMBS, element_bits>, StandardBoundedElement>(m2_recast, m1_recast);
        } else {
            auto m1_mont_ele_reduced = montgomery_reduce_m1.template reduce_auto<MULT_OK_BOUND>(element.m1);
            auto m1_recast = m1_mont_ele_reduced.template recast<0, MULT_OK_BOUND>();
            return ReadyToReduce<WideElement<Bounds<0, A>, LIMBS, element_bits>, MultOkBoundedElement>(m2_recast, m1_recast);
        }
    }

    // Expand, maintains RNS bound and allows for multiplications
    template<int batch, int64_t RNS_Bound>
    INLINE std::array<ExpandedElement<StandardBoundedElement, Bounds<0, RNS_Bound>>, batch>
    batch_expand(const std::array<SmallElement<StandardBoundedElement, Bounds<0, RNS_Bound>>, batch> &c_m2) const {
        auto wide_zero = WideZero<LIMBS, element_bits>();
        using ZeroAccElement = ReadyToReduce<decltype(wide_zero), StandardBoundedElement>;
        std::array<ZeroAccElement, batch> zero_acc_elements;
        for (int i = 0; i < batch; i++) {
            zero_acc_elements[i] = ZeroAccElement(wide_zero, c_m2[i].element);
        }
        auto c_m1_acc = batch_change_base_r2<batch>(zero_acc_elements);
        std::array<ExpandedElement<StandardBoundedElement, Bounds<0, RNS_Bound>>, batch> result;
        for (int i = 0; i < batch; i++) {
            result[i].m1 = montgomery_reduce_m1.template reduce_auto<STANDARD_BOUND>(c_m1_acc[i]);
            result[i].m2 = c_m2[i].element;
        }
        return result;
    }

    template<int batch>
    INLINE std::array<StandardBoundedElement, batch>
    batch_expand(const std::array<StandardBoundedElement, batch> &c_m2) const {
        auto wide_zero = WideZero<LIMBS, element_bits>();
        using ZeroAccElement = ReadyToReduce<decltype(wide_zero), StandardBoundedElement>;
        std::array<ZeroAccElement, batch> zero_acc_elements;
        for (int i = 0; i < batch; i++) {
            zero_acc_elements[i] = ZeroAccElement(wide_zero, c_m2[i]);
        }
        auto c_m1_acc = batch_change_base_r2<batch>(zero_acc_elements);
        std::array<StandardBoundedElement, batch> result;
        for (int i = 0; i < batch; i++) {
            result[i] = montgomery_reduce_m1.template reduce_auto<STANDARD_BOUND>(c_m1_acc[i]);
        }
        return result;
    }


    template<class E1, class E2, class E3, class E4, class E5, class E6>
    INLINE auto batch_expand(const E1 &e1, const E2 &e2, const E3 &e3, const E4 &e4, const E5 &e5, const E6 &e6) const {
        require_non_negative_rns<decltype(e1.rns_bounds)>();
        require_non_negative_rns<decltype(e2.rns_bounds)>();
        require_non_negative_rns<decltype(e3.rns_bounds)>();
        require_non_negative_rns<decltype(e4.rns_bounds)>();
        require_non_negative_rns<decltype(e5.rns_bounds)>();
        require_non_negative_rns<decltype(e6.rns_bounds)>();
        std::array<StandardBoundedElement, 6> data = {
            e1.element, e2.element, e3.element, e4.element, e5.element, e6.element};
        auto expanded = batch_expand<6>(data);
        // Reapply bounds
        auto e1_expanded = ExpandedElement<StandardBoundedElement, std::remove_cv_t<decltype(e1.rns_bounds)>>(expanded[0], e1);
        auto e2_expanded = ExpandedElement<StandardBoundedElement, std::remove_cv_t<decltype(e2.rns_bounds)>>(expanded[1], e2);
        auto e3_expanded = ExpandedElement<StandardBoundedElement, std::remove_cv_t<decltype(e3.rns_bounds)>>(expanded[2], e3);
        auto e4_expanded = ExpandedElement<StandardBoundedElement, std::remove_cv_t<decltype(e4.rns_bounds)>>(expanded[3], e4);
        auto e5_expanded = ExpandedElement<StandardBoundedElement, std::remove_cv_t<decltype(e5.rns_bounds)>>(expanded[4], e5);
        auto e6_expanded = ExpandedElement<StandardBoundedElement, std::remove_cv_t<decltype(e6.rns_bounds)>>(expanded[5], e6);
        return std::tuple<decltype(e1_expanded), decltype(e2_expanded), decltype(e3_expanded), decltype(e4_expanded), decltype(e5_expanded), decltype(e6_expanded)>{e1_expanded, e2_expanded, e3_expanded, e4_expanded, e5_expanded, e6_expanded};
    }

    // Reduce, sets RNS bound to 1
    template<int batch, int64_t A=MAX_ADD>
    INLINE auto batch_reduce(const std::array<ReadyElement<A>, batch> &c_in) const {
        if constexpr (change_base == ChangeBaseVariation::FixedPerm) {
            static_assert(TRUE_LIMBS == 6 && LIMBS == 8, "FixedPerm requires 6×8 ring");
            static_assert(MUL_3B_COEFF == 9, "FixedPerm requires bn254 (3b=9)");
            return bn254_batch_reduce_perm_bound<BoundedRing, batch, A>(*this, c_in);
        } else if constexpr (change_base == ChangeBaseVariation::MatrixNoK) {
            return matrix_nok::dispatch_batch_reduce_bound<BoundedRing, batch, A>(*this, c_in);
        } else {
            auto c_m2_acc = batch_change_base<batch>(c_in, sys.intrns2.r1);
            // We split reduce wide and reduce auto to enable scalar operations on unexpanded form.
            using MontHalfType = decltype(montgomery_reduce_m2.reduce_wide(c_m2_acc[0]));
            // Apply RNS Bounds manually since now RNS reduced
            using OutType = SmallElement<MontHalfType, Bounds<0, 1>>;
            std::array<OutType, batch> result;
            for (int i = 0; i < batch; i++) {
                result[i] = OutType(montgomery_reduce_m2.reduce_wide(c_m2_acc[i]));
            }
            return result;
        }
    }

    // One 3p step in RNS (modulus_in_rns); PointAdd uses k× this before prep/expand.
    INLINE SmallStandardElement small_modulus_in_rns() const {
        return SmallStandardElement(modulus_in_rns.m2.template recast<0, STANDARD_BOUND>());
    }

    template<class ElementType, class RNSBoundsType>
    INLINE auto prep_expand(const SmallElement<ElementType, RNSBoundsType> &e) const {
        if constexpr (RNSBoundsType::lower < 0) {
            // Add k×3p (value + RNS tag) so metadata is non-negative before batch_expand.
            static_assert(RNSBoundsType::lower >= -RNS_HALF_MAX_NEG, "RNS lower bound must be >= -RNS_HALF_MAX_NEG");
            if constexpr (RNSBoundsType::lower >= -4) {
                return prep_expand(small_rns_offset_4 + e);
            } else {
                return prep_expand(sub_rns_offset + e);
            }
        } else {
            return SmallElement<StandardBoundedElement, RNSBoundsType>(
                montgomery_reduce_m2.template reduce_auto<STANDARD_BOUND>(e.element));
        }
    }

    // Helper to do reduction followed by expansion.
    template<int batch>
    INLINE std::array<StandardElement, batch> batch_reduce_expand(const std::array<ReadyElement<MAX_ADD>, batch> &elements) const {
        if constexpr (change_base == ChangeBaseVariation::FixedPerm) {
            static_assert(TRUE_LIMBS == 6 && LIMBS == 8, "FixedPerm requires 6×8 ring");
            static_assert(MUL_3B_COEFF == 9, "FixedPerm requires bn254 (3b=9)");
            return bn254_batch_reduce_expand_perm<BoundedRing, batch>(*this, elements);
        } else if constexpr (change_base == ChangeBaseVariation::MatrixNoK) {
            return matrix_nok::dispatch_batch_reduce_expand<BoundedRing, batch>(*this, elements);
        }
        auto cm2_acc = batch_reduce<batch>(elements);
        std::array<SmallStandardElement, batch> c_m2;
        for (int i = 0; i < batch; i++) {
            c_m2[i] = SmallElement<StandardBoundedElement, Bounds<0, 1>>(montgomery_reduce_m2.template reduce_auto<STANDARD_BOUND>(cm2_acc[i].element));
        }
        return batch_expand<batch>(c_m2);
    }

    
    template<int batch>
    INLINE std::array<StandardElement, batch> from_bigint_batch(const std::array<BigInt, batch> &a) const {
        auto [c_m2_hi, c_m2_lo] = sys.convert_to.template batch_convert_to_rns<batch>(a);
        std::array<SmallStandardElement, batch> result;
        using WideBoundedElement = BoundedElement<Bounds<0, MAX_ADD>, LIMBS, element_bits>;
        for (int i = 0; i < batch; i++) {
            result[i] = SmallElement<StandardBoundedElement, Bounds<0, 1>>(montgomery_reduce_m2.template reduce_auto<STANDARD_BOUND>(WideElement<Bounds<0, MAX_ADD>, LIMBS, element_bits>(WideBoundedElement(c_m2_hi[i]), WideBoundedElement(c_m2_lo[i]))));
        }
        return batch_expand<batch>(result);
    }

    template<int batch>
    INLINE std::array<StandardElement, batch> convert_to_rns_batch(const std::array<AVXVector<LIMBS>, batch> &c_m2) const {
        std::array<AVXVector<LIMBS>, batch> c_m2_hi{};
        std::array<AVXVector<LIMBS>, batch> c_m2_lo{};
        sys.convert_to.conversion_matrix.template rns_reduce_raw_batch<batch>(&c_m2, c_m2_hi, c_m2_lo);
        std::array<SmallStandardElement, batch> result;
        using WideBoundedElement = BoundedElement<Bounds<0, MAX_ADD>, LIMBS, element_bits>;
        for (int i = 0; i < batch; i++) {
            result[i] = SmallElement<StandardBoundedElement, Bounds<0, 1>>(montgomery_reduce_m2.template reduce_auto<STANDARD_BOUND>(WideElement<Bounds<0, MAX_ADD>, LIMBS, element_bits>(WideBoundedElement(c_m2_hi[i]), WideBoundedElement(c_m2_lo[i]))));
        }
        return batch_expand<batch>(result);
    }
    
    INLINE StandardElement from_bigint(const BigInt &a) const {
        std::array<BigInt, 1> a_array = {a};
        return from_bigint_batch<1>(a_array)[0];
    }

    // Exact Montgomery encode/decode (bypasses fast convert_to/from path)
    INLINE StandardElement from_bigint_exact(const BigInt &a) const {
        return wrap(sys.to_mont_exact(a, false));
    }

    INLINE BigInt from_mont_exact_rns(
        const std::array<int64_t, LIMBS> &m1,
        const std::array<int64_t, LIMBS> &m2,
        bool wide) const {
        return std::get<0>(sys.from_mont_exact(m1, m2, wide));
    }

    // Expanded (m1 + m2 bounded limbs) — e.g. StandardElement after batch_reduce_expand.
    template<int64_t lower, int64_t upper, int64_t rns_lower, int64_t rns_upper>
    INLINE BigInt to_bigint_exact(
        const ExpandedElement<BoundedElement<Bounds<lower, upper>, LIMBS, element_bits>, Bounds<rns_lower, rns_upper>> &element,
        bool wide = false) const {
        (void)rns_lower;
        (void)rns_upper;
        return from_mont_exact_rns(element.m1.to_array(), element.m2.to_array(), wide);
    }

    INLINE BigInt to_bigint_exact(const StandardElement &element, bool wide = false) const {
        return to_bigint_exact<0, STANDARD_BOUND, 0, 1>(element, wide);
    }

    // Unexpanded: moduli2 only (after batch_reduce); recover m1 via batch_expand.
    template<class RNSBoundsType>
    INLINE BigInt to_bigint_exact(
        const SmallElement<StandardBoundedElement, RNSBoundsType> &element,
        bool wide = false) const {
        SmallElement<StandardBoundedElement, Bounds<0, 1>> norm(element.element);
        std::array<SmallElement<StandardBoundedElement, Bounds<0, 1>>, 1> arr = {norm};
        return to_bigint_exact(batch_expand<1, 1>(arr)[0], wide);
    }

    template<class ElementType, class RNSBoundsType>
    INLINE BigInt to_bigint_exact(const SmallElement<ElementType, RNSBoundsType> &element, bool wide = false) const {
        constexpr int64_t upper = ElementType::bounds.upper;
        constexpr int64_t lower = ElementType::bounds.lower;
        StandardBoundedElement m2_std =
            montgomery_reduce_m2.template reduce_auto<STANDARD_BOUND, upper, lower>(element.element);
        return to_bigint_exact(SmallElement<StandardBoundedElement, RNSBoundsType>(m2_std), wide);
    }

    // Wide expanded (high/low limbs); pass wide=true for full reconstruction.
    template<class WideElementType, class RNSBoundsType>
    INLINE BigInt to_bigint_exact(const ExpandedElement<WideElementType, RNSBoundsType> &element, bool wide = true) const {
        std::array<int64_t, LIMBS> m1_reduced;
        std::array<int64_t, LIMBS> m2_reduced;
        auto m1_arr = element.m1.low.to_array();
        auto m2_arr = element.m2.low.to_array();
        auto m1_high_arr = element.m1.high.to_array();
        auto m2_high_arr = element.m2.high.to_array();
        for (int i = 0; i < LIMBS; i++) {
            m1_reduced[i] = (((int128_t)(m1_high_arr[i]) << 52) + (int128_t)(m1_arr[i])) % sys.moduli1[i];
            m2_reduced[i] = (((int128_t)(m2_high_arr[i]) << 52) + (int128_t)(m2_arr[i])) % sys.moduli2[i];
        }
        return from_mont_exact_rns(m1_reduced, m2_reduced, wide);
    }

    // Helper to convert variadic arguments to array using index_sequence
    template<int64_t A, typename... Args, std::size_t... I>
    INLINE auto make_ready_array_impl(const std::tuple<Args...> &args, std::index_sequence<I...>) const {
        return std::array<ReadyElement<A>, sizeof...(Args)>{ready<A>(std::get<I>(args))...};
    }

    template<int64_t A, typename... Args, std::size_t... I>
    INLINE auto make_prep_array_impl(const std::tuple<Args...> &args, std::index_sequence<I...>) const {
        return std::array<SmallStandardElement, sizeof...(Args)>{prep(std::get<I>(args))...};
    }

     // Variadic template version that captures all the different arities
     template<int64_t A, typename... Args>
     INLINE auto batch_ready(const Args&... args) const {
         constexpr int batch = static_cast<int>(sizeof...(Args));
         static_assert(batch > 0, "batch_reduce_expand requires at least one argument");
         auto tuple = std::tie(args...);
         auto ready_array = make_ready_array_impl<A>(tuple, std::make_index_sequence<sizeof...(Args)>{});
         return ready_array;
    }

    // Variadic template version that captures all the different arities
    template<typename... Args>
    INLINE auto batch_reduce_expand(const Args&... args) const {
        constexpr int batch = static_cast<int>(sizeof...(Args));
        auto reduced = batch_ready<MAX_ADD>(args...);
        return batch_reduce_expand<batch>(reduced);
    }

    template<typename... Args>
    INLINE auto batch_reduce(const Args&... args) const {
        constexpr int batch = static_cast<int>(sizeof...(Args));
        auto reduced = batch_ready<MULT_OK_BOUND*MULT_OK_BOUND>(args...);
        return batch_reduce<batch>(reduced);
    }

    INLINE StandardElement modmul(const StandardElement &a, const StandardElement &b) const {
        auto product = (a * b);
        auto ready_product = ready<MAX_ADD>(product);
        std::array<ReadyElement<MAX_ADD>, 1> ready_elements = {ready_product};
        return batch_reduce_expand<1>(ready_elements)[0];
    }

    template<int batch>
    INLINE std::array<StandardElement, batch> batch_modmul(const std::array<StandardElement, batch> &a, const std::array<StandardElement, batch> &b) const {
        std::array<ReadyElement<MAX_ADD>, batch> ready_elements;
        for (int i = 0; i < batch; i++) {
            ready_elements[i] = ready<MAX_ADD>(a[i] * b[i]);
        }
        return batch_reduce_expand<batch>(ready_elements);
    }

    // Paper 3b·x; BLS12-381 b=4 → 12, bn254 b=3 → 9.
    template<class Element>
    INLINE auto mul_3b(const Element &x) const {
        if constexpr (MUL_3B_COEFF == 9) {
            if constexpr(element_bits <= 49) {
                return x.template mul_scalar<9>();
            } else {
                return x * std::integral_constant<int, 9>{};
            }
        }
        return x * std::integral_constant<int, MUL_3B_COEFF>{};
    }

    
    // Find a number matching the bounds.
    template<int ele_bound, int rns_bound>
    INLINE WideOffsetElement<ele_bound, rns_bound> compute_wide_offset() const {
        static_assert(ele_bound >= 0, "Element bound must be non-negative");
        static_assert(ele_bound <= 1ull << (64 - 52), "Element bound must be less than or equal to 2**(64 - 52)");
        BigInt desired_value = WIDE_REDUNDANT_MULTIPLE * rns_bound * modulus * modulus;
        auto desired_value_rns = sys.to_mont_exact(desired_value, true);
        std::array<uint64_t, LIMBS> val_m1_high;
        std::array<uint64_t, LIMBS> val_m2_high;
        std::array<uint64_t, LIMBS> val_m1_low;
        std::array<uint64_t, LIMBS> val_m2_low;
        uint64_t low_part_min = (1ull << 52) * ele_bound;
        uint64_t low_mask = (1ull << 52) - 1;
        for (int j = 0; j < LIMBS; j++) {
            uint128_t val_m1 = (uint128_t)(ele_bound) * (uint128_t)(sys.moduli1[j]) * (uint128_t)(sys.moduli1[j]);
            uint128_t val_m2 = (uint128_t)(ele_bound) * (uint128_t)(sys.moduli2[j]) * (uint128_t)(sys.moduli2[j]);
            val_m1_low[j] = val_m1 & low_mask;
            val_m2_low[j] = val_m2 & low_mask;
            val_m1_high[j] = val_m1 >> 52;
            val_m2_high[j] = val_m2 >> 52;

            val_m1_low[j] += desired_value_rns.first[j];
            val_m2_low[j] += desired_value_rns.second[j];
            // The low bits must be at least 2**52 * ele_bound to ensure the low part is positive when subtracting low and high without borrowing.  Add multiples of the modulus to make this happen
            val_m1_low[j] += sys.moduli1[j] * ele_bound;
            val_m2_low[j] += sys.moduli2[j] * ele_bound;
            // Because 2**52 - t < 2**52, may need to add a few more.
            while (val_m1_low[j] < low_part_min) {
                val_m1_low[j] += sys.moduli1[j];
            }
            while (val_m2_low[j] < low_part_min) {
                val_m2_low[j] += sys.moduli2[j];
            }
        }
        

        auto val_m1_high_avx = AVXVector<LIMBS>(val_m1_high);
        auto val_m2_high_avx = AVXVector<LIMBS>(val_m2_high);
        auto val_m1_low_avx = AVXVector<LIMBS>(val_m1_low);
        auto val_m2_low_avx = AVXVector<LIMBS>(val_m2_low);
        using ElementType = BoundedElement<Bounds<ele_bound, ele_bound + 1>, LIMBS, element_bits>;
        using WType = WideElement<Bounds<ele_bound, ele_bound + 1>, LIMBS, element_bits>;
        auto w1 = WType(ElementType(val_m1_high_avx), ElementType(val_m1_low_avx));
        auto w2 = WType(ElementType(val_m2_high_avx), ElementType(val_m2_low_avx));
        return WideOffsetElement<ele_bound, rns_bound>(w1, w2);
    }

    template<int64_t k>
    INLINE SmallRNSOffsetAt<k> compute_small_rns_offset() const {
        BigInt desired_value = BASE_REDUNDANT_MULTIPLE * k * modulus;
        auto desired_value_rns = sys.to_mont_exact(desired_value, false);
        std::array<uint64_t, LIMBS> val_m2;
        for (int i = 0; i < LIMBS; i++) {
            val_m2[i] = desired_value_rns.second[i];
        }
        auto val_m2_avx = AVXVector<LIMBS>(val_m2);
        return SmallRNSOffsetAt<k>(StandardBoundedElement(val_m2_avx));
    }

    template<int64_t lower, int64_t upper, int64_t rns_lower, int64_t rns_upper>
    INLINE auto check_bounds(const ExpandedElement<BoundedElement<Bounds<lower, upper>, LIMBS, element_bits>, Bounds<rns_lower, rns_upper>> &element, const std::string &context = "") const {
        bool ok = true;
        std::array<int64_t, LIMBS> moduli1;
        std::array<int64_t, LIMBS> moduli2;
        for (int i = 0; i < LIMBS; i++) {
            moduli1[i] = sys.moduli1[i];
            moduli2[i] = sys.moduli2[i];
        }

        if (!element.m1.check_bounds(moduli1, context)) {
            std::cout << context << ": M1 bounds do not match" << std::endl;
            ok = false;
        }
        if (!element.m2.check_bounds(moduli2, context)) {
            std::cout << context << ": M2 bounds do not match" << std::endl;
            ok = false;
        }
        auto result = sys.from_mont_exact(element.m1.to_array(), element.m2.to_array(), false);
        BigInt redundance = std::get<1>(result) / BASE_REDUNDANT_MULTIPLE;
        if (redundance < rns_lower) {
            std::cout << context << ": RNS redundance out of bounds too low: " << redundance << " < " << rns_lower << std::endl;
            ok = false;
        }
        if (redundance > rns_upper) {
            std::cout << context << ": RNS redundance out of bounds too high: " << redundance << " > " << rns_upper << std::endl;
            ok = false;
        }
        return std::make_pair(std::get<0>(result), ok);
    }


    template<class WideElementType, class RNSBoundsType>
    INLINE auto check_bounds(const ExpandedElement<WideElementType, RNSBoundsType> &element, const std::string &context = "") const {
        // Extract bounds from the types
        constexpr int64_t rns_lower = RNSBoundsType::lower;
        constexpr int64_t rns_upper = RNSBoundsType::upper;
        bool ok = true;
        std::array<int64_t, LIMBS> moduli1;
        std::array<int64_t, LIMBS> moduli2;
        for (int i = 0; i < LIMBS; i++) {
            moduli1[i] = sys.moduli1[i];
            moduli2[i] = sys.moduli2[i];
        }
        if (!element.m1.check_bounds(moduli1, context + ".m1")) {
            std::cout << context << ": M1 bounds do not match" << std::endl;
            ok = false;
        }
        if (!element.m2.check_bounds(moduli2, context + ".m2")) {
            std::cout << context << ": M2 bounds do not match" << std::endl;
            ok = false;
        }
        std::array<int64_t, LIMBS> m1_reduced;
        std::array<int64_t, LIMBS> m2_reduced;
        auto m1_arr = element.m1.low.to_array();
        auto m2_arr = element.m2.low.to_array();
        auto m1_high_arr = element.m1.high.to_array();
        auto m2_high_arr = element.m2.high.to_array();
        for (int i = 0; i < LIMBS; i++) {
            m1_reduced[i] = (((int128_t)(m1_high_arr[i]) << 52) + (int128_t)(m1_arr[i])) % sys.moduli1[i];
            m2_reduced[i] = (((int128_t)(m2_high_arr[i]) << 52) + (int128_t)(m2_arr[i])) % sys.moduli2[i];
        }
        auto result = sys.from_mont_exact(m1_reduced, m2_reduced, true);
        auto RNSVal = std::get<0>(result);
        auto RNSRedundance = std::get<1>(result);
        // For wide elements representing offset * modulus^2, from_mont_exact returns redundance = offset * modulus
        // So we need to divide by target to get the offset
        BigInt offset_value = RNSRedundance / (WIDE_REDUNDANT_MULTIPLE * sys.target);
        if (offset_value < rns_lower) {
            std::cout << context << ": RNS redundance out of bounds too low: " << RNSRedundance << " / " << WIDE_REDUNDANT_MULTIPLE * sys.target << " = " << offset_value << " < " << rns_lower << std::endl;
            ok = false;
        }
        if (offset_value > rns_upper) {
            std::cout << context << ": RNS redundance out of bounds too high: " << RNSRedundance << " / " << WIDE_REDUNDANT_MULTIPLE * sys.target << " = " << offset_value << " > " << rns_upper << std::endl;
            ok = false;
        }
        return std::make_pair(RNSVal, ok);
    }

    template<class RNSMatrix>
    INLINE bool check_matrix_bound(const RNSMatrix& rns_matrix, std::array<int64_t, LIMBS> &moduli) const {
        std::array<uint64_t, LIMBS> max_input;
        for (int i = 0; i < LIMBS; i++) {
            max_input[i] = (1ull << 52) - 1;
        }
        auto wide_zero_val = WideZero<LIMBS, element_bits>();
        // Recast wide_zero to match ReadyElement bounds - ReadyElement uses MAX_ADD bounds
        using ReadyWideType = WideElement<Bounds<0, MAX_ADD>, LIMBS, element_bits>;
        auto wide_zero_recast = ReadyWideType(
            BoundedElement<Bounds<0, MAX_ADD>, LIMBS, element_bits>(wide_zero_val.low.data),
            BoundedElement<Bounds<0, MAX_ADD>, LIMBS, element_bits>(wide_zero_val.high.data)
        );
        // ReadyElement expects StandardBoundedElement (bounds [0, standard_bound]) for the second parameter
        auto max_input_bounded = StandardBoundedElement(AVXVector<LIMBS>(max_input));
        auto ready_elem = ReadyElement<MAX_ADD>(wide_zero_recast, max_input_bounded);
        auto result = batch_change_base<1>(std::array<ReadyElement<MAX_ADD>, 1> {ready_elem}, rns_matrix)[0];
        return result.check_bounds(moduli, "check_matrix_bound result");
    }


    INLINE bool check_all() const {
        bool ok = true;

        std::array<int64_t, LIMBS> moduli1;
        std::array<int64_t, LIMBS> moduli2;
        std::array<int64_t, LIMBS> max_value_array;
        for (int i = 0; i < LIMBS; i++) {
            moduli1[i] = sys.moduli1[i];
            moduli2[i] = sys.moduli2[i];
            max_value_array[i] = (1ull << 52) - 1;
        }
        ok &= check_matrix_bound(sys.intrns2.r1, moduli2);
        ok &= check_matrix_bound(sys.intrns2.r2, moduli1);
        ok &= check_matrix_bound(sys.convert_to.conversion_matrix, moduli2);
        // Convert from doesn't do modular reduction
        ok &= check_matrix_bound(sys.convert_from.conversion_matrix, max_value_array);
        
        // Check negative_range_offset value is 0
        auto neg_offset_result = check_bounds(negative_range_offset(), "negative_range_offset");
        ok &= (neg_offset_result.first == 0) && neg_offset_result.second;
        
        // Check elementwise_moduli values are 0
        auto moduli1_result = check_bounds(elementwise_moduli<1>(), "elementwise_moduli<1>");
        ok &= (moduli1_result.first == 0) && moduli1_result.second;
        
        auto moduli2_result = check_bounds(elementwise_moduli<2>(), "elementwise_moduli<2>");
        ok &= (moduli2_result.first == 0) && moduli2_result.second;
        
        auto moduli4_result = check_bounds(elementwise_moduli<4>(), "elementwise_moduli<4>");
        ok &= (moduli4_result.first == 0) && moduli4_result.second;
        
        return ok;
    }

private:
    template<int batch, class ZeroAccElement>
    INLINE auto batch_change_base_r2(const std::array<ZeroAccElement, batch> &zero_acc_elements) const {
        if constexpr (change_base == ChangeBaseVariation::FixedPerm) {
            static_assert(TRUE_LIMBS == 6 && LIMBS == 8, "FixedPerm requires 6×8 ring");
            static_assert(MUL_3B_COEFF == 9, "FixedPerm requires bn254 (3b=9)");
            return batch_change_base_perm<batch>(zero_acc_elements, bn254_r1_perm::r2_fixed_perm_wide_reducer);
        }
        return batch_change_base<batch>(zero_acc_elements, sys.intrns2.r2);
    }
};
