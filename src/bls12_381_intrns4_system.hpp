#pragma once

#include "../cpu/vector/multiplication.hpp"
#include "../scripts/bls12_381_intrns4_params.hpp"

/** Build IntRNS2System from Python IntRNS4 BLS12-381 export (8×50-bit moduli, qr rotation). */
template <int bits, int padded_limbs>
inline IntRNS2System<bits, padded_limbs, MontgomeryReduce, bls12_381_intrns4::BLS12_381_CONVERT_TO_RADIX_BITS>
make_BLS12_381_IntRNS4SystemFromExport(const BigInt &target) {
    using namespace bls12_381_intrns4;
    static_assert(padded_limbs == BLS12_381_PADDED_DIM, "IntRNS4 BLS12-381 export uses 8-wide moduli");
    const BigInt M1_big(M1);
    const BigInt M2_big(M2);
    const BigInt M1M2_big(M1M2);
    const BigInt total_rotation_big(total_rotation);
    const BigInt total_inv_rotation_big(total_inv_rotation);
    const BigInt inv_factor_target_big(inv_factor_target);

    return IntRNS2System<bits, padded_limbs, MontgomeryReduce, bls12_381_intrns4::BLS12_381_CONVERT_TO_RADIX_BITS>(
        IntRNS2<MontgomeryReduce, padded_limbs, padded_limbs>(
            bls12_381_intrns4::moduli1_padded,
            bls12_381_intrns4::moduli2_padded,
            target,
            BigInt(r1_premult),
            BigInt(r1_postmult),
            BigInt(r2_premult),
            BigInt(r2_postmult)),
        ConvertToRNS<bits, padded_limbs, bls12_381_intrns4::BLS12_381_CONVERT_TO_RADIX_BITS>(
            BigInt(convert_to_premult),
            BigInt(convert_to_postmult),
            target,
            bls12_381_intrns4::moduli2_padded),
        ConvertFromRNS<bits, padded_limbs>(
            BigInt(convert_from_premult),
            BigInt(convert_from_postmult),
            target,
            bls12_381_intrns4::moduli2_padded),
        bls12_381_intrns4::moduli1_padded,
        bls12_381_intrns4::moduli2_padded,
        target,
        M1_big,
        M2_big,
        M1M2_big,
        total_rotation_big,
        total_inv_rotation_big,
        inv_factor_target_big);
}
