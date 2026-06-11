#include "bn254_ring_params.hpp"
#include "ec.hpp"
#include "../cpu/precompute/gmp_wrapper.hpp"
#include "test_data/test_values_bn254_ec.hpp"
#include <cassert>
#include <iostream>

// BN254 G1 EC tests: IntRNS4 fixed perm wide, 6×46 residues, b=3.
using Bn254Ring = bn254_params::Bn254RingProduction;
using Bn254G1 = G1<Bn254Ring>;

static void assert_affine_on_curve(const char* label, int case_idx, const BigInt& x, const BigInt& y, const BigInt& p) {
    BigInt yy = (y * y) % p;
    BigInt rhs = (x * x * x) % p;
    rhs = (rhs + BigInt(bn254_params::BN254_B)) % p;
    if (yy != rhs) {
        std::cerr << "case " << case_idx << ": " << label << " not on Y²=X³+" << bn254_params::BN254_B << "\n"
                  << "  (" << x.to_string() << ", " << y.to_string() << ")\n";
        assert(false);
    }
}

template<class Ring>
AffinePoint<BigInt> to_bigint_g1_exact(
    const ProjectivePoint<typename Ring::StandardElement>& point,
    const Ring& ring,
    const BigInt& modulus) {
    BigInt x = ring.to_bigint_exact(point.x);
    BigInt y = ring.to_bigint_exact(point.y);
    BigInt z = ring.to_bigint_exact(point.z);
    BigInt z_inv = z.mod_inverse(modulus);
    BigInt affine_x = (x * z_inv) % modulus;
    if (affine_x < 0) affine_x = affine_x + modulus;
    BigInt affine_y = (y * z_inv) % modulus;
    if (affine_y < 0) affine_y = affine_y + modulus;
    return AffinePoint<BigInt>{affine_x, affine_y};
}

static void check_affine(
    int case_idx,
    const char* label,
    const AffinePoint<BigInt>& got,
    const char* exp_x,
    const char* exp_y) {
    BigInt ex(exp_x);
    BigInt ey(exp_y);
    if (got.x != ex || got.y != ey) {
        std::cerr << "case " << case_idx << " " << label << " mismatch\n"
                  << "  got: (" << got.x.to_string() << ", " << got.y.to_string() << ")\n"
                  << "  exp: (" << ex.to_string() << ", " << ey.to_string() << ")\n";
        assert(false);
    }
}

template<class Ring>
void run_case(int case_idx, const Bn254EcCase& c, const Ring& ring, const Bn254G1& curve, const BigInt& p) {
    BigInt p_x(c.p_x), p_y(c.p_y), q_x(c.q_x), q_y(c.q_y);
    assert_affine_on_curve("P", case_idx, p_x, p_y, p);
    assert_affine_on_curve("Q", case_idx, q_x, q_y, p);
    assert_affine_on_curve("add expected", case_idx, BigInt(c.add_x), BigInt(c.add_y), p);
    assert_affine_on_curve("mixed expected", case_idx, BigInt(c.mix_x), BigInt(c.mix_y), p);
    assert_affine_on_curve("double expected", case_idx, BigInt(c.dbl_x), BigInt(c.dbl_y), p);

    auto px = ring.from_bigint_exact(p_x);
    auto py = ring.from_bigint_exact(p_y);
    auto qx = ring.from_bigint_exact(q_x);
    auto qy = ring.from_bigint_exact(q_y);

    using Proj = typename Bn254G1::ProjPoint;
    using Aff = typename Bn254G1::AffPoint;

    Proj P(px, py, ring.one());
    Proj Q(qx, qy, ring.one());
    Aff Qaff(qx, qy);

    check_affine(case_idx, "add", to_bigint_g1_exact(curve.add_point(P, Q, ring), ring, p), c.add_x, c.add_y);
    check_affine(case_idx, "mixed", to_bigint_g1_exact(curve.add_mixed_point(P, Qaff, ring), ring, p), c.mix_x, c.mix_y);
    check_affine(case_idx, "double", to_bigint_g1_exact(curve.double_point(P, ring), ring, p), c.dbl_x, c.dbl_y);

    BigInt pz_x(c.pz_x), pz_y(c.pz_y), q2_x(c.q2_x), q2_y(c.q2_y);
    assert_affine_on_curve("Pz", case_idx, pz_x, pz_y, p);
    assert_affine_on_curve("Q2", case_idx, q2_x, q2_y, p);
    assert_affine_on_curve("mixed_z expected", case_idx, BigInt(c.mz_x), BigInt(c.mz_y), p);
    auto pzx = ring.from_bigint_exact(pz_x);
    auto pzy = ring.from_bigint_exact(pz_y);
    Aff Q2aff(ring.from_bigint_exact(q2_x), ring.from_bigint_exact(q2_y));
    Proj Pz(pzx, pzy, ring.one());
    check_affine(
        case_idx,
        "mixed_z",
        to_bigint_g1_exact(curve.add_mixed_point(Pz, Q2aff, ring), ring, p),
        c.mz_x,
        c.mz_y);
}

namespace test_values_bn254 {
#include "test_data/test_values_bn254.hpp"
}

void test_field_mul_exact(const Bn254Ring& ring) {
    BigInt a(test_values_bn254::a);
    BigInt b(test_values_bn254::b);
    BigInt c_exp(test_values_bn254::c_correct);

    auto ra = ring.from_bigint_exact(a);
    auto rb = ring.from_bigint_exact(b);
    std::array<typename Bn254Ring::StandardElement, 1> as = {ra};
    std::array<typename Bn254Ring::StandardElement, 1> bs = {rb};
    BigInt c_exact = ring.to_bigint_exact(ring.batch_modmul<1>(as, bs)[0]);
    if (c_exact != c_exp) {
        std::cerr << "field mul (exact conv) mismatch\n  got " << c_exact.to_string()
                  << "\n  exp " << c_exp.to_string() << "\n";
        assert(false);
    }

    auto ra_fast = ring.from_bigint(a);
    auto rb_fast = ring.from_bigint(b);
    std::array<typename Bn254Ring::StandardElement, 1> as_fast = {ra_fast};
    std::array<typename Bn254Ring::StandardElement, 1> bs_fast = {rb_fast};
    BigInt c_fast = ring.to_bigint(ring.batch_modmul<1>(as_fast, bs_fast)[0]);
    if (c_fast != c_exp) {
        std::cerr << "field mul (fast conv) mismatch\n  got " << c_fast.to_string()
                  << "\n  exp " << c_exp.to_string() << "\n";
        assert(false);
    }

    std::cout << "field mul OK (exact and fast convert agree)" << std::endl;
}

int main() {
    BigInt p(bn254_params::BN254_P);
    Bn254Ring ring(p);
    test_field_mul_exact(ring);
    Bn254G1 curve;

    for (int i = 0; i < EC_NUM_CASES; i++) {
        run_case(i, EC_CASES[i], ring, curve, p);
    }
    std::cout << "BN254 G1 EC (FixedPerm / IntRNS4): " << EC_NUM_CASES
              << " cases passed (6×46, b=3)" << std::endl;
    return 0;
}
