#pragma once
#include <array>
#include <cstdint>
/* 
Polynomial 0(22) : 1 x0y0 -1 x1y1 +1 x2y4 -1 x2y5 -1 x3y4 -1 x3y5 +1 x4y2 -1 x4y3 -1 x5y2 -1 x5y3 +1 x6y10 -1 x6y11 -1 x7y10 -1 x7y11 +1 x8y8 -1 x8y9 -1 x9y8 -1 x9y9 +1 x10y6 -1 x10y7 -1 x11y6 -1 x11y7

Polynomial 1(22) : 1 x0y1 +1 x1y0 +1 x2y4 +1 x2y5 +1 x3y4 -1 x3y5 +1 x4y2 +1 x4y3 +1 x5y2 -1 x5y3 +1 x6y10 +1 x6y11 +1 x7y10 -1 x7y11 +1 x8y8 +1 x8y9 +1 x9y8 -1 x9y9 +1 x10y6 +1 x10y7 +1 x11y6 -1 x11y7

Polynomial 2(18) : 1 x0y2 -1 x1y3 +1 x2y0 -1 x3y1 +1 x4y4 -1 x4y5 -1 x5y4 -1 x5y5 +1 x6y6 -1 x7y7 +1 x8y10 -1 x8y11 -1 x9y10 -1 x9y11 +1 x10y8 -1 x10y9 -1 x11y8 -1 x11y9

Polynomial 3(18) : 1 x0y3 +1 x1y2 +1 x2y1 +1 x3y0 +1 x4y4 +1 x4y5 +1 x5y4 -1 x5y5 +1 x6y7 +1 x7y6 +1 x8y10 +1 x8y11 +1 x9y10 -1 x9y11 +1 x10y8 +1 x10y9 +1 x11y8 -1 x11y9

Polynomial 4(14) : 1 x0y4 -1 x1y5 +1 x2y2 -1 x3y3 +1 x4y0 -1 x5y1 +1 x6y8 -1 x7y9 +1 x8y6 -1 x9y7 +1 x10y10 -1 x10y11 -1 x11y10 -1 x11y11

Polynomial 5(14) : 1 x0y5 +1 x1y4 +1 x2y3 +1 x3y2 +1 x4y1 +1 x5y0 +1 x6y9 +1 x7y8 +1 x8y7 +1 x9y6 +1 x10y10 +1 x10y11 +1 x11y10 -1 x11y11

Polynomial 6(20) : 1 x0y6 -1 x1y7 +1 x2y10 -1 x2y11 -1 x3y10 -1 x3y11 +1 x4y8 -1 x4y9 -1 x5y8 -1 x5y9 +1 x6y0 -1 x7y1 +1 x8y4 -1 x8y5 -1 x9y4 -1 x9y5 +1 x10y2 -1 x10y3 -1 x11y2 -1 x11y3

Polynomial 7(20) : 1 x0y7 +1 x1y6 +1 x2y10 +1 x2y11 +1 x3y10 -1 x3y11 +1 x4y8 +1 x4y9 +1 x5y8 -1 x5y9 +1 x6y1 +1 x7y0 +1 x8y4 +1 x8y5 +1 x9y4 -1 x9y5 +1 x10y2 +1 x10y3 +1 x11y2 -1 x11y3

Polynomial 8(16) : 1 x0y8 -1 x1y9 +1 x2y6 -1 x3y7 +1 x4y10 -1 x4y11 -1 x5y10 -1 x5y11 +1 x6y2 -1 x7y3 +1 x8y0 -1 x9y1 +1 x10y4 -1 x10y5 -1 x11y4 -1 x11y5

Polynomial 9(16) : 1 x0y9 +1 x1y8 +1 x2y7 +1 x3y6 +1 x4y10 +1 x4y11 +1 x5y10 -1 x5y11 +1 x6y3 +1 x7y2 +1 x8y1 +1 x9y0 +1 x10y4 +1 x10y5 +1 x11y4 -1 x11y5

Polynomial 10(12) : 1 x0y10 -1 x1y11 +1 x2y8 -1 x3y9 +1 x4y6 -1 x5y7 +1 x6y4 -1 x7y5 +1 x8y2 -1 x9y3 +1 x10y0 -1 x11y1

Polynomial 11(12) : 1 x0y11 +1 x1y10 +1 x2y9 +1 x3y8 +1 x4y7 +1 x5y6 +1 x6y5 +1 x7y4 +1 x8y3 +1 x9y2 +1 x10y1 +1 x11y0
 */
template<class Ring>
inline auto fp12_flat_mul(const std::array<typename Ring::StandardElement, 12> &x, const std::array<typename Ring::StandardElement, 12> &y, Ring &ring) {
	auto minus_y1 = ring.negate(y[1]);
	auto minus_y3 = ring.negate(y[3]);
	auto minus_y5 = ring.negate(y[5]);
	auto minus_y7 = ring.negate(y[7]);
	auto minus_y9 = ring.negate(y[9]);
	auto minus_y11 = ring.negate(y[11]);
	// 0 1
	auto poly_0_pos_shared = offset + x[2] * y[4] + x[3] * minus_y5 + x[4] * y[2] + x[5] * minus_y3 + x[6] * y[10] + x[7] * minus_y11 + x[8] * y[8] + x[9] * minus_y9 + x[10] * y[6] + x[11] * minus_y7;
	auto poly_0_neg_shared = x[2] * y[5] + x[3] * y[4] + x[4] * y[3] + x[5] * y[2] + x[6] * y[11] + x[7] * y[10] + x[8] * y[9] + x[9] * y[8] + x[10] * y[7] + x[11] * y[6];
	auto poly_0_unique = poly_0_pos_shared - poly_0_neg_shared + x[0] * y[0] + x[1] * minus_y1;
	auto poly_1_unique = poly_0_pos_shared + poly_0_neg_shared + x[0] * y[1] + x[1] * y[0];
	// 2 3
	auto poly_2_pos_shared = offset + x[4] * y[4] + x[5] * minus_y5 + x[8] * y[10] + x[9] * minus_y11 + x[10] * y[8] + x[11] * minus_y9;
	auto poly_2_neg_shared = x[4] * y[5] + x[5] * y[4] + x[8] * y[11] + x[9] * y[10] + x[10] * y[9] + x[11] * y[8];
	auto poly_2_unique = poly_2_pos_shared - poly_2_neg_shared + x[0] * y[2] + x[1] * minus_y3 + x[2] * y[0] + x[3] * minus_y1 + x[6] * y[6] + x[7] * minus_y7;
	auto poly_3_unique = poly_2_pos_shared + poly_2_neg_shared + x[0] * y[3] + x[1] * y[2] + x[2] * y[1] + x[3] * y[0] + x[6] * y[7] + x[7] * y[6];
	// 4 5
	auto poly_4_pos_shared = offset + x[10] * y[10] + x[11] * minus_y11;
	auto poly_4_neg_shared = x[10] * y[11] + x[11] * y[10];
	auto poly_4_unique = poly_4_pos_shared - poly_4_neg_shared + x[0] * y[4] + x[1] * minus_y5 + x[2] * y[2] + x[3] * minus_y3 + x[4] * y[0] + x[5] * minus_y1 + x[6] * y[8] + x[7] * minus_y9 + x[8] * y[6] + x[9] * minus_y7;
	auto poly_5_unique = poly_4_pos_shared + poly_4_neg_shared + x[0] * y[5] + x[1] * y[4] + x[2] * y[3] + x[3] * y[2] + x[4] * y[1] + x[5] * y[0] + x[6] * y[9] + x[7] * y[8] + x[8] * y[7] + x[9] * y[6];
	// 6 7
	auto poly_6_pos_shared = offset + x[2] * y[10] + x[3] * minus_y11 + x[4] * y[8] + x[5] * minus_y9 + x[8] * y[4] + x[9] * minus_y5 + x[10] * y[2] + x[11] * minus_y3;
	auto poly_6_neg_shared = x[2] * y[11] + x[3] * y[10] + x[4] * y[9] + x[5] * y[8] + x[8] * y[5] + x[9] * y[4] + x[10] * y[3] + x[11] * y[2];
	auto poly_6_unique = poly_6_pos_shared - poly_6_neg_shared + x[0] * y[6] + x[1] * minus_y7 + x[6] * y[0] + x[7] * minus_y1;
	auto poly_7_unique = poly_6_pos_shared + poly_6_neg_shared + x[0] * y[7] + x[1] * y[6] + x[6] * y[1] + x[7] * y[0];
	// 8 9
	auto poly_8_pos_shared = offset + x[4] * y[10] + x[5] * minus_y11 + x[10] * y[4] + x[11] * minus_y5;
	auto poly_8_neg_shared = x[4] * y[11] + x[5] * y[10] + x[10] * y[5] + x[11] * y[4];
	auto poly_8_unique = poly_8_pos_shared - poly_8_neg_shared + x[0] * y[8] + x[1] * minus_y9 + x[2] * y[6] + x[3] * minus_y7 + x[6] * y[2] + x[7] * minus_y3 + x[8] * y[0] + x[9] * minus_y1;
	auto poly_9_unique = poly_8_pos_shared + poly_8_neg_shared + x[0] * y[9] + x[1] * y[8] + x[2] * y[7] + x[3] * y[6] + x[6] * y[3] + x[7] * y[2] + x[8] * y[1] + x[9] * y[0];
	// 10 11
	auto poly_10_unique = x[0] * y[10] + x[1] * minus_y11 + x[2] * y[8] + x[3] * minus_y9 + x[4] * y[6] + x[5] * minus_y7 + x[6] * y[4] + x[7] * minus_y5 + x[8] * y[2] + x[9] * minus_y3 + x[10] * y[0] + x[11] * minus_y1;
	auto poly_11_unique = x[0] * y[11] + x[1] * y[10] + x[2] * y[9] + x[3] * y[8] + x[4] * y[7] + x[5] * y[6] + x[6] * y[5] + x[7] * y[4] + x[8] * y[3] + x[9] * y[2] + x[10] * y[1] + x[11] * y[0];
	return ring.template batch_reduce_expand(poly_0_unique,poly_1_unique,poly_2_unique,poly_3_unique,poly_4_unique,poly_5_unique,poly_6_unique,poly_7_unique,poly_8_unique,poly_9_unique,poly_10_unique,poly_11_unique);
}
/* 
Polynomial 0(10) : 1 x0y0 -1 x1y1 +1 x4y2 -1 x4y3 -1 x5y2 -1 x5y3 +1 x8y4 -1 x8y5 -1 x9y4 -1 x9y5

Polynomial 1(10) : 1 x0y1 +1 x1y0 +1 x4y2 +1 x4y3 +1 x5y2 -1 x5y3 +1 x8y4 +1 x8y5 +1 x9y4 -1 x9y5

Polynomial 2(8) : 1 x0y2 -1 x1y3 +1 x2y0 -1 x3y1 +1 x10y4 -1 x10y5 -1 x11y4 -1 x11y5

Polynomial 3(8) : 1 x0y3 +1 x1y2 +1 x2y1 +1 x3y0 +1 x10y4 +1 x10y5 +1 x11y4 -1 x11y5

Polynomial 4(6) : 1 x2y2 -1 x3y3 +1 x4y0 -1 x5y1 +1 x6y4 -1 x7y5

Polynomial 5(6) : 1 x2y3 +1 x3y2 +1 x4y1 +1 x5y0 +1 x6y5 +1 x7y4

Polynomial 6(10) : 1 x4y4 -1 x4y5 -1 x5y4 -1 x5y5 +1 x6y0 -1 x7y1 +1 x10y2 -1 x10y3 -1 x11y2 -1 x11y3

Polynomial 7(10) : 1 x4y4 +1 x4y5 +1 x5y4 -1 x5y5 +1 x6y1 +1 x7y0 +1 x10y2 +1 x10y3 +1 x11y2 -1 x11y3

Polynomial 8(6) : 1 x0y4 -1 x1y5 +1 x6y2 -1 x7y3 +1 x8y0 -1 x9y1

Polynomial 9(6) : 1 x0y5 +1 x1y4 +1 x6y3 +1 x7y2 +1 x8y1 +1 x9y0

Polynomial 10(6) : 1 x2y4 -1 x3y5 +1 x8y2 -1 x9y3 +1 x10y0 -1 x11y1

Polynomial 11(6) : 1 x2y5 +1 x3y4 +1 x8y3 +1 x9y2 +1 x10y1 +1 x11y0
 */
template<class Ring>
inline auto flat_mulxy0z00(const std::array<typename Ring::StandardElement, 12> &x, const std::array<typename Ring::StandardElement, 6> &y, Ring &ring) {
	auto minus_y1 = ring.negate(y[1]);
	auto minus_y3 = ring.negate(y[3]);
	auto minus_y5 = ring.negate(y[5]);
	// 0 1
	auto poly_0_pos_shared = offset + x[4] * y[2] + x[5] * minus_y3 + x[8] * y[4] + x[9] * minus_y5;
	auto poly_0_neg_shared = x[4] * y[3] + x[5] * y[2] + x[8] * y[5] + x[9] * y[4];
	auto poly_0_unique = poly_0_pos_shared - poly_0_neg_shared + x[0] * y[0] + x[1] * minus_y1;
	auto poly_1_unique = poly_0_pos_shared + poly_0_neg_shared + x[0] * y[1] + x[1] * y[0];
	// 2 3
	auto poly_2_pos_shared = offset + x[10] * y[4] + x[11] * minus_y5;
	auto poly_2_neg_shared = x[10] * y[5] + x[11] * y[4];
	auto poly_2_unique = poly_2_pos_shared - poly_2_neg_shared + x[0] * y[2] + x[1] * minus_y3 + x[2] * y[0] + x[3] * minus_y1;
	auto poly_3_unique = poly_2_pos_shared + poly_2_neg_shared + x[0] * y[3] + x[1] * y[2] + x[2] * y[1] + x[3] * y[0];
	// 4 5
	auto poly_4_unique = x[2] * y[2] + x[3] * minus_y3 + x[4] * y[0] + x[5] * minus_y1 + x[6] * y[4] + x[7] * minus_y5;
	auto poly_5_unique = x[2] * y[3] + x[3] * y[2] + x[4] * y[1] + x[5] * y[0] + x[6] * y[5] + x[7] * y[4];
	// 6 7
	auto poly_6_pos_shared = offset + x[4] * y[4] + x[5] * minus_y5 + x[10] * y[2] + x[11] * minus_y3;
	auto poly_6_neg_shared = x[4] * y[5] + x[5] * y[4] + x[10] * y[3] + x[11] * y[2];
	auto poly_6_unique = poly_6_pos_shared - poly_6_neg_shared + x[6] * y[0] + x[7] * minus_y1;
	auto poly_7_unique = poly_6_pos_shared + poly_6_neg_shared + x[6] * y[1] + x[7] * y[0];
	// 8 9
	auto poly_8_unique = x[0] * y[4] + x[1] * minus_y5 + x[6] * y[2] + x[7] * minus_y3 + x[8] * y[0] + x[9] * minus_y1;
	auto poly_9_unique = x[0] * y[5] + x[1] * y[4] + x[6] * y[3] + x[7] * y[2] + x[8] * y[1] + x[9] * y[0];
	// 10 11
	auto poly_10_unique = x[2] * y[4] + x[3] * minus_y5 + x[8] * y[2] + x[9] * minus_y3 + x[10] * y[0] + x[11] * minus_y1;
	auto poly_11_unique = x[2] * y[5] + x[3] * y[4] + x[8] * y[3] + x[9] * y[2] + x[10] * y[1] + x[11] * y[0];
	return ring.template batch_reduce_expand(poly_0_unique,poly_1_unique,poly_2_unique,poly_3_unique,poly_4_unique,poly_5_unique,poly_6_unique,poly_7_unique,poly_8_unique,poly_9_unique,poly_10_unique,poly_11_unique);
}
/* 
Polynomial 0(13) : 1 x0x0 -1 x1x1 +1 x8x8 -1 x9x9 +2 x2x4 -2 x2x5 -2 x3x4 -2 x3x5 +2 x6x10 -2 x6x11 -2 x7x10 -2 x7x11 -2 x8x9

Polynomial 1(12) : 1 x8x8 -1 x9x9 +2 x0x1 +2 x2x4 +2 x2x5 +2 x3x4 -2 x3x5 +2 x6x10 +2 x6x11 +2 x7x10 -2 x7x11 +2 x8x9

Polynomial 2(11) : 1 x4x4 -1 x5x5 +1 x6x6 -1 x7x7 +2 x0x2 -2 x1x3 -2 x4x5 +2 x8x10 -2 x8x11 -2 x9x10 -2 x9x11

Polynomial 3(10) : 1 x4x4 -1 x5x5 +2 x0x3 +2 x1x2 +2 x4x5 +2 x6x7 +2 x8x10 +2 x8x11 +2 x9x10 -2 x9x11

Polynomial 4(9) : 1 x2x2 -1 x3x3 +1 x10x10 -1 x11x11 +2 x0x4 -2 x1x5 +2 x6x8 -2 x7x9 -2 x10x11

Polynomial 5(8) : 1 x10x10 -1 x11x11 +2 x0x5 +2 x1x4 +2 x2x3 +2 x6x9 +2 x7x8 +2 x10x11

Polynomial 6(10) : 2 x0x6 -2 x1x7 +2 x2x10 -2 x2x11 -2 x3x10 -2 x3x11 +2 x4x8 -2 x4x9 -2 x5x8 -2 x5x9

Polynomial 7(10) : 2 x0x7 +2 x1x6 +2 x2x10 +2 x2x11 +2 x3x10 -2 x3x11 +2 x4x8 +2 x4x9 +2 x5x8 -2 x5x9

Polynomial 8(8) : 2 x0x8 -2 x1x9 +2 x2x6 -2 x3x7 +2 x4x10 -2 x4x11 -2 x5x10 -2 x5x11

Polynomial 9(8) : 2 x0x9 +2 x1x8 +2 x2x7 +2 x3x6 +2 x4x10 +2 x4x11 +2 x5x10 -2 x5x11

Polynomial 10(6) : 2 x0x10 -2 x1x11 +2 x2x8 -2 x3x9 +2 x4x6 -2 x5x7

Polynomial 11(6) : 2 x0x11 +2 x1x10 +2 x2x9 +2 x3x8 +2 x4x7 +2 x5x6
 */
template<class Ring>
inline auto fp12_flat_sqr(const std::array<typename Ring::StandardElement, 12> &x, Ring &ring) {
	auto minus_x1 = ring.negate(x[1]);
	auto minus_x3 = ring.negate(x[3]);
	auto minus_x5 = ring.negate(x[5]);
	auto minus_x7 = ring.negate(x[7]);
	auto minus_x9 = ring.negate(x[9]);
	auto minus_x11 = ring.negate(x[11]);
	// 0 1
	auto poly_0_pos_shared = offset + ((x[2] * x[4] + x[3] * minus_x5 + x[6] * x[10] + x[7] * minus_x11) * 2) + x[8] * x[8] + x[9] * minus_x9;
	auto poly_0_neg_shared = (x[2] * x[5] + x[3] * x[4] + x[6] * x[11] + x[7] * x[10] + x[8] * x[9]) * 2;
	auto poly_0_unique = poly_0_pos_shared - poly_0_neg_shared + x[0] * x[0] + x[1] * minus_x1;
	auto poly_1_unique = poly_0_pos_shared + poly_0_neg_shared + (x[0] * x[1]) * 2;
	// 2 3
	auto poly_2_pos_shared = offset + ((x[8] * x[10] + x[9] * minus_x11) * 2) + x[4] * x[4] + x[5] * minus_x5;
	auto poly_2_neg_shared = (x[4] * x[5] + x[8] * x[11] + x[9] * x[10]) * 2;
	auto poly_2_unique = poly_2_pos_shared - poly_2_neg_shared + ((x[0] * x[2] + x[1] * minus_x3) * 2) + x[6] * x[6] + x[7] * minus_x7;
	auto poly_3_unique = poly_2_pos_shared + poly_2_neg_shared + (x[0] * x[3] + x[1] * x[2] + x[6] * x[7]) * 2;
	// 4 5
	auto poly_4_pos_shared = offset + x[10] * x[10] + x[11] * minus_x11;
	auto poly_4_unique = poly_4_pos_shared + ((x[0] * x[4] + x[1] * minus_x5 + x[6] * x[8] + x[7] * minus_x9 + x[10] * minus_x11) * 2) + x[2] * x[2] + x[3] * minus_x3;
	auto poly_5_unique = poly_4_pos_shared + (x[0] * x[5] + x[1] * x[4] + x[2] * x[3] + x[6] * x[9] + x[7] * x[8] + x[10] * x[11]) * 2;
	// 6 7
	auto poly_6_pos_shared = offset + x[2] * x[10] + x[3] * minus_x11 + x[4] * x[8] + x[5] * minus_x9;
	auto poly_6_neg_shared = x[2] * x[11] + x[3] * x[10] + x[4] * x[9] + x[5] * x[8];
	auto poly_6_unique = (poly_6_pos_shared - poly_6_neg_shared + x[0] * x[6] + x[1] * minus_x7) * 2;
	auto poly_7_unique = (poly_6_pos_shared + poly_6_neg_shared + x[0] * x[7] + x[1] * x[6]) * 2;
	// 8 9
	auto poly_8_pos_shared = offset + x[4] * x[10] + x[5] * minus_x11;
	auto poly_8_neg_shared = x[4] * x[11] + x[5] * x[10];
	auto poly_8_unique = (poly_8_pos_shared - poly_8_neg_shared + x[0] * x[8] + x[1] * minus_x9 + x[2] * x[6] + x[3] * minus_x7) * 2;
	auto poly_9_unique = (poly_8_pos_shared + poly_8_neg_shared + x[0] * x[9] + x[1] * x[8] + x[2] * x[7] + x[3] * x[6]) * 2;
	// 10 11
	auto poly_10_unique = (x[0] * x[10] + x[1] * minus_x11 + x[2] * x[8] + x[3] * minus_x9 + x[4] * x[6] + x[5] * minus_x7) * 2;
	auto poly_11_unique = (x[0] * x[11] + x[1] * x[10] + x[2] * x[9] + x[3] * x[8] + x[4] * x[7] + x[5] * x[6]) * 2;
	return ring.template batch_reduce_expand(poly_0_unique,poly_1_unique,poly_2_unique,poly_3_unique,poly_4_unique,poly_5_unique,poly_6_unique,poly_7_unique,poly_8_unique,poly_9_unique,poly_10_unique,poly_11_unique);
}
/* 
Polynomial 0(6) : 3 x0x0 -2 x0 -3 x1x1 +3 x8x8 -3 x9x9 -6 x8x9

Polynomial 1(5) : -2 x1 +3 x8x8 -3 x9x9 +6 x0x1 +6 x8x9

Polynomial 2(6) : -2 x2 +3 x4x4 -3 x5x5 +3 x6x6 -3 x7x7 -6 x4x5

Polynomial 3(5) : -2 x3 +3 x4x4 -3 x5x5 +6 x4x5 +6 x6x7

Polynomial 4(6) : 3 x2x2 -3 x3x3 -2 x4 +3 x10x10 -3 x11x11 -6 x10x11

Polynomial 5(5) : -2 x5 +3 x10x10 -3 x11x11 +6 x2x3 +6 x10x11

Polynomial 6(5) : 2 x6 +6 x2x10 -6 x2x11 -6 x3x10 -6 x3x11

Polynomial 7(5) : 2 x7 +6 x2x10 +6 x2x11 +6 x3x10 -6 x3x11

Polynomial 8(3) : 2 x8 +6 x0x8 -6 x1x9

Polynomial 9(3) : 2 x9 +6 x0x9 +6 x1x8

Polynomial 10(3) : 2 x10 +6 x4x6 -6 x5x7

Polynomial 11(3) : 2 x11 +6 x4x7 +6 x5x6
 */
template<class Ring>
inline auto fp12_flat_cyclotomic_sqr(const std::array<typename Ring::StandardElement, 12> &x, Ring &ring) {
	auto minus_x1 = ring.negate(x[1]);
	auto minus_x3 = ring.negate(x[3]);
	auto minus_x5 = ring.negate(x[5]);
	auto minus_x7 = ring.negate(x[7]);
	auto minus_x9 = ring.negate(x[9]);
	auto minus_x11 = ring.negate(x[11]);
	auto two_thirds = ring.wrap(ring.template from_hex("0x8ab05f8bdd54cde190937e76bc3e447cc27c3d6fbd7063fcd104635a790520c0a395554e5c6aaaa9354ffffffffe38f"));
	auto neg_two_thirds = ring.wrap(ring.template from_hex("0x11560bf17baa99bc32126fced787c88f984f87adf7ae0c7f9a208c6b4f20a4181472aaa9cb8d555526a9ffffffffc71c"));
	// 0 1
	auto poly_0_pos_shared = offset + x[8] * x[8] + x[9] * minus_x9;
	auto poly_0_unique = (poly_0_pos_shared + ((x[8] * minus_x9) * 2) + x[0] * x[0] + x[0] * neg_two_thirds + x[1] * minus_x1) * 3;
	auto poly_1_unique = (poly_0_pos_shared + ((x[0] * x[1] + x[8] * x[9]) * 2) + x[1] * neg_two_thirds) * 3;
	// 2 3
	auto poly_2_pos_shared = offset + x[4] * x[4] + x[5] * minus_x5;
	auto poly_2_unique = (poly_2_pos_shared + ((x[4] * minus_x5) * 2) + x[2] * neg_two_thirds + x[6] * x[6] + x[7] * minus_x7) * 3;
	auto poly_3_unique = (poly_2_pos_shared + ((x[6] * x[7] + x[4] * x[5]) * 2) + x[3] * neg_two_thirds) * 3;
	// 4 5
	auto poly_4_pos_shared = offset + x[10] * x[10] + x[11] * minus_x11;
	auto poly_4_unique = (poly_4_pos_shared + ((x[10] * minus_x11) * 2) + x[2] * x[2] + x[3] * minus_x3 + x[4] * neg_two_thirds) * 3;
	auto poly_5_unique = (poly_4_pos_shared + ((x[2] * x[3] + x[10] * x[11]) * 2) + x[5] * neg_two_thirds) * 3;
	// 6 7
	auto poly_6_pos_shared = offset + (x[2] * x[10] + x[3] * minus_x11) * 2;
	auto poly_6_neg_shared = (x[2] * x[11] + x[3] * x[10]) * 2;
	auto poly_6_unique = (poly_6_pos_shared - poly_6_neg_shared + x[6] * two_thirds) * 3;
	auto poly_7_unique = (poly_6_pos_shared + poly_6_neg_shared + x[7] * two_thirds) * 3;
	// 8 9
	auto poly_8_unique = (((x[0] * x[8] + x[1] * minus_x9) * 2) + x[8] * two_thirds) * 3;
	auto poly_9_unique = (((x[0] * x[9] + x[1] * x[8]) * 2) + x[9] * two_thirds) * 3;
	// 10 11
	auto poly_10_unique = (((x[4] * x[6] + x[5] * minus_x7) * 2) + x[10] * two_thirds) * 3;
	auto poly_11_unique = (((x[4] * x[7] + x[5] * x[6]) * 2) + x[11] * two_thirds) * 3;
	return ring.template batch_reduce_expand(poly_0_unique,poly_1_unique,poly_2_unique,poly_3_unique,poly_4_unique,poly_5_unique,poly_6_unique,poly_7_unique,poly_8_unique,poly_9_unique,poly_10_unique,poly_11_unique);
}
template<class Ring>
inline auto fp12_flat_frobenius_map1(Ring &ring, const std::array<typename Ring::StandardElement, 12> &x) {
	auto frob1_2 = ring.wrap(ring.from_hex("0x1a0111ea397fe699ec02408663d4de85aa0d857d89759ad4897d29650fb85f9b409427eb4f49fffd8bfd00000000aaac"));
	auto frob1_2_2 = x[2] * frob1_2;
	auto frob1_2_3 = x[3] * frob1_2;
	auto frob1_3 = ring.wrap(ring.from_hex("0x1a0111ea397fe699ec02408663d4de85aa0d857d89759ad4897d29650fb85f9b409427eb4f49fffd8bfd00000000aaad"));
	auto frob1_3_4 = x[4] * frob1_3;
	auto frob1_4 = ring.wrap(ring.from_hex("0x5f19672fdf76ce51ba69c6076a0f77eaddb3a93be6f89688de17d813620a00022e01fffffffefffe"));
	auto frob1_4_5 = x[5] * frob1_4;
	auto frob1_5 = ring.wrap(ring.from_hex("0x1904d3bf02bb0667c231beb4202c0d1f0fd603fd3cbd5f4f7b2443d784bab9c4f67ea53d63e7813d8d0775ed92235fb8"));
	auto frob1_5_6 = x[6] * frob1_5;
	auto frob1_6 = ring.wrap(ring.from_hex("0xfc3e2b36c4e03288e9e902231f9fb854a14787b6c7b36fec0c8ec971f63c5f282d5ac14d6c7ec22cf78a126ddc4af3"));
	auto frob1_6_6 = x[6] * frob1_6;
	auto frob1_6_7 = x[7] * frob1_6;
	auto frob1_7 = ring.wrap(ring.from_hex("0x6af0e0437ff400b6831e36d6bd17ffe48395dabc2d3435e77f76e17009241c5ee67992f72ec05f4c81084fbede3cc09"));
	auto frob1_7_8 = x[8] * frob1_7;
	auto frob1_7_9 = x[9] * frob1_7;
	auto frob1_8 = ring.wrap(ring.from_hex("0x135203e60180a68ee2e9c448d77a2cd91c3dedd930b1cf60ef396489f61eb45e304466cf3e67fa0af1ee7b04121bdea2"));
	auto frob1_8_9 = x[9] * frob1_8;
	auto frob1_9 = ring.wrap(ring.from_hex("0x5b2cfd9013a5fd8df47fa6b48b1e045f39816240c0b8fee8beadf4d8e9c0566c63a3e6e257f87329b18fae980078116"));
	auto frob1_9_10 = x[10] * frob1_9;
	auto frob1_10 = ring.wrap(ring.from_hex("0x144e4211384586c16bd3ad4afa99cc9170df3560e77982d0db45f3536814f0bd5871c1908bd478cd1ee605167ff82995"));
	auto frob1_10_10 = x[10] * frob1_10;
	auto frob1_10_11 = x[11] * frob1_10;
	auto poly2_unique = frob1_2_3;
	auto poly3_unique = frob1_2_2;
	auto poly4_unique = frob1_3_4;
	auto poly5_unique = frob1_4_5;
	auto shared_frob1_6_7 = frob1_6_7;
	auto shared_frob1_7_8 = frob1_7_8;
	auto shared_frob1_10_11 = frob1_10_11;
	auto poly6_unique = shared_frob1_6_7 + frob1_5_6;
	auto poly7_unique = shared_frob1_6_7 + frob1_6_6;
	auto poly8_unique = shared_frob1_7_8 + frob1_7_9;
	auto poly9_unique = shared_frob1_7_8 + frob1_8_9;
	auto poly10_unique = shared_frob1_10_11 + frob1_9_10;
	auto poly11_unique = shared_frob1_10_11 + frob1_10_10;
	auto [reduced_poly2, reduced_poly3, reduced_poly4, reduced_poly5, reduced_poly6, reduced_poly7, reduced_poly8, reduced_poly9, reduced_poly10, reduced_poly11] = ring.template batch_reduce_expand(poly2_unique, poly3_unique, poly4_unique, poly5_unique, poly6_unique, poly7_unique, poly8_unique, poly9_unique, poly10_unique, poly11_unique);
	return std::array<typename Ring::StandardElement, 12>{x[0], ring.negate(x[1]), reduced_poly2, reduced_poly3, reduced_poly4, reduced_poly5, reduced_poly6, reduced_poly7, reduced_poly8, reduced_poly9, reduced_poly10, reduced_poly11};
}
template<class Ring>
inline auto fp12_flat_frobenius_map2(Ring &ring, const std::array<typename Ring::StandardElement, 12> &x) {
	auto frob2_1 = ring.wrap(ring.from_hex("0x5f19672fdf76ce51ba69c6076a0f77eaddb3a93be6f89688de17d813620a00022e01fffffffefffe"));
	auto frob2_1_2 = x[2] * frob2_1;
	auto frob2_1_3 = x[3] * frob2_1;
	auto frob2_2 = ring.wrap(ring.from_hex("0x1a0111ea397fe699ec02408663d4de85aa0d857d89759ad4897d29650fb85f9b409427eb4f49fffd8bfd00000000aaac"));
	auto frob2_2_4 = x[4] * frob2_2;
	auto frob2_2_5 = x[5] * frob2_2;
	auto frob2_3 = ring.wrap(ring.from_hex("0x5f19672fdf76ce51ba69c6076a0f77eaddb3a93be6f89688de17d813620a00022e01fffffffeffff"));
	auto frob2_3_6 = x[6] * frob2_3;
	auto frob2_3_7 = x[7] * frob2_3;
	auto frob2_5 = ring.wrap(ring.from_hex("0x1a0111ea397fe699ec02408663d4de85aa0d857d89759ad4897d29650fb85f9b409427eb4f49fffd8bfd00000000aaad"));
	auto frob2_5_10 = x[10] * frob2_5;
	auto frob2_5_11 = x[11] * frob2_5;
	auto poly2_unique = frob2_1_2;
	auto poly3_unique = frob2_1_3;
	auto poly4_unique = frob2_2_4;
	auto poly5_unique = frob2_2_5;
	auto poly6_unique = frob2_3_6;
	auto poly7_unique = frob2_3_7;
	auto poly10_unique = frob2_5_10;
	auto poly11_unique = frob2_5_11;
	auto [reduced_poly2, reduced_poly3, reduced_poly4, reduced_poly5, reduced_poly6, reduced_poly7, reduced_poly10, reduced_poly11] = ring.template batch_reduce_expand(poly2_unique, poly3_unique, poly4_unique, poly5_unique, poly6_unique, poly7_unique, poly10_unique, poly11_unique);
	return std::array<typename Ring::StandardElement, 12>{x[0], x[1], reduced_poly2, reduced_poly3, reduced_poly4, reduced_poly5, reduced_poly6, reduced_poly7, ring.negate(x[8]), ring.negate(x[9]), reduced_poly10, reduced_poly11};
}
template<class Ring>
inline auto fp12_flat_frobenius_map3(Ring &ring, const std::array<typename Ring::StandardElement, 12> &x) {
	auto frob3_2 = ring.wrap(ring.from_hex("0x135203e60180a68ee2e9c448d77a2cd91c3dedd930b1cf60ef396489f61eb45e304466cf3e67fa0af1ee7b04121bdea2"));
	auto frob3_2_6 = x[6] * frob3_2;
	auto frob3_2_8 = x[8] * frob3_2;
	auto frob3_2_9 = x[9] * frob3_2;
	auto frob3_2_10 = x[10] * frob3_2;
	auto frob3_2_11 = x[11] * frob3_2;
	auto frob3_3 = ring.wrap(ring.from_hex("0x6af0e0437ff400b6831e36d6bd17ffe48395dabc2d3435e77f76e17009241c5ee67992f72ec05f4c81084fbede3cc09"));
	auto frob3_3_9 = x[9] * frob3_3;
	auto frob3_3_10 = x[10] * frob3_3;
	auto frob3_3_6 = x[6] * frob3_3;
	auto frob3_3_7 = x[7] * frob3_3;
	auto shared_frob3_3_7 = frob3_3_7;
	auto shared_frob3_2_8 = frob3_2_8;
	auto shared_frob3_2_11 = frob3_2_11;
	auto poly6_unique = shared_frob3_3_7 + frob3_2_6;
	auto poly7_unique = shared_frob3_3_7 + frob3_3_6;
	auto poly8_unique = shared_frob3_2_8 + frob3_2_9;
	auto poly9_unique = shared_frob3_2_8 + frob3_3_9;
	auto poly10_unique = shared_frob3_2_11 + frob3_3_10;
	auto poly11_unique = shared_frob3_2_11 + frob3_2_10;
	auto [reduced_poly6, reduced_poly7, reduced_poly8, reduced_poly9, reduced_poly10, reduced_poly11] = ring.template batch_reduce_expand(poly6_unique, poly7_unique, poly8_unique, poly9_unique, poly10_unique, poly11_unique);
	return std::array<typename Ring::StandardElement, 12>{x[0], ring.negate(x[1]), x[3], x[2], ring.negate(x[4]), x[5], reduced_poly6, reduced_poly7, reduced_poly8, reduced_poly9, reduced_poly10, reduced_poly11};
}
template<class Ring>
inline auto fp12_flat_conjugate(const std::array<typename Ring::StandardElement, 12> &x, Ring &ring) {
	return std::array<typename Ring::StandardElement, 12>{x[0], x[1], x[2], x[3], x[4], x[5], ring.negate(x[6]), ring.negate(x[7]), ring.negate(x[8]), ring.negate(x[9]), ring.negate(x[10]), ring.negate(x[11])};
}
/* 
Polynomial 0(13) : 1 x0x0 -1 x1x1 -1 x8x8 +1 x9x9 +2 x2x4 -2 x2x5 -2 x3x4 -2 x3x5 -2 x6x10 +2 x6x11 +2 x7x10 +2 x7x11 +2 x8x9

Polynomial 1(12) : -1 x8x8 +1 x9x9 +2 x0x1 +2 x2x4 +2 x2x5 +2 x3x4 -2 x3x5 -2 x6x10 -2 x6x11 -2 x7x10 +2 x7x11 -2 x8x9

Polynomial 2(11) : 1 x4x4 -1 x5x5 -1 x6x6 +1 x7x7 +2 x0x2 -2 x1x3 -2 x4x5 -2 x8x10 +2 x8x11 +2 x9x10 +2 x9x11

Polynomial 3(10) : 1 x4x4 -1 x5x5 +2 x0x3 +2 x1x2 +2 x4x5 -2 x6x7 -2 x8x10 -2 x8x11 -2 x9x10 +2 x9x11

Polynomial 4(9) : 1 x2x2 -1 x3x3 -1 x10x10 +1 x11x11 +2 x0x4 -2 x1x5 -2 x6x8 +2 x7x9 +2 x10x11

Polynomial 5(8) : -1 x10x10 +1 x11x11 +2 x0x5 +2 x1x4 +2 x2x3 -2 x6x9 -2 x7x8 -2 x10x11
 */
template<class Ring>
inline auto fp12_flat_inverse_to_fp6(const std::array<typename Ring::StandardElement, 12> &x, Ring &ring) {
	auto minus_x1 = ring.negate(x[1]);
	auto minus_x3 = ring.negate(x[3]);
	auto minus_x5 = ring.negate(x[5]);
	auto minus_x6 = ring.negate(x[6]);
	auto minus_x7 = ring.negate(x[7]);
	auto minus_x8 = ring.negate(x[8]);
	auto minus_x9 = ring.negate(x[9]);
	auto minus_x10 = ring.negate(x[10]);
	auto minus_x11 = ring.negate(x[11]);
	// 0 1
	auto poly_0_pos_shared = offset + ((x[2] * x[4] + x[3] * minus_x5 + x[6] * minus_x10 + x[7] * x[11]) * 2) + x[8] * minus_x8 + x[9] * x[9];
	auto poly_0_neg_shared = (x[2] * x[5] + x[3] * x[4] + x[6] * minus_x11 + x[7] * minus_x10 + x[8] * minus_x9) * 2;
	auto poly_0_unique = poly_0_pos_shared - poly_0_neg_shared + x[0] * x[0] + x[1] * minus_x1;
	auto poly_1_unique = poly_0_pos_shared + poly_0_neg_shared + (x[0] * x[1]) * 2;
	// 2 3
	auto poly_2_pos_shared = offset + ((x[8] * minus_x10 + x[9] * x[11]) * 2) + x[4] * x[4] + x[5] * minus_x5;
	auto poly_2_neg_shared = (x[4] * x[5] + x[8] * minus_x11 + x[9] * minus_x10) * 2;
	auto poly_2_unique = poly_2_pos_shared - poly_2_neg_shared + ((x[0] * x[2] + x[1] * minus_x3) * 2) + x[6] * minus_x6 + x[7] * x[7];
	auto poly_3_unique = poly_2_pos_shared + poly_2_neg_shared + (x[0] * x[3] + x[1] * x[2] + x[6] * minus_x7) * 2;
	// 4 5
	auto poly_4_pos_shared = offset + x[10] * minus_x10 + x[11] * x[11];
	auto poly_4_unique = poly_4_pos_shared + ((x[0] * x[4] + x[1] * minus_x5 + x[6] * minus_x8 + x[7] * x[9] + x[10] * x[11]) * 2) + x[2] * x[2] + x[3] * minus_x3;
	auto poly_5_unique = poly_4_pos_shared + (x[0] * x[5] + x[1] * x[4] + x[2] * x[3] + x[6] * minus_x9 + x[7] * minus_x8 + x[10] * minus_x11) * 2;
	return ring.template batch_reduce_expand(poly_0_unique,poly_1_unique,poly_2_unique,poly_3_unique,poly_4_unique,poly_5_unique);
}
/* 
Polynomial 0(10) : 1 x0y0 -1 x1y1 +1 x2y4 -1 x2y5 -1 x3y4 -1 x3y5 +1 x4y2 -1 x4y3 -1 x5y2 -1 x5y3

Polynomial 1(10) : 1 x0y1 +1 x1y0 +1 x2y4 +1 x2y5 +1 x3y4 -1 x3y5 +1 x4y2 +1 x4y3 +1 x5y2 -1 x5y3

Polynomial 2(8) : 1 x0y2 -1 x1y3 +1 x2y0 -1 x3y1 +1 x4y4 -1 x4y5 -1 x5y4 -1 x5y5

Polynomial 3(8) : 1 x0y3 +1 x1y2 +1 x2y1 +1 x3y0 +1 x4y4 +1 x4y5 +1 x5y4 -1 x5y5

Polynomial 4(6) : 1 x0y4 -1 x1y5 +1 x2y2 -1 x3y3 +1 x4y0 -1 x5y1

Polynomial 5(6) : 1 x0y5 +1 x1y4 +1 x2y3 +1 x3y2 +1 x4y1 +1 x5y0

Polynomial 6(10) : -1 x6y0 +1 x7y1 -1 x8y4 +1 x8y5 +1 x9y4 +1 x9y5 -1 x10y2 +1 x10y3 +1 x11y2 +1 x11y3

Polynomial 7(10) : -1 x6y1 -1 x7y0 -1 x8y4 -1 x8y5 -1 x9y4 +1 x9y5 -1 x10y2 -1 x10y3 -1 x11y2 +1 x11y3

Polynomial 8(8) : -1 x6y2 +1 x7y3 -1 x8y0 +1 x9y1 -1 x10y4 +1 x10y5 +1 x11y4 +1 x11y5

Polynomial 9(8) : -1 x6y3 -1 x7y2 -1 x8y1 -1 x9y0 -1 x10y4 -1 x10y5 -1 x11y4 +1 x11y5

Polynomial 10(6) : -1 x6y4 +1 x7y5 -1 x8y2 +1 x9y3 -1 x10y0 +1 x11y1

Polynomial 11(6) : -1 x6y5 -1 x7y4 -1 x8y3 -1 x9y2 -1 x10y1 -1 x11y0
 */
template<class Ring>
inline auto fp6_flat_inverse_to_fp12(const std::array<typename Ring::StandardElement, 12> &x, const std::array<typename Ring::StandardElement, 6> &y, Ring &ring) {
	auto minus_x0 = ring.negate(x[0]);
	auto minus_x1 = ring.negate(x[1]);
	auto minus_x2 = ring.negate(x[2]);
	auto minus_x3 = ring.negate(x[3]);
	auto minus_x4 = ring.negate(x[4]);
	auto minus_x5 = ring.negate(x[5]);
	// 0 1
	auto poly_0_pos_shared = offset + x[2] * x[4] + x[3] * minus_x5 + x[4] * x[2] + x[5] * minus_x3;
	auto poly_0_neg_shared = x[2] * x[5] + x[3] * x[4] + x[4] * x[3] + x[5] * x[2];
	auto poly_0_unique = poly_0_pos_shared - poly_0_neg_shared + x[0] * x[0] + x[1] * minus_x1;
	auto poly_1_unique = poly_0_pos_shared + poly_0_neg_shared + x[0] * x[1] + x[1] * x[0];
	// 2 3
	auto poly_2_pos_shared = offset + x[4] * x[4] + x[5] * minus_x5;
	auto poly_2_neg_shared = x[4] * x[5] + x[5] * x[4];
	auto poly_2_unique = poly_2_pos_shared - poly_2_neg_shared + x[0] * x[2] + x[1] * minus_x3 + x[2] * x[0] + x[3] * minus_x1;
	auto poly_3_unique = poly_2_pos_shared + poly_2_neg_shared + x[0] * x[3] + x[1] * x[2] + x[2] * x[1] + x[3] * x[0];
	// 4 5
	auto poly_4_unique = x[0] * x[4] + x[1] * minus_x5 + x[2] * x[2] + x[3] * minus_x3 + x[4] * x[0] + x[5] * minus_x1;
	auto poly_5_unique = x[0] * x[5] + x[1] * x[4] + x[2] * x[3] + x[3] * x[2] + x[4] * x[1] + x[5] * x[0];
	// 6 7
	auto poly_6_pos_shared = offset + x[8] * minus_x4 + x[9] * x[5] + x[10] * minus_x2 + x[11] * x[3];
	auto poly_6_neg_shared = x[8] * minus_x5 + x[9] * minus_x4 + x[10] * minus_x3 + x[11] * minus_x2;
	auto poly_6_unique = poly_6_pos_shared - poly_6_neg_shared + x[6] * minus_x0 + x[7] * x[1];
	auto poly_7_unique = poly_6_pos_shared + poly_6_neg_shared + x[6] * minus_x1 + x[7] * minus_x0;
	// 8 9
	auto poly_8_pos_shared = offset + x[10] * minus_x4 + x[11] * x[5];
	auto poly_8_neg_shared = x[10] * minus_x5 + x[11] * minus_x4;
	auto poly_8_unique = poly_8_pos_shared - poly_8_neg_shared + x[6] * minus_x2 + x[7] * x[3] + x[8] * minus_x0 + x[9] * x[1];
	auto poly_9_unique = poly_8_pos_shared + poly_8_neg_shared + x[6] * minus_x3 + x[7] * minus_x2 + x[8] * minus_x1 + x[9] * minus_x0;
	// 10 11
	auto poly_10_unique = x[6] * minus_x4 + x[7] * x[5] + x[8] * minus_x2 + x[9] * x[3] + x[10] * minus_x0 + x[11] * x[1];
	auto poly_11_unique = x[6] * minus_x5 + x[7] * minus_x4 + x[8] * minus_x3 + x[9] * minus_x2 + x[10] * minus_x1 + x[11] * minus_x0;
	return ring.template batch_reduce_expand(poly_0_unique,poly_1_unique,poly_2_unique,poly_3_unique,poly_4_unique,poly_5_unique,poly_6_unique,poly_7_unique,poly_8_unique,poly_9_unique,poly_10_unique,poly_11_unique);
}
/* 
Polynomial 0(6) : 1 x0x0 -1 x1x1 -1 x2x4 +1 x2x5 +1 x3x4 +1 x3x5

Polynomial 1(5) : 2 x0x1 -1 x2x4 -1 x2x5 -1 x3x4 +1 x3x5

Polynomial 2(5) : 1 x4x4 -1 x5x5 -1 x0x2 +1 x1x3 -2 x4x5

Polynomial 3(5) : 1 x4x4 -1 x5x5 -1 x0x3 -1 x1x2 +2 x4x5

Polynomial 4(4) : 1 x2x2 -1 x3x3 -1 x0x4 +1 x1x5

Polynomial 5(3) : -1 x0x5 -1 x1x4 +2 x2x3
 */
template<class Ring>
inline auto fp6_flat_inverse_to_fp2_c0_c1_c2(const std::array<typename Ring::StandardElement, 6> &x, Ring &ring) {
	auto minus_x1 = ring.negate(x[1]);
	auto minus_x2 = ring.negate(x[2]);
	auto minus_x3 = ring.negate(x[3]);
	auto minus_x4 = ring.negate(x[4]);
	auto minus_x5 = ring.negate(x[5]);
	// 0 1
	auto poly_0_pos_shared = offset + x[2] * minus_x4 + x[3] * x[5];
	auto poly_0_neg_shared = x[2] * minus_x5 + x[3] * minus_x4;
	auto poly_0_unique = poly_0_pos_shared - poly_0_neg_shared + x[0] * x[0] + x[1] * minus_x1;
	auto poly_1_unique = poly_0_pos_shared + poly_0_neg_shared + (x[0] * x[1]) * 2;
	// 2 3
	auto poly_2_pos_shared = offset + x[4] * x[4] + x[5] * minus_x5;
	auto poly_2_unique = poly_2_pos_shared + ((x[4] * minus_x5) * 2) + x[0] * minus_x2 + x[1] * x[3];
	auto poly_3_unique = poly_2_pos_shared + ((x[4] * x[5]) * 2) + x[0] * minus_x3 + x[1] * minus_x2;
	// 4 5
	auto poly_4_unique = x[2] * x[2] + x[3] * minus_x3 + x[0] * minus_x4 + x[1] * x[5];
	auto poly_5_unique = ((x[2] * x[3]) * 2) + x[0] * minus_x5 + x[1] * minus_x4;
	return ring.template batch_reduce_expand(poly_0_unique,poly_1_unique,poly_2_unique,poly_3_unique,poly_4_unique,poly_5_unique);
}
/* 
Polynomial 0(10) : 1 x0x6 -1 x1x7 +1 x2x10 -1 x2x11 -1 x3x10 -1 x3x11 +1 x4x8 -1 x4x9 -1 x5x8 -1 x5x9

Polynomial 1(10) : 1 x0x7 +1 x1x6 +1 x2x10 +1 x2x11 +1 x3x10 -1 x3x11 +1 x4x8 +1 x4x9 +1 x5x8 -1 x5x9
 */
template<class Ring>
inline auto fp6_flat_inverse_to_fp2_ret(const std::array<typename Ring::StandardElement, 6> &x, Ring &ring) {
	auto minus_x1 = ring.negate(x[1]);
	auto minus_x7 = ring.negate(x[7]);
	auto minus_x9 = ring.negate(x[9]);
	auto minus_x11 = ring.negate(x[11]);
	// 0 1
	auto poly_0_pos_shared = offset + x[2] * x[10] + x[3] * minus_x11 + x[4] * x[8] + x[5] * minus_x9;
	auto poly_0_neg_shared = x[2] * x[11] + x[3] * x[10] + x[4] * x[9] + x[5] * x[8];
	auto poly_0_unique = poly_0_pos_shared - poly_0_neg_shared + x[0] * x[6] + x[1] * minus_x7;
	auto poly_1_unique = poly_0_pos_shared + poly_0_neg_shared + x[0] * x[7] + x[1] * x[6];
	return ring.template batch_reduce_expand(poly_0_unique,poly_1_unique);
}
/* 
Polynomial 0(2) : 1 x0x2 -1 x1x3

Polynomial 1(2) : 1 x0x3 +1 x1x2

Polynomial 2(2) : 1 x0x4 -1 x1x5

Polynomial 3(2) : 1 x0x5 +1 x1x4

Polynomial 4(2) : 1 x0y0 -1 x1y1

Polynomial 5(2) : 1 x0y1 +1 x1y0
 */
template<class Ring>
inline auto fp2_flat_inverse_to_fp6(const std::array<typename Ring::StandardElement, 6> &x, const std::array<typename Ring::StandardElement, 2> &y, Ring &ring) {
	auto minus_x1 = ring.negate(x[1]);
	auto minus_x3 = ring.negate(x[3]);
	auto minus_x5 = ring.negate(x[5]);
	// 0 1
	auto poly_0_unique = x[0] * x[2] + x[1] * minus_x3;
	auto poly_1_unique = x[0] * x[3] + x[1] * x[2];
	// 2 3
	auto poly_2_unique = x[0] * x[4] + x[1] * minus_x5;
	auto poly_3_unique = x[0] * x[5] + x[1] * x[4];
	// 4 5
	auto poly_4_unique = x[0] * x[0] + x[1] * minus_x1;
	auto poly_5_unique = x[0] * x[1] + x[1] * x[0];
	return ring.template batch_reduce_expand(poly_0_unique,poly_1_unique,poly_2_unique,poly_3_unique,poly_4_unique,poly_5_unique);
}
/* 
Polynomial 0(2) : 1 x0x0 +1 x1x1
 */
template<class Ring>
inline auto fp2_flat_inverse_to_fp(const std::array<typename Ring::StandardElement, 2> &x, Ring &ring) {
	auto poly_0_unique = x[0] * x[0] + x[1] * x[1];
	return ring.template batch_reduce_expand(poly_0_unique);
}
/* 
Polynomial 0(1) : 1 x0y0

Polynomial 1(1) : -1 x1y0
 */
template<class Ring>
inline auto fp_flat_inverse_to_fp2(const std::array<typename Ring::StandardElement, 2> &x, const std::array<typename Ring::StandardElement, 1> &y, Ring &ring) {
	auto minus_x0 = ring.negate(x[0]);
	// 0 1
	auto poly_0_unique = x[0] * x[0];
	auto poly_1_unique = x[1] * minus_x0;
	return ring.template batch_reduce_expand(poly_0_unique,poly_1_unique);
}
