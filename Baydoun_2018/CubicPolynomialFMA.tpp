//
//  CubicPolynomialFMA.tpp
//  Baydoun_2018
//
//  Created by Анесян Погос Артурович on 16.11.2022.
//

#include "CubicPolynomialFMA.hpp"
#include <iostream>
#include <vector>
#include "Global.hpp"

template<typename f_pt>
CubicPolynomialFMA<f_pt>::CubicPolynomialFMA(f_pt b, f_pt c, f_pt d) {
	this -> b = b;
	this -> c = c;
	this -> d = d;
	calculateAllCoefficients(b, c, d);
}

template<typename f_pt>
void CubicPolynomialFMA<f_pt>::calculateAllCoefficients(f_pt b, f_pt c, f_pt d) {
	squareOfb = b * b;
	cubeOfb = pow(b, 3);
	biquadrateOfB = pow(b, 4);

	squareOfC = c * c;
	cubeOfC = pow(c, 3);
	biquadrateOfC = pow(c, 4);

	squareOfD = d * d;
	cubeOfd = pow(d, 3);
	const f_pt delta0 = calculateDelta0();

	sqrtOfDelta0 = (delta0 > 0) ? sqrt(delta0) : std::complex<f_pt>(0, sqrt(abs(delta0)));

	d_0 = calculateD0();
	coefficient = std::complex<f_pt>(static_cast<f_pt>(0.0L),
									 sqrt(static_cast<f_pt>(3.0L)) / static_cast<f_pt>(9.0L));

	smallDeltaL = calculateSmallDeltaL();
	A1 = calculateA1();
	A2 = calculateA2();
	R1 = calculateR(true);
	R2 = calculateR(false);
	alphaCoefficient = pow(static_cast<f_pt>(4.0L), static_cast<f_pt>(1.0L) / static_cast<f_pt>(3.0L)) / static_cast<f_pt>(2.0L);
	alpha1 = calculateAlpha1();
	alpha2 = calculateAlpha2();
}

template<typename f_pt>
std::vector<std::complex<f_pt>> CubicPolynomialFMA<f_pt>::calculateRoots() {
	const std::complex<f_pt> ms1 = { static_cast<f_pt>(-1.0L), static_cast<f_pt>(0.0L) };
	const std::complex<f_pt> ms2 = { static_cast<f_pt>(1.0L) / static_cast<f_pt>(2.0L), static_cast<f_pt>(-sqrt(3.0L) / static_cast<f_pt>(2.0L)) };
	const std::complex<f_pt> ms3 = std::conj(ms2);

	const std::vector<std::complex<f_pt>> ms = {ms1, ms2, ms3};
	std::vector<std::complex<f_pt>> result;

	for(std::complex<f_pt> m: ms) {
		std::complex<f_pt> brackets = complex_pr_product_difference(alpha1, R1, -m * alpha2, R2);

		std::complex<f_pt> xm = cfma(m, brackets, -std::complex<f_pt>(b, static_cast<f_pt>(0.0L)) / std::complex<f_pt>(static_cast<f_pt>(3.0L), static_cast<f_pt>(0.0L)));
		result.push_back(xm);
	}

	return result;
}

// MARK: - Definition 3.1

// ∆o
template<typename f_pt>
f_pt CubicPolynomialFMA<f_pt>::calculateDelta0() {
	const f_pt bc18d = pr_product_difference(b, c, static_cast<f_pt>(-18.0L), d);
	const f_pt cubeOfbdCubeOfC = std::fma(cubeOfb, d, cubeOfC);
	const f_pt bc = b * c;
	const f_pt helper = pr_product_difference(bc, bc18d, static_cast<f_pt>(4.0L), cubeOfbdCubeOfC);
	const f_pt result = std::fma(static_cast<f_pt>(-27.0L), squareOfD, helper);
//	std::cout << "∆o fma = " << result << std::endl;
	return result;
}

// d_0
template<typename f_pt>
f_pt CubicPolynomialFMA<f_pt>::calculateD0() {
	const f_pt bcd = std::fma(b, c, -d);
	const f_pt bc2d = pr_product_difference(b, c, static_cast<f_pt>(2.0L), d);

	const f_pt _12cbb = std::fma(static_cast<f_pt>(-12.0L), c, squareOfb);
	const f_pt _2squareOfB_BcD_7c_bc2d = pr_product_difference(static_cast<f_pt>(2.0L) * squareOfb, bcd, static_cast<f_pt>(7.0L) * c, bc2d);
	const f_pt squareOfD_12cbbBiquadrateOfC = std::fma(squareOfD, _12cbb, biquadrateOfC);

	const f_pt result = std::fma(static_cast<f_pt>(2.0L) * b * c, _2squareOfB_BcD_7c_bc2d, squareOfD_12cbbBiquadrateOfC);

//	std::cout << "d0 fma = " << result << std::endl;
	return result;
}

// ∆l
template<typename f_pt>
f_pt CubicPolynomialFMA<f_pt>::calculateDeltaL() {

	const f_pt _165dCubeOfB = std::fma(static_cast<f_pt>(16.5L), d, cubeOfb);
	const f_pt _36dCubeOfC = std::fma(static_cast<f_pt>(36.0L), squareOfD, cubeOfC);

	const f_pt helper1 = std::fma(static_cast<f_pt>(8.0L) * cubeOfb, _165dCubeOfB, _36dCubeOfC);

	const f_pt bc_2d = std::fma(b, c, static_cast<f_pt>(-2.0L) * d);
	const f_pt squareOfD_7CubeOfC = std::fma(static_cast<f_pt>(-7.0L), cubeOfC, squareOfD);

	const f_pt helper2 = pr_product_difference(static_cast<f_pt>(11.0L) * cubeOfC,
											   bc_2d,
											   static_cast<f_pt>(-2.0L) * cubeOfb,
											   squareOfD_7CubeOfC);

	const f_pt helper3 = pr_product_difference(squareOfC, helper1, static_cast<f_pt>(-3.0L) * b, helper2);

	const f_pt biquadrateOfBc_6squareOfD = pr_product_difference(biquadrateOfB, c, static_cast<f_pt>(6.0L), squareOfD);
	const f_pt helper4 = std::fma(static_cast<f_pt>(12.0L) * c, biquadrateOfBc_6squareOfD, squareOfb * squareOfD);
	const f_pt helper5 = std::fma(static_cast<f_pt>(97.0L) * squareOfb, squareOfC, static_cast<f_pt>(9.0L) * squareOfD);

	const f_pt helper6 = pr_product_difference(static_cast<f_pt>(2.0L) * b, helper4, static_cast<f_pt>(-3.0L) * d, helper5);

	const f_pt result = pr_product_difference(static_cast<f_pt>(2.0L) * c, helper3, d, helper6);

	return result;
}

// MARK: - Definition 3.2

// δl
template<typename f_pt>
std::complex<f_pt> CubicPolynomialFMA<f_pt>::calculateSmallDeltaL() {
	const std::complex<f_pt> _bcd(std::fma(-b, c, d), static_cast<f_pt>(0.0L));
	const std::complex<f_pt> firstBracketSqrtOfDelta0 = _bcd * sqrtOfDelta0;

	const f_pt bcd = std::fma(b, c, -d);
	const f_pt _2cubeOfCd = std::fma(static_cast<f_pt>(2.0L), cubeOfC, squareOfD);
	const f_pt secondBracket = std::fma(static_cast<f_pt>(4.0L) * b * c, bcd, _2cubeOfCd);

	const std::complex<f_pt> result = complex_pr_product_difference(firstBracketSqrtOfDelta0,
																	std::complex<f_pt>(secondBracket, 0),
																	-coefficient,
																	std::complex<f_pt>(calculateDeltaL(), 0));

//	std::cout << "δl fma = " << result << std::endl;
	return result;
}

// A1
template<typename f_pt>
f_pt CubicPolynomialFMA<f_pt>::calculateA1() {
	const std::complex<f_pt> firstMember_1 = std::complex<f_pt>(static_cast<f_pt>(0.0L), (static_cast<f_pt>(-2.0L) * sqrt(static_cast<f_pt>(3.0L))) / static_cast<f_pt>(3.0L));
	const f_pt _4SquareBOf_13c = pr_product_difference(static_cast<f_pt>(4.0L), squareOfb, static_cast<f_pt>(13.0L), c);
	const f_pt _2SquareBOf_15c = pr_product_difference(static_cast<f_pt>(2.0L), squareOfb, static_cast<f_pt>(15.0L), c);
	const f_pt bc_4SquareBOf_13c_d_2SquareBOf_15c = pr_product_difference(b * c, _4SquareBOf_13c, d, _2SquareBOf_15c);

	const std::complex<f_pt> firstMember_2 = std::complex<f_pt>(bc_4SquareBOf_13c_d_2SquareBOf_15c, static_cast<f_pt>(0.0L));
	const std::complex<f_pt> secondMember =  static_cast<f_pt>(2.0L) * static_cast<f_pt>(c) * sqrtOfDelta0;
	const std::complex<f_pt> result = cfma(firstMember_1, firstMember_2, secondMember);

//	std::cout << "A1 fma = " << result << std::endl;
	return result.imag();
}

// A2
template<typename f_pt>
f_pt CubicPolynomialFMA<f_pt>::calculateA2() {
	const f_pt bc_d = std::fma(b, c, -d);
	const f_pt _20cubeOfC_squareOfd = std::fma(static_cast<f_pt>(20.0L), cubeOfC, -squareOfD);
	const f_pt helper1 = static_cast<f_pt>(2.0L) * squareOfb * std::fma(static_cast<f_pt>(4.0L) * b * c, bc_d, -_20cubeOfC_squareOfd);

	const f_pt _116bd_23squareOfC = pr_product_difference(static_cast<f_pt>(116.0L) * b, d, static_cast<f_pt>(-23.0L), squareOfC);

	const f_pt _11bc_3d = pr_product_difference(static_cast<f_pt>(11.0L) * b, c, static_cast<f_pt>(3.0L), d);

	const f_pt helper2 = pr_product_difference(static_cast<f_pt>(3.0L) * d, _11bc_3d, static_cast<f_pt>(-7.0L), cubeOfC);
	const f_pt helper3 = pr_product_difference(b * squareOfC, _116bd_23squareOfC, static_cast<f_pt>(3.0L) * d, helper2);

	const f_pt firstHalf = std::fma(b, helper1, helper3);

	const f_pt _8squareOfB_c = std::fma(static_cast<f_pt>(8.0L), squareOfb, c);
	const f_pt _10bc_3d = pr_product_difference(static_cast<f_pt>(10.0L) * b, c, static_cast<f_pt>(3.0L), d);
	const f_pt secondHalf = pr_product_difference(squareOfC, _8squareOfB_c, d, _10bc_3d);

	const f_pt result = cfma(-std::complex<f_pt>(static_cast<f_pt>(0.0L), static_cast<f_pt>(sqrt(3.0L))) * sqrtOfDelta0,
												std::complex<f_pt>(secondHalf, static_cast<f_pt>(0.0L)),
												std::complex<f_pt>(firstHalf, static_cast<f_pt>(0.0L))).real();

//	std::cout << "A2 fma = " << result1 << std::endl;

	return result;
}

// MARK: - Theorem 3.3

// R1 && R2
template<typename f_pt>
std::complex<f_pt> CubicPolynomialFMA<f_pt>::calculateR(bool isR1) {
	const std::complex<f_pt> cb_3d = pr_product_difference(c, b, static_cast<f_pt>(3.0L), d);
	const std::complex<f_pt> _2cubeOfB_9Cb_3d = complex_pr_product_difference(std::complex<f_pt>(static_cast<f_pt>(2.0L), 0),
																			  std::complex<f_pt>(cubeOfb, 0),
																			  std::complex<f_pt>(static_cast<f_pt>(9.0L), 0),
																			  cb_3d);

	const std::complex<f_pt> secondMember = _2cubeOfB_9Cb_3d / static_cast<f_pt>(27.0L);

	const std::complex<f_pt> result = isR1
	? pow(cfma(coefficient, sqrtOfDelta0, secondMember), static_cast<f_pt>(1.0L / 3.0L))
	: pow(cfma(coefficient, sqrtOfDelta0, -secondMember), static_cast<f_pt>(1.0L / 3.0L));
//	std::string r = isR1 ? "R1 fma = " : "R2 fma = ";

//	std::cout << r << result << std::endl;
	return result;
}

// α1
template<typename f_pt>
std::complex<f_pt> CubicPolynomialFMA<f_pt>::calculateAlpha1() {
	const std::complex<f_pt> firstBracket = std::complex<f_pt>(static_cast<f_pt>(0.0L), A1) * pow(smallDeltaL, static_cast<f_pt>(1.0L / 3.0L));
	const std::complex<f_pt> secondBracket = -d_0 * R1;

	const f_pt argsDiff = argp(firstBracket) - argp(secondBracket);

	const std::complex<f_pt> result = alphaCoefficient * exp(std::complex<f_pt>(static_cast<f_pt>(0.0L), argsDiff));

//	std::cout << "α1 fma = " << result << std::endl;
	return result;
}

// α2
template<typename f_pt>
std::complex<f_pt> CubicPolynomialFMA<f_pt>::calculateAlpha2() {
	const std::complex<f_pt> firstBracket = A2 * pow(smallDeltaL, static_cast<f_pt>(2.0L / 3.0L));
	const std::complex<f_pt> secondBracket = d_0 * d_0 * R2;

	const f_pt argsDiff = argp(firstBracket) - argp(secondBracket);

	const std::complex<f_pt> result = alphaCoefficient * exp(std::complex<f_pt>(static_cast<f_pt>(0.0L), argsDiff));

//	std::cout << "α2 fma = " << result << std::endl;
	return result;
}
