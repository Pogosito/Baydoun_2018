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
	coefficient = std::complex<f_pt>(static_cast<f_pt>(0.0), sqrt(static_cast<f_pt>(3.0))) / std::complex<f_pt>(static_cast<f_pt>(9.0), static_cast<f_pt>(0.0));
	smallDeltaL = calculateSmallDeltaL();
	A1 = calculateA1();
	A2 = calculateA2();
	R1 = calculateR(true);
	R2 = calculateR(false);
	alphaCoefficient = pow(static_cast<f_pt>(4.0), static_cast<f_pt>(1.0) / static_cast<f_pt>(3.0)) / static_cast<f_pt>(2.0);
	alpha1 = calculateAlpha1();
	alpha2 = calculateAlpha2();
}

template<typename f_pt>
std::vector<std::complex<f_pt>> CubicPolynomialFMA<f_pt>::calculateRoots() {
	const std::complex<f_pt> ms1 = { static_cast<f_pt>(-1.0), static_cast<f_pt>(0.0) };
	const std::complex<f_pt> ms2 = { static_cast<f_pt>(1.0) / static_cast<f_pt>(2.0), static_cast<f_pt>(-sqrt(3.0) / static_cast<f_pt>(2.0)) };
	const std::complex<f_pt> ms3 = std::conj(ms2);

	const std::vector<std::complex<f_pt>> ms = {ms1, ms2, ms3};
	std::vector<std::complex<f_pt>> result;

	for(std::complex<f_pt> m: ms) {
		std::complex<f_pt> mAlpha1 = multiplyComplexNumbersFMA(m, alpha1);
		std::complex<f_pt> mAlpha1R1 = multiplyComplexNumbersFMA(mAlpha1, R1);
		std::complex<f_pt> squareOfM = multiplyComplexNumbersFMA(m, m);
		std::complex<f_pt> alpha2R2 = multiplyComplexNumbersFMA(alpha2, R2);
		std::complex<f_pt> squareOfMAlpha2R2 = multiplyComplexNumbersFMA(squareOfM, alpha2R2);

		std::complex<f_pt> xm = mAlpha1R1 + squareOfMAlpha2R2 - std::complex<f_pt>(b, static_cast<f_pt>(0.0)) / std::complex<f_pt>(static_cast<f_pt>(3.0), static_cast<f_pt>(0.0));
		result.push_back(xm);
	}

	return result;
}

// MARK: - Definition 3.1

// ∆o
template<typename f_pt>
f_pt CubicPolynomialFMA<f_pt>::calculateDelta0() {
	const f_pt bc18d = pr_product_difference(b, c, static_cast<f_pt>(-18.0), d);
	const f_pt cubeOfbdCubeOfC = std::fma(cubeOfb, d, cubeOfC);
	const f_pt bc = b * c;
	const f_pt result = pr_product_difference(bc, bc18d, static_cast<f_pt>(4.0), cubeOfbdCubeOfC) - static_cast<f_pt>(27.0) * squareOfD;

//	std::cout << "∆o fma = " << result << std::endl;
	return result;
}

// d_0
template<typename f_pt>
f_pt CubicPolynomialFMA<f_pt>::calculateD0() {
	const f_pt bcd = std::fma(b, c, -d);
	const f_pt bc2d = pr_product_difference(b, c, static_cast<f_pt>(2.0), d);

	const f_pt _12cbb = std::fma(static_cast<f_pt>(-12.0), c, squareOfb);
	const f_pt _2squareOfB_BcD_7c_bc2d = pr_product_difference(static_cast<f_pt>(2.0) * squareOfb, bcd, static_cast<f_pt>(7.0) * c, bc2d);
	const f_pt squareOfD_12cbbBiquadrateOfC = std::fma(squareOfD, _12cbb, biquadrateOfC);

	const f_pt result = std::fma(static_cast<f_pt>(2.0) * b * c, _2squareOfB_BcD_7c_bc2d, squareOfD_12cbbBiquadrateOfC);

//	std::cout << "d0 fma = " << result << std::endl;
	return result;
}

// ∆l
template<typename f_pt>
f_pt CubicPolynomialFMA<f_pt>::calculateDeltaL() {
	const f_pt dCubeOfb = std::fma(static_cast<f_pt>(16.5), d, cubeOfb);
	const f_pt _33squareOfbC = std::fma(static_cast<f_pt>(33.0), squareOfb, c);

	const f_pt _6d11bc = pr_product_difference(static_cast<f_pt>(6.0), d, static_cast<f_pt>(11.0) * b, c);
	const f_pt _72cSquareOfB = std::fma(static_cast<f_pt>(72.0), c, -squareOfb);

	const f_pt firstBracket = pr_product_difference(static_cast<f_pt>(8.0) * cubeOfb, dCubeOfb, -squareOfC, _33squareOfbC) + static_cast<f_pt>(6.0) * d * _6d11bc;

	const f_pt secondMember = static_cast<f_pt>(12.0) * biquadrateOfB * c * std::fma(static_cast<f_pt>(-7.0), cubeOfC, squareOfD);

	const f_pt thirdBracket = pr_product_difference(static_cast<f_pt>(24.0), cubeOfb, static_cast<f_pt>(-291.0), d);

	const f_pt fourthMember = cubeOfd * pr_product_difference(static_cast<f_pt>(2.0) * b, _72cSquareOfB, static_cast<f_pt>(27.0), d);

	const f_pt result = std::fma(static_cast<f_pt>(2.0) * cubeOfC, firstBracket, secondMember) + std::fma(-squareOfb * squareOfC * d, thirdBracket, fourthMember);

//	std::cout << "∆l fma = " << result << std::endl;
	return result;
}

// MARK: - Definition 3.2

// δl
template<typename f_pt>
std::complex<f_pt> CubicPolynomialFMA<f_pt>::calculateSmallDeltaL() {
	const std::complex<f_pt> _bcd(std::fma(-b, c, d), static_cast<f_pt>(0.0));
	const std::complex<f_pt> firstBracketSqrtOfDelta0 = _bcd * sqrtOfDelta0;

	const f_pt bcd = std::fma(b, c, -d);
	const f_pt _2cubeOfCd = std::fma(static_cast<f_pt>(2.0), cubeOfC, squareOfD);
	const f_pt secondBracket = std::fma(static_cast<f_pt>(4.0) * b * c, bcd, _2cubeOfCd);

	const std::complex<f_pt> result = firstBracketSqrtOfDelta0 * secondBracket + coefficient * calculateDeltaL();

//	std::cout << "δl fma = " << result << std::endl;
	return result;
}

// A1
template<typename f_pt>
f_pt CubicPolynomialFMA<f_pt>::calculateA1() {
	const std::complex<f_pt> firstMember_1 = std::complex<f_pt>(static_cast<f_pt>(0.0), (static_cast<f_pt>(-2.0) * sqrt(static_cast<f_pt>(3.0))) / static_cast<f_pt>(3.0));
	const f_pt _4SquareBOf_13c = pr_product_difference(static_cast<f_pt>(4.0), squareOfb, static_cast<f_pt>(13.0), c);
	const f_pt _2SquareBOf_15c = pr_product_difference(static_cast<f_pt>(2.0), squareOfb, static_cast<f_pt>(15.0), c);
	const f_pt bc_4SquareBOf_13c_d_2SquareBOf_15c = pr_product_difference(b * c, _4SquareBOf_13c, d, _2SquareBOf_15c);

	const std::complex<f_pt> firstMember_2 = std::complex<f_pt>(bc_4SquareBOf_13c_d_2SquareBOf_15c, static_cast<f_pt>(0.0));
	const std::complex<f_pt> firstMember = firstMember_1 * firstMember_2;
	const std::complex<f_pt> secondMember =  static_cast<f_pt>(2.0) * static_cast<f_pt>(c) * sqrtOfDelta0;
	const std::complex<f_pt> result = firstMember + secondMember;

//	std::cout << "A1 fma = " << result << std::endl;
	return result.imag();
}

// A2
template<typename f_pt>
f_pt CubicPolynomialFMA<f_pt>::calculateA2() {
	const f_pt _5cSquareOfb = std::fma(static_cast<f_pt>(-5.0), c, squareOfb);
	const f_pt _4bcd = pr_product_difference(static_cast<f_pt>(-4.0), b, -c, d);

	const f_pt firstHalf = std::fma(static_cast<f_pt>(8.0) * cubeOfb, _5cSquareOfb * squareOfC, static_cast<f_pt>(2.0) * cubeOfb * d * _4bcd) + static_cast<f_pt>(116.0) * squareOfb * squareOfC * d;

	const f_pt _23cubeOfC_99SquareOfD = pr_product_difference(static_cast<f_pt>(23.0), cubeOfC, static_cast<f_pt>(99.0), squareOfD);
	const f_pt _21CubeOfC_27SquareOfD = pr_product_difference(static_cast<f_pt>(21.0), cubeOfC, static_cast<f_pt>(27.0), squareOfD);

	const f_pt secondHalf = pr_product_difference(b * c, _23cubeOfC_99SquareOfD, d, _21CubeOfC_27SquareOfD);

	const f_pt _8SquareOfbC = std::fma(static_cast<f_pt>(8.0), squareOfb, c);

	const f_pt _3d_10bc = pr_product_difference(static_cast<f_pt>(3.0), d, static_cast<f_pt>(10.0) * b, c);

	const f_pt bracket = pr_product_difference(squareOfC, _8SquareOfbC, -d, _3d_10bc);
	const f_pt tenthMember = (std::complex<f_pt>(0.0, static_cast<f_pt>(sqrt(3.0))) * bracket * sqrtOfDelta0).real();

	const f_pt result = firstHalf + secondHalf - tenthMember;

//	std::cout << "A2 fma = " << result << std::endl;
	return result;
}

// MARK: - Theorem 3.3

// R1 && R2
template<typename f_pt>
std::complex<f_pt> CubicPolynomialFMA<f_pt>::calculateR(bool isR1) {
	const std::complex<f_pt> firstMember = coefficient * sqrtOfDelta0;
	const std::complex<f_pt> cb_3d = pr_product_difference(c, b, static_cast<f_pt>(3.0), d);
	const std::complex<f_pt> _2cubeOfB_9Cb_3d = static_cast<f_pt>(2.0) * cubeOfb - static_cast<f_pt>(9.0) * cb_3d;
	const std::complex<f_pt> secondMember = _2cubeOfB_9Cb_3d / static_cast<f_pt>(27.0);

	const std::complex<f_pt> result = isR1
	? pow(firstMember + secondMember, static_cast<f_pt>(1.0 / 3.0))
	: pow(firstMember - secondMember, static_cast<f_pt>(1.0 / 3.0));
	std::string r = isR1 ? "R1 fma = " : "R2 fma = ";

//	std::cout << r << result << std::endl;
	return result;
}

// α1
template<typename f_pt>
std::complex<f_pt> CubicPolynomialFMA<f_pt>::calculateAlpha1() {
	const std::complex<f_pt> firstBracket = std::complex<f_pt>(static_cast<f_pt>(0.0), A1) * pow(smallDeltaL, static_cast<f_pt>(1.0 / 3.0));
	const std::complex<f_pt> secondBracket = -d_0 * R1;

	const f_pt argsDiff = std::arg(firstBracket) - std::arg(secondBracket);

	const std::complex<f_pt> result = alphaCoefficient * exp(std::complex<f_pt>(static_cast<f_pt>(0.0), argsDiff));

//	std::cout << "α1 fma = " << result << std::endl;
	return result;
}

// α2
template<typename f_pt>
std::complex<f_pt> CubicPolynomialFMA<f_pt>::calculateAlpha2() {
	const std::complex<f_pt> firstBracket = A2 * pow(smallDeltaL, static_cast<f_pt>(2.0 / 3.0));
	const std::complex<f_pt> secondBracket = d_0 * d_0 * R2;

	const f_pt argsDiff = std::arg(firstBracket) - std::arg(secondBracket);

	const std::complex<f_pt> result = alphaCoefficient * exp(std::complex<f_pt>(static_cast<f_pt>(0.0), argsDiff));

//	std::cout << "α2 fma = " << result << std::endl;
	return result;
}
