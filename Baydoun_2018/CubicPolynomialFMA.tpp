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

	squareOfC = pow(c, 2);
	cubeOfC = pow(c, 3);

	squareOfD = pow(d, 2);
	cubeOfd = pow(d, 3);
	const f_pt delta0 = calculateDelta0();

	sqrtOfDelta0 = (delta0 > 0) ? sqrt(delta0) : std::complex<f_pt>(0, sqrt(abs(delta0)));

	d_0 = calculateD0();
	coefficient = std::complex<f_pt>(0, sqrt(3.0)) / std::complex<f_pt>(9.0, 0);
	smallDeltaL = calculateSmallDeltaL();
	A1 = calculateA1();
	A2 = calculateA2();
	R1 = calculateR(true);
	R2 = calculateR(false);
	alphaCoefficient = pow(4, 1.0 / 3.0) / 2.0;
	alpha1 = calculateAlpha1();
	alpha2 = calculateAlpha2();
}

template<typename f_pt>
std::vector<std::complex<f_pt>> CubicPolynomialFMA<f_pt>::calculateRoots() {
	const std::complex<f_pt> ms1 = { -1, 0 };
	const std::complex<f_pt> ms2 = { 1.0 / 2.0, static_cast<f_pt>(-sqrt(3.0) / static_cast<f_pt>(2.0)) };
	const std::complex<f_pt> ms3 = std::conj(ms2);

	const std::vector<std::complex<f_pt>> ms = {ms1, ms2, ms3};
	std::vector<std::complex<f_pt>> result;

	for(std::complex<f_pt> m: ms) {
		std::complex<f_pt> mAlpha1 = multiplyComplexNumbersFMA(m, alpha1);
		std::complex<f_pt> mAlpha1R1 = multiplyComplexNumbersFMA(mAlpha1, R1);
		std::complex<f_pt> squareOfM = multiplyComplexNumbersFMA(m, m);
		std::complex<f_pt> alpha2R2 = multiplyComplexNumbersFMA(alpha2, R2);
		std::complex<f_pt> squareOfMAlpha2R2 = multiplyComplexNumbersFMA(squareOfM, alpha2R2);

		std::complex<f_pt> xm = mAlpha1R1 + squareOfMAlpha2R2 - std::complex<f_pt>(b, 0) / std::complex<f_pt>(3.0, 0);
		result.push_back(xm);
	}

	return result;
}

// MARK: - Definition 3.1

// ∆o
template<typename f_pt>
f_pt CubicPolynomialFMA<f_pt>::calculateDelta0() {
	const f_pt bc18d = std::fmal(b, c, 18 * d);
	const f_pt cubeOfbdCubeOfC = std::fmal(cubeOfb, d, cubeOfC);
	const f_pt result = std::fmal(b * c, bc18d, - 4 * cubeOfbdCubeOfC) - 27 * squareOfD;

//	std::cout << "∆o fma = " << result << std::endl;
	return result;
}

// d_0
template<typename f_pt>
f_pt CubicPolynomialFMA<f_pt>::calculateD0() {
	const f_pt bcd = std::fmal(b, c, -d);
	const f_pt bc2d = std::fmal(b, c, -2 * d);
	const f_pt _12cbb = std::fmal(-12, c, squareOfb);

	const f_pt result = 2 * b * c * (std::fmal(2 * squareOfb, bcd, -7 * c * bc2d)) + std::fmal(pow(d, 2), _12cbb, pow(squareOfC, 2));

//	std::cout << "d0 fma = " << result << std::endl;
	return result;
}

// ∆l
template<typename f_pt>
f_pt CubicPolynomialFMA<f_pt>::calculateDeltaL() {
	const f_pt dCubeOfb = std::fmal(16.5, d, cubeOfb);
	const f_pt _33squareOfbC = std::fmal(33, squareOfb, c);
	const f_pt _6d11bc = std::fmal(6, d, -11 * b * c);
	const f_pt _72cSquareOfB = std::fmal(72, c, -squareOfb);

	const f_pt firstBracket = std::fmal(8 * cubeOfb, dCubeOfb, squareOfC * _33squareOfbC) + 6 * d * _6d11bc;

	const f_pt secondMember = 12 * pow(squareOfb, 2) * c * std::fmal(-7, cubeOfC, squareOfD);
	const f_pt thirdBracket = std::fmal(24, cubeOfb, 291 * d);
	const f_pt fourthMember = cubeOfd * std::fmal(2 * b, _72cSquareOfB, -27 * d);

	const f_pt result = std::fmal(2 * cubeOfC, firstBracket, secondMember) + std::fmal(-squareOfb * squareOfC * d, thirdBracket, fourthMember);

//	std::cout << "∆l fma = " << result << std::endl;
	return result;
}

// MARK: - Definition 3.2

// δl
template<typename f_pt>
std::complex<f_pt> CubicPolynomialFMA<f_pt>::calculateSmallDeltaL() {
	const std::complex<f_pt> _bcd(std::fmal(-b, c, d), 0);
	const std::complex<f_pt> firstBracketSqrtOfDelta0 = _bcd * sqrtOfDelta0;
	
	const f_pt bcd = std::fmal(b, c, -d);
	const f_pt _2cubeOfCd = std::fmal(2, cubeOfC, squareOfD);
	const f_pt secondBracket = std::fmal(4 * b * c, bcd, _2cubeOfCd);
	const std::complex<f_pt> result = firstBracketSqrtOfDelta0 * secondBracket + coefficient * calculateDeltaL();

//	std::cout << "δl fma = " << result << std::endl;
	return result;
}

// A1
template<typename f_pt>
f_pt CubicPolynomialFMA<f_pt>::calculateA1() {
	const std::complex<f_pt> firstMember_1 = std::complex<f_pt>(0, (-2.0 * static_cast<f_pt>(sqrt(3.0))) / 3.0);
	const std::complex<f_pt> firstMember_2	= std::complex<f_pt>(std::fmal(b * c, std::fmal(4, squareOfb, -13 * c), -d * fmal(2, squareOfb, -15 * c)), 0);
	const std::complex<f_pt> firstMember = firstMember_1 * firstMember_2;
	const std::complex<f_pt> secondMember =  static_cast<f_pt>(2.0) * static_cast<f_pt>(c) * sqrtOfDelta0;
	const std::complex<f_pt> result = firstMember + secondMember;

//	std::cout << "A1 fma = " << result << std::endl;
	return result.imag();
}

// A2
template<typename f_pt>
f_pt CubicPolynomialFMA<f_pt>::calculateA2() {
	const f_pt firstHalf = std::fmal(8 * cubeOfb * squareOfC, std::fmal(-5, c, squareOfb), 2 * cubeOfb * d * fmal(-4 * b, c, d)) + 116 * squareOfb * squareOfC * d;

	const f_pt secondHalf = std::fmal(b * c, std::fmal(23, cubeOfC, -99 * squareOfD), -d * fmal(21, cubeOfC, -27 * squareOfD));
	const f_pt bracket = std::fmal(squareOfC, std::fmal(8, squareOfb, c), d * std::fmal(3, d, -10 * b * c));
	const f_pt tenthMember = (std::complex<f_pt>(0, sqrt(3.0)) * bracket * sqrtOfDelta0).real();

	const f_pt result = firstHalf + secondHalf - tenthMember;

//	std::cout << "A2 fma = " << result << std::endl;
	return result;
}

// MARK: - Theorem 3.3

// R1 && R2
template<typename f_pt>
std::complex<f_pt> CubicPolynomialFMA<f_pt>::calculateR(bool isR1) {
	const std::complex<f_pt> firstMember = coefficient * sqrtOfDelta0;
	const std::complex<f_pt> secondMember = (2 * cubeOfb - 9 * std::fmal(c, b, -3 * d)) / 27.0;

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
	const std::complex<f_pt> firstBracket = std::complex<f_pt>(0, A1) * pow(smallDeltaL, static_cast<f_pt>(1.0 / 3.0));
	const std::complex<f_pt> secondBracket = -d_0 * R1;

	const f_pt arg1 = std::arg(firstBracket);
	const f_pt arg2 = std::arg(secondBracket);
	const f_pt argsDiff = arg1 - arg2;

	const std::complex<f_pt> result = alphaCoefficient * exp(std::complex<f_pt>(0, argsDiff));

//	std::cout << "α1 fma = " << result << std::endl;
	return result;
}

// α2
template<typename f_pt>
std::complex<f_pt> CubicPolynomialFMA<f_pt>::calculateAlpha2() {
	const std::complex<f_pt> firstBracket = A2 * pow(smallDeltaL, static_cast<f_pt>(2.0 / 3.0));
	const std::complex<f_pt> secondBracket = d_0 * d_0 * R2;

	const f_pt arg1 = std::arg(firstBracket);
	const f_pt arg2 = std::arg(secondBracket);

	const f_pt argsDiff = arg1 - arg2;

	const std::complex<f_pt> result = alphaCoefficient * exp(std::complex<f_pt>(0, argsDiff));

//	std::cout << "α2 fma = " << result << std::endl;
	return result;
}
