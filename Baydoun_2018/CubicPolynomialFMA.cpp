//
//  CubicPolynomialFMA.cpp
//  Baydoun_2018
//
//  Created by Анесян Погос Артурович on 16.11.2022.
//

#include "CubicPolynomialFMA.hpp"
#include <iostream>
#include <vector>
#include "Global.hpp"

CubicPolynomialFMA::CubicPolynomialFMA(long double b,
									   long double c,
									   long double d) {
	this -> b = b;
	this -> c = c;
	this -> d = d;

	calculateAllCoefficients(b, c, d);
}

void CubicPolynomialFMA::calculateAllCoefficients(long double b,
												  long double c,
												  long double d) {
	squareOfb = pow(b, 2);
	cubeOfb = pow(b, 3);

	squareOfC = pow(c, 2);
	cubeOfC = pow(c, 3);

	squareOfD = pow(d, 2);
	cubeOfd = pow(d, 3);
	const long double delta0 = calculateDelta0();

	sqrtOfDelta0 = (delta0 > 0) ? sqrt(delta0) : std::complex<long double>(0, sqrt(abs(delta0)));

	d_0 = calculateD0();
	coefficient = std::complex<long double>(0, sqrt(3.0L)) / 9.0L;
	smallDeltaL = calculateSmallDeltaL();
	A1 = calculateA1();
	A2 = calculateA2();
	R1 = calculateR(true);
	R2 = calculateR(false);
	alphaCoefficient = pow(4, 1.0L / 3.0L) / 2.0L;
}

std::vector<std::complex<long double>> CubicPolynomialFMA::calculateRoots() {
	const std::complex<long double> ms1 = { -1, 0 };
	const std::complex<long double> ms2 = { 1.0 / 2.0, -sqrt(3) / 2.0 };
	const std::complex<long double> ms3 = { 1.0 / 2.0, sqrt(3) / 2.0 };

	const std::vector<std::complex<long double>> ms = {ms1, ms2, ms3};
	std::vector<std::complex<long double>> result;

	for(std::complex<long double> m: ms) {
		std::complex<long double> mAlpha1 = multyplyComplexNumbersFMA(m, alpha1);
		std::complex<long double> mAlpha1R1 = multyplyComplexNumbersFMA(mAlpha1, R1);
		std::complex<long double> squareOfM = multyplyComplexNumbersFMA(m, m);
		std::complex<long double> alpha2R2 = multyplyComplexNumbersFMA(alpha2, R2);
		std::complex<long double> squareOfMAlpha2R2 = multyplyComplexNumbersFMA(squareOfM, alpha2R2);

		std::complex<long double> xm = mAlpha1R1 + squareOfMAlpha2R2 - b / 3.0L;
		m == ms3 ? result.push_back(xm.real()) : result.push_back(xm);
	}

	return result;
}

// MARK: - Definition 3.1

// ∆o
long double CubicPolynomialFMA::calculateDelta0() {
	const long double bc18d = std::fmal(b, c, 18 * d);
	const long double cubeOfbdCubeOfC = std::fmal(cubeOfb, d, cubeOfC);
	const long double result = std::fmal(b * c, bc18d, - 4 * cubeOfbdCubeOfC) - 27 * squareOfD;

//	std::cout << "∆o fma = " << result << std::endl;
	return result;
}

// d_0
long double CubicPolynomialFMA::calculateD0() {
	const long double bcd = std::fmal(b, c, -d);
	const long double bc2d = std::fmal(b, c, -2 * d);
	const long double _12cbb = std::fmal(-12, c, squareOfb);

	const long double result = 2 * b * c * (std::fmal(2 * squareOfb, bcd, -7 * c * bc2d)) + std::fmal(pow(d, 2), _12cbb, pow(squareOfC, 2));

//	std::cout << "d0 fma = " << result << std::endl;
	return result;
}

// ∆l
long double CubicPolynomialFMA::calculateDeltaL() {
	const long double dCubeOfb = std::fmal(16.5, d, cubeOfb);
	const long double _33squareOfbC = std::fmal(33, squareOfb, c);
	const long double _6d11bc = std::fmal(6, d, -11 * b * c);
	const long double _72cSquareOfB = std::fmal(72, c, -squareOfb);

	const long double firstBracket = std::fmal(8 * cubeOfb, dCubeOfb, squareOfC * _33squareOfbC) + 6 * d * _6d11bc;

	const long double secondMember = 12 * pow(squareOfb, 2) * c * std::fmal(-7, cubeOfC, squareOfD);
	const long double thirdBracket = std::fmal(24, cubeOfb, 291 * d);
	const long double fourthMember = cubeOfd * std::fmal(2 * b, _72cSquareOfB, -27 * d);

	long double result = std::fmal(2 * cubeOfC, firstBracket, secondMember) + std::fmal(-squareOfb * squareOfC * d, thirdBracket, fourthMember);

//	std::cout << "∆l fma = " << result << std::endl;
	return result;
}

// MARK: - Definition 3.2

// δl
long double CubicPolynomialFMA::calculateSmallDeltaL() {
	const std::complex<long double> _bcd(std::fmal(-b, c, d), 0);
	const std::complex<long double> firstBracketSqrtOfDelta0 = _bcd * sqrtOfDelta0;
	
	const long double bcd = std::fmal(b, c, -d);
	const long double _2cubeOfCd = std::fmal(2, cubeOfC, squareOfD);
	const long double secondBracket = std::fmal(4 * b * c, bcd, _2cubeOfCd);
	const std::complex<long double> result = firstBracketSqrtOfDelta0 * secondBracket + coefficient * calculateDeltaL();

//	std::cout << "δl fma = " << result << std::endl;
	return result.imag();
}

// A1
long double CubicPolynomialFMA::calculateA1() {
	const std::complex<long double> firstMember = std::complex<long double>(0, (-2.0L * sqrt(3L)) / 3.0L) * std::fmal(b * c, std::fmal(4, squareOfb, -13 * c), -d * fmal(2, squareOfb, -15 * c));
	const std::complex<long double> secondMember = 2.0L * c * sqrtOfDelta0;
	const std::complex<long double> result = firstMember + secondMember;

//	std::cout << "A1 fma = " << result << std::endl;
	return result.imag();
}

// A2
long double CubicPolynomialFMA::calculateA2() {
	const long double firstHalf = std::fmal(8 * cubeOfb * squareOfC, std::fmal(-5, c, squareOfb), 2 * cubeOfb * d * fmal(-4 * b, c, d)) + 116 * squareOfb * squareOfC * d;
	const long double secondHalf = std::fmal(b * c, std::fmal(23, cubeOfC, -99 * squareOfD), -d * fmal(21, cubeOfC, -27 * squareOfD));
	const long double bracket = std::fmal(squareOfC, std::fmal(8, squareOfb, c), d * std::fmal(3, d, -10 * b * c));
	const long double tenthMember = (std::complex<long double>(0, sqrt(3L)) * bracket * sqrtOfDelta0).real();

	const long double result = firstHalf + secondHalf - tenthMember;

//	std::cout << "A2 fma = " << result << std::endl;
	return result;
}

// MARK: - Theorem 3.3

// R1 && R2
std::complex<long double> CubicPolynomialFMA::calculateR(bool isR1) {
	const std::complex<long double> firstMember = coefficient * sqrtOfDelta0;
	const std::complex<long double> secondMember = (2 * cubeOfb - 9 * std::fmal(c, b, -3 * d)) / 27.0L;

	const std::complex<long double> result = isR1 ? pow(firstMember + secondMember, 1.0L / 3.0L) : pow(firstMember - secondMember, 1.0L / 3.0L);
	std::string r = isR1 ? "R1 fma = " : "R2 fma = ";

//	std::cout << r << result << std::endl;
	return result;
}

// α1
std::complex<long double> CubicPolynomialFMA::calculateAlpha1() {
	const std::complex<long double> firstBracket = A1 * pow(smallDeltaL, 1.0L / 3.0L);
	const std::complex<long double> secondBracket = -d_0 * R1;

	const long double arg1 = std::arg(firstBracket);
	const long double arg2 = std::arg(secondBracket);
	const long double argsDiff = arg1 - arg2;

	const std::complex<long double> result = alphaCoefficient * exp(std::complex<long double>(0, argsDiff));

//	std::cout << "α1 fma = " << result << std::endl;
	return result;
}

// α2
std::complex<long double> CubicPolynomialFMA::calculateAlpha2() {
	const std::complex<long double> firstBracket = A2 * pow(smallDeltaL, 2.0L / 3.0L);
	const std::complex<long double> secondBracket = d_0 * d_0 * R2;

	const long double arg1 = std::arg(firstBracket);
	const long double arg2 = std::arg(secondBracket);

	const long double argsDiff = arg1 - arg2;
	const std::complex<long double> result = alphaCoefficient * exp(std::complex<long double>(0, argsDiff));

//	std::cout << "α2 fma = " << result << std::endl;
	return result;
}
