//
//  ComplexCubicPolynomial.cpp
//  Baydoun_2018
//
//  Created by Pogos Anesyan on 27.10.2022.
//

#include "CubicPolynomial.hpp"
#include <iostream>
#include <vector>

CubicPolynomial::CubicPolynomial(std::complex<long double> b,
								 std::complex<long double> c,
								 std::complex<long double> d) {
	this -> b = b;
	this -> c = c;
	this -> d = d;

	calculeteAllCoiffecents(b, c, d);
}

void CubicPolynomial::calculeteAllCoiffecents(std::complex<long double> b,
											  std::complex<long double> c,
											  std::complex<long double> d) {
	squareOfB = b * b;
	cubeOfB = squareOfB * b;
	fourthDegreeOfB = pow(squareOfB, 2);

	squareOfC = pow(c, 2);
	cubeOfC = squareOfC * c;

	squareOfD = pow(d, 2);
	cubeOfD = squareOfD * d;

	multiplicationsOfBAndC = b * c;
	multiplicationsOfSquaresBAndC = squareOfB * squareOfC;

	sqrtOfDelta0 = std::complex<long double>(0, sqrt(calculateDelta0()).imag());
	d_0 = calculateD0();
	coefficient = std::complex<long double>(0, sqrt(3.0L)) / 9.0L;
	smallDeltaL = calculateSmallDeltaL();
	A1 = calculateA1();
	A2 = calculateA2();
	R1 = calculateR(true);
	R2 = calculateR(false);
	alphaCoefficient = pow(4, 1.0L / 3.0L) / 2.0L;
	alpha1 = calculateAlpha1();
	alpha2 = calculateAlpha2();
}

std::vector<std::complex<long double>> CubicPolynomial::calculateRoots() {
	const std::complex<long double> ms1 = { -1, 0 };
	const std::complex<long double> ms2 = { 1.0 / 2.0, -sqrt(3) / 2.0 };
	const std::complex<long double> ms3 = std::conj(ms2);

	const std::vector<std::complex<long double>> ms = {ms1, ms2, ms3};
	std::vector<std::complex<long double>> result;

	for(std::complex<long double> m: ms) {
		std::complex<long double> xm = m * alpha1 * R1 + m * m * alpha2 * R2 - b / 3.0L;
		result.push_back(xm);
	}

	return result;
}

// MARK: - Definition 3.1

// ∆o
std::complex<long double> CubicPolynomial::calculateDelta0() {
	const std::complex<long double> firstMember = -4.0L * cubeOfB * d;
	const std::complex<long double> secondMember = multiplicationsOfSquaresBAndC;
	const std::complex<long double> thirdMember = 18.0L * multiplicationsOfBAndC * d;
	const std::complex<long double> fourthMember = -4.0L * cubeOfC;
	const std::complex<long double> fifthMember = -27.0L * squareOfD;
	const std::complex<long double> result = std::complex<double>((firstMember + secondMember + thirdMember + fourthMember + fifthMember).real(), 0);

//	std::cout << "∆o = " << result << std::endl;
	return result;
}

// d_0
std::complex<long double> CubicPolynomial::calculateD0() {
	const std::complex<long double> firstMember = 4.0L * fourthDegreeOfB * squareOfC;
	const std::complex<long double> secondMember = 4.0L * cubeOfB * c * d;
	const std::complex<long double> thirdMember = 14.0L * squareOfB * cubeOfC;
	const std::complex<long double> fourthMember = squareOfB * squareOfD;
	const std::complex<long double> fifthMember = 28.0L * b * squareOfC * d;
	const std::complex<long double> sixthMember = pow(squareOfC, 2);
	const std::complex<long double> seventhMember = 12.0L * c * squareOfD;

	const std::complex<long double> result = std::complex<double>((firstMember - secondMember - thirdMember + fourthMember + fifthMember + sixthMember - seventhMember).real(), 0);

//	std::cout << "d0 = " << result << std::endl;
	return result;
}

// ∆l
std::complex<long double> CubicPolynomial::calculateDeltaL() {

	const std::complex<long double> firstBracket = (8.0L * cubeOfB * cubeOfB + 132.0L * cubeOfB * d + 36.0L * squareOfD + cubeOfC + 33.0L * multiplicationsOfSquaresBAndC - 66.0L * multiplicationsOfBAndC * d);
	const std::complex<long double> firstMember = 2.0L * cubeOfC * firstBracket;
	const std::complex<long double> secondMember = 12.0L * fourthDegreeOfB * c * (squareOfD - 7.0L * cubeOfC);
	const std::complex<long double> thirdMember = multiplicationsOfSquaresBAndC * d * (24.0L * cubeOfB + 291.0L * d);

	const std::complex<long double> fourthMember = cubeOfD * (144.0L * multiplicationsOfBAndC - 2.0L * cubeOfB - 27.0L * d);
	long double helper = (firstMember + secondMember - thirdMember + fourthMember).real();
	std::complex<long double> result(helper, 0);

//	std::cout << "∆l = " << result << std::endl;
	return result;
}

// MARK: - Definition 3.2

// δl
std::complex<long double> CubicPolynomial::calculateSmallDeltaL() {
	const std::complex<long double> firstBracket = d - multiplicationsOfBAndC;

	const std::complex<long double> secondBracket = 4.0L * multiplicationsOfSquaresBAndC - 4.0L * multiplicationsOfBAndC * d + 2.0L * cubeOfC + squareOfD;
	std::complex<long double> helper(secondBracket.real(), 0);
	const std::complex<long double> firstMember = firstBracket * sqrtOfDelta0 * helper;
	const std::complex<long double> secondMember = coefficient * calculateDeltaL();
	const std::complex<long double> result = firstMember + secondMember;

//	std::cout << "δl = " << result << std::endl;
	return result;
}

// A1
std::complex<long double> CubicPolynomial::calculateA1() {
	const std::complex<long double> firstBracket = 4.0L * cubeOfB * c - 2.0L * d * squareOfB - 13.0L * b * squareOfC + 15.0L * d * c;
	const std::complex<long double> helper(firstBracket.real(), 0);
	const std::complex<long double> firstMember = std::complex<long double>(0, (-2.0L * sqrt(3.0L)) / 3.0L) * helper;
	const std::complex<long double> secondMember = 2.0L * c * sqrtOfDelta0;
	const std::complex<long double> result = firstMember + secondMember;

//	std::cout << "A1 = " << result << std::endl;
	return result;
}

// A2
std::complex<long double> CubicPolynomial::calculateA2() {
	const std::complex<long double> firstMember = 8.0L * multiplicationsOfSquaresBAndC * cubeOfB;
	const std::complex<long double> secondMember = 8.0L * multiplicationsOfBAndC * cubeOfB * d;
	const std::complex<long double> thirdMember = 40.0L * multiplicationsOfSquaresBAndC * multiplicationsOfBAndC;
	const std::complex<long double> fourthMember = 2.0L * cubeOfB * squareOfD;
	const std::complex<long double> fifthMember = 116.0L * multiplicationsOfSquaresBAndC * d;

	const std::complex<long double> helper((firstMember - secondMember - thirdMember + fourthMember + fifthMember).real(), 0);

	const std::complex<long double> sixthMember = 23.0L * multiplicationsOfBAndC * cubeOfC;
	const std::complex<long double> seventhMember = 99.0L * multiplicationsOfBAndC * squareOfD;
	const std::complex<long double> eightMember = 21.0L * cubeOfC * d;
	const std::complex<long double> ninthMember = 27.0L * cubeOfD;

	const std::complex<long double> secondHalf = sixthMember - seventhMember - eightMember + ninthMember;
	
	const std::complex<long double> helper2((8.0L * multiplicationsOfSquaresBAndC - 10.0L * multiplicationsOfBAndC * d + cubeOfC + 3.0L * squareOfD).real(), 0);
	const std::complex<long double> tenthMember = std::complex<long double>(0, sqrt(3.0L)) * helper2 * sqrtOfDelta0;

	const std::complex<long double> result = helper + secondHalf - tenthMember;

	const std::complex<long double> helper3 = std::complex<long double>(result.real(), 0);

//	std::cout << "A2 = " << helper3 << std::endl;
	return helper3;
}

// MARK: - Theorem 3.3

// R1 && R2
std::complex<long double> CubicPolynomial::calculateR(bool isR1) {
	const std::complex<long double> firstMember = coefficient * sqrtOfDelta0;
	const std::complex<long double> secondMember = (2.0L * cubeOfB - 9.0L * c * b + 27.0L * d) / 27.0L;
	const std::complex<long double> helper(secondMember.real(), 0);

	const std::complex<long double> result = isR1 ? pow(firstMember + helper, 1.0L / 3.0L) : pow(firstMember - helper, 1.0L / 3.0L);

	std::string r = isR1 ? "R1 = " : "R2 = ";
//	std::cout << r << result << std::endl;
	return result;
}

// α1
std::complex<long double> CubicPolynomial::calculateAlpha1() {
	const std::complex<long double> firstBracket = A1 * pow(smallDeltaL, 1.0L / 3.0L);
	const std::complex<long double> secondBracket = -d_0 * R1;

	const long double arg1 = std::arg(firstBracket);
	const long double arg2 = std::arg(secondBracket);
	const long double argsDiff = arg1 - arg2;

	const std::complex<long double> result = alphaCoefficient * exp(std::complex<long double>(0, argsDiff));

//	std::cout << "α1 = " << result << std::endl;
	return result;
}

// α2
std::complex<long double> CubicPolynomial::calculateAlpha2() {
	const std::complex<long double> firstBracket = A2 * pow(smallDeltaL, 2.0L / 3.0L);

	const std::complex<long double> secondBracket = d_0 * d_0 * R2;

	const long double arg1 = std::arg(firstBracket);
	const long double arg2 = std::arg(secondBracket);

	const long double argsDiff = arg1 - arg2;
	const std::complex<long double> result = alphaCoefficient * exp(std::complex<long double>(0, argsDiff));

//	std::cout << "α2 = " << result << std::endl;
	return result;
}
