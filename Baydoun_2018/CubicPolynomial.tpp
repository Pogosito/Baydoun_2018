//
//  ComplexCubicPolynomial.cpp
//  Baydoun_2018
//
//  Created by Pogos Anesyan on 27.10.2022.
//

#include "CubicPolynomial.hpp"
#include <iostream>
#include <vector>
#include "Global.hpp"

template<typename f_pt>
CubicPolynomial<f_pt>::CubicPolynomial(std::complex<f_pt> b,
									   std::complex<f_pt> c,
									   std::complex<f_pt> d) {
	this -> b = b;
	this -> c = c;
	this -> d = d;

	calculeteAllCoiffecents(b, c, d);
}

template<typename f_pt>
void CubicPolynomial<f_pt>::calculeteAllCoiffecents(std::complex<f_pt> b,
													std::complex<f_pt> c,
													std::complex<f_pt> d) {
	squareOfB = b * b;
	cubeOfB = squareOfB * b;
	fourthDegreeOfB = pow(b, 4);

	squareOfC = c * c;
	cubeOfC = squareOfC * c;

	squareOfD = d * d;
	cubeOfD = pow(d, 3);

	multiplicationsOfBAndC = b * c;
	multiplicationsOfSquaresBAndC = squareOfB * squareOfC;

	sqrtOfDelta0 = std::complex<f_pt>(0, sqrt(calculateDelta0()).imag());
	d_0 = calculateD0();
	coefficient = std::complex<f_pt>(0, static_cast<f_pt>(sqrt(3.0L)) / static_cast<f_pt>(9.0L));
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
std::vector<std::complex<f_pt>> CubicPolynomial<f_pt>::calculateRoots() {
	const std::complex<f_pt> ms1 = { static_cast<f_pt>(-1.0L), static_cast<f_pt>(0.0L) };
	const std::complex<f_pt> ms2 = { static_cast<f_pt>(1.0L) / static_cast<f_pt>(2.0L), static_cast<f_pt>(-sqrt(3.0L)) / static_cast<f_pt>(2.0L) };
	const std::complex<f_pt> ms3 = std::conj(ms2);

	const std::vector<std::complex<f_pt>> ms = {ms1, ms2, ms3};
	std::vector<std::complex<f_pt>> result;

	for(std::complex<f_pt> m: ms) {
		std::complex<f_pt> xm = m * alpha1 * R1 + m * m * alpha2 * R2 - b / static_cast<f_pt>(3.0L);
		result.push_back(xm);
	}

	return result;
}

// MARK: - Definition 3.1

// ∆o
template<typename f_pt>
std::complex<f_pt> CubicPolynomial<f_pt>::calculateDelta0() {
	const std::complex<f_pt> firstMember = static_cast<f_pt>(-4.0L) * cubeOfB * d;
	const std::complex<f_pt> secondMember = multiplicationsOfSquaresBAndC;
	const std::complex<f_pt> thirdMember = static_cast<f_pt>(18.0L) * multiplicationsOfBAndC * d;
	const std::complex<f_pt> fourthMember = static_cast<f_pt>(-4.0L) * cubeOfC;
	const std::complex<f_pt> fifthMember = static_cast<f_pt>(-27.0L) * squareOfD;
	const std::complex<f_pt> result = std::complex<f_pt>((firstMember + secondMember + thirdMember + fourthMember + fifthMember).real(), 0);

//	std::cout << "∆o = " << result << std::endl;
	return result;
}

// d_0
template<typename f_pt>
std::complex<f_pt> CubicPolynomial<f_pt>::calculateD0() {
	const std::complex<f_pt> firstMember = static_cast<f_pt>(4.0L) * fourthDegreeOfB * squareOfC;
	const std::complex<f_pt> secondMember = static_cast<f_pt>(4.0L) * cubeOfB * c * d;
	const std::complex<f_pt> thirdMember = static_cast<f_pt>(14.0L) * squareOfB * cubeOfC;
	const std::complex<f_pt> fourthMember = squareOfB * squareOfD;
	const std::complex<f_pt> fifthMember = static_cast<f_pt>(28.0L) * b * squareOfC * d;
	const std::complex<f_pt> sixthMember = pow(squareOfC, static_cast<f_pt>(2.0L));
	const std::complex<f_pt> seventhMember = static_cast<f_pt>(12.0L) * c * squareOfD;

	const std::complex<f_pt> result = std::complex<f_pt>((firstMember - secondMember - thirdMember + fourthMember + fifthMember + sixthMember - seventhMember).real(), 0);

//	std::cout << "d0 = " << result << std::endl;
	return result;
}

// ∆l
template<typename f_pt>
std::complex<f_pt> CubicPolynomial<f_pt>::calculateDeltaL() {

	const std::complex<f_pt> firstBracket = (static_cast<f_pt>(8.0L) * cubeOfB * cubeOfB + static_cast<f_pt>(132.0L) * cubeOfB * d + static_cast<f_pt>(36.0L) * squareOfD + cubeOfC + static_cast<f_pt>(33.0L) * multiplicationsOfSquaresBAndC - static_cast<f_pt>(66.0L) * multiplicationsOfBAndC * d);
	const std::complex<f_pt> firstMember = static_cast<f_pt>(2.0L) * cubeOfC * firstBracket;
	const std::complex<f_pt> secondMember = static_cast<f_pt>(12.0L) * fourthDegreeOfB * c * (squareOfD - static_cast<f_pt>(7.0L) * cubeOfC);
	const std::complex<f_pt> thirdMember = multiplicationsOfSquaresBAndC * d * (static_cast<f_pt>(24.0L) * cubeOfB + static_cast<f_pt>(291.0L) * d);

	const std::complex<f_pt> fourthMember = cubeOfD * (static_cast<f_pt>(144.0L) * multiplicationsOfBAndC - static_cast<f_pt>(2.0L) * cubeOfB - static_cast<f_pt>(27.0L) * d);
	f_pt helper = (firstMember + secondMember - thirdMember + fourthMember).real();
	std::complex<f_pt> result(helper, 0);

//	std::cout << "∆l = " << result << std::endl;
	return result;
}

// MARK: - Definition 3.2

// δl
template<typename f_pt>
std::complex<f_pt> CubicPolynomial<f_pt>::calculateSmallDeltaL() {
	const std::complex<f_pt> firstBracket = d - multiplicationsOfBAndC;

	const std::complex<f_pt> secondBracket = static_cast<f_pt>(4.0L) * multiplicationsOfSquaresBAndC - static_cast<f_pt>(4.0L) * multiplicationsOfBAndC * d + static_cast<f_pt>(2.0L) * cubeOfC + squareOfD;
	std::complex<f_pt> helper(secondBracket.real(), 0);
	const std::complex<f_pt> firstMember = firstBracket * sqrtOfDelta0 * helper;
	const std::complex<f_pt> secondMember = coefficient * calculateDeltaL();
	const std::complex<f_pt> result = firstMember + secondMember;

//	std::cout << "δl = " << result << std::endl;
	return result;
}

// A1
template<typename f_pt>
std::complex<f_pt> CubicPolynomial<f_pt>::calculateA1() {
	const std::complex<f_pt> firstBracket = static_cast<f_pt>(4.0L) * cubeOfB * c - static_cast<f_pt>(2.0L) * d * squareOfB - static_cast<f_pt>(13.0L) * b * squareOfC + static_cast<f_pt>(15.0L) * d * c;
	const std::complex<f_pt> helper(firstBracket.real(), 0);
	const std::complex<f_pt> firstMember = std::complex<f_pt>(0, (static_cast<f_pt>(-2.0L) * static_cast<f_pt>(sqrt(3.0L))) / static_cast<f_pt>(3.0L)) * helper;
	const std::complex<f_pt> secondMember = static_cast<f_pt>(2.0L) * c * sqrtOfDelta0;
	const std::complex<f_pt> result = firstMember + secondMember;

//	std::cout << "A1 = " << result << std::endl;
	return result;
}

// A2
template<typename f_pt>
std::complex<f_pt> CubicPolynomial<f_pt>::calculateA2() {
	const std::complex<f_pt> firstMember = static_cast<f_pt>(8.0L) * multiplicationsOfSquaresBAndC * cubeOfB;
	const std::complex<f_pt> secondMember = static_cast<f_pt>(8.0L) * multiplicationsOfBAndC * cubeOfB * d;
	const std::complex<f_pt> thirdMember = static_cast<f_pt>(40.0L) * multiplicationsOfSquaresBAndC * multiplicationsOfBAndC;
	const std::complex<f_pt> fourthMember = static_cast<f_pt>(2.0L) * cubeOfB * squareOfD;
	const std::complex<f_pt> fifthMember = static_cast<f_pt>(116.0L) * multiplicationsOfSquaresBAndC * d;

	const std::complex<f_pt> helper((firstMember - secondMember - thirdMember + fourthMember + fifthMember).real(), 0);

	const std::complex<f_pt> sixthMember = static_cast<f_pt>(23.0L) * multiplicationsOfBAndC * cubeOfC;
	const std::complex<f_pt> seventhMember = static_cast<f_pt>(99.0L) * multiplicationsOfBAndC * squareOfD;
	const std::complex<f_pt> eightMember = static_cast<f_pt>(21.0L) * cubeOfC * d;
	const std::complex<f_pt> ninthMember = static_cast<f_pt>(27.0L) * cubeOfD;

	const std::complex<f_pt> secondHalf = sixthMember - seventhMember - eightMember + ninthMember;
	
	const std::complex<f_pt> helper2((static_cast<f_pt>(8.0L) * multiplicationsOfSquaresBAndC - static_cast<f_pt>(10.0L) * multiplicationsOfBAndC * d + cubeOfC + static_cast<f_pt>(3.0L) * squareOfD).real(), 0);
	const std::complex<f_pt> tenthMember = std::complex<f_pt>(0, static_cast<f_pt>(sqrt(3.0L))) * helper2 * sqrtOfDelta0;

	const std::complex<f_pt> result = helper + secondHalf - tenthMember;

	const std::complex<f_pt> helper3 = std::complex<f_pt>(result.real(), 0);

//	std::cout << "A2 = " << helper3 << std::endl;
	return helper3;
}

// MARK: - Theorem 3.3

// R1 && R2
template<typename f_pt>
std::complex<f_pt> CubicPolynomial<f_pt>::calculateR(bool isR1) {
	const std::complex<f_pt> firstMember = coefficient * sqrtOfDelta0;
	const std::complex<f_pt> secondMember = (static_cast<f_pt>(2.0L) * cubeOfB - static_cast<f_pt>(9.0L) * c * b + static_cast<f_pt>(27.0L) * d) / static_cast<f_pt>(27.0L);
	const std::complex<f_pt> helper(secondMember.real(), 0);

	const std::complex<f_pt> result = isR1
	? pow(firstMember + helper, static_cast<f_pt>(1.0L) / static_cast<f_pt>(3.0L))
	: pow(firstMember - helper, static_cast<f_pt>(1.0L) / static_cast<f_pt>(3.0L));

//	std::string r = isR1 ? "R1 = " : "R2 = ";
//	std::cout << r << result << std::endl;
	return result;
}

// α1
template<typename f_pt>
std::complex<f_pt> CubicPolynomial<f_pt>::calculateAlpha1() {
	const std::complex<f_pt> firstBracket = A1 * pow(smallDeltaL, static_cast<f_pt>(1.0L) / static_cast<f_pt>(3.0L));
	const std::complex<f_pt> secondBracket = -d_0 * R1;

	const f_pt arg1 = std::arg(firstBracket);
	const f_pt arg2 = std::arg(secondBracket);
	const f_pt argsDiff = arg1 - arg2;

	const std::complex<f_pt> result = alphaCoefficient * exp(std::complex<f_pt>(0, argsDiff));

//	std::cout << "α1 = " << result << std::endl;
	return result;
}

// α2
template<typename f_pt>
std::complex<f_pt> CubicPolynomial<f_pt>::calculateAlpha2() {
	const std::complex<f_pt> firstBracket = A2 * pow(smallDeltaL, static_cast<f_pt>(2.0L) / static_cast<f_pt>(3.0L));

	const std::complex<f_pt> secondBracket = d_0 * d_0 * R2;

	const f_pt arg1 = std::arg(firstBracket);
	const f_pt arg2 = std::arg(secondBracket);

	const f_pt argsDiff = arg1 - arg2;
	const std::complex<f_pt> result = alphaCoefficient * exp(std::complex<f_pt>(0, argsDiff));

//	std::cout << "α2 = " << result << std::endl;
	return result;
}
