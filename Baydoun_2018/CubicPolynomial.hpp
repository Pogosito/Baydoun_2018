//
//  ComplexCubicPolynomial.hpp
//  Baydoun_2018
//
//  Created by Pogos Anesyan on 27.10.2022.
//

#ifndef ComplexCubicPolynomial_hpp
#define ComplexCubicPolynomial_hpp

#include <stdio.h>
#include <complex>

class CubicPolynomial final {

// MARK: - Public

public:

	CubicPolynomial(std::complex<long double> b,
					std::complex<long double> c,
					std::complex<long double> d);

	std::vector<std::complex<long double>> calculateRoots();

// MARK: - Private

private:

	// MARK: - Coefficients

	void calculeteAllCoiffecents(std::complex<long double> b,
								 std::complex<long double> c,
								 std::complex<long double> d);

	// B
	std::complex<long double> b;
	std::complex<long double> squareOfB;
	std::complex<long double> cubeOfB;
	std::complex<long double> fourthDegreeOfB;

	// C
	std::complex<long double> c;
	std::complex<long double> squareOfC;
	std::complex<long double> cubeOfC;

	// D
	std::complex<long double> d;
	std::complex<long double> squareOfD;
	std::complex<long double> cubeOfD;

	// Helpers
	std::complex<long double> multiplicationsOfBAndC;
	std::complex<long double> multiplicationsOfSquaresBAndC;

	std::complex<long double> sqrtOfDelta0;
	std::complex<long double> d_0;
	std::complex<long double> coefficient;
	std::complex<long double> smallDeltaL;
	std::complex<long double> A1;
	std::complex<long double> A2;
	std::complex<long double> R1;
	std::complex<long double> R2;
	std::complex<long double> alphaCoefficient;
	std::complex<long double> alpha1;
	std::complex<long double> alpha2;

	// MARK: - Definition 3.1

	// ∆o
	std::complex<long double> calculateDelta0();

	// d_0
	std::complex<long double> calculateD0();

	// ∆l
	std::complex<long double> calculateDeltaL();

	// MARK: - Definition 3.2

	// δl
	std::complex<long double> calculateSmallDeltaL();

	// A1
	std::complex<long double> calculateA1();

	// A2
	std::complex<long double> calculateA2();

	// MARK: - Theorem 3.3

	// R1 && R2
	std::complex<long double> calculateR(bool isR1);

	// α1
	std::complex<long double> calculateAlpha1();

	// α2
	std::complex<long double> calculateAlpha2();
};

#endif /* ComplexCubicPolynomial_hpp */
