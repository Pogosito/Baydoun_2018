//
//  ComplexCubicPolynomial.hpp
//  Baydoun_2018
//
//  Created by Pogos Anesyan on 27.10.2022.
//

#ifndef CubicPolynomial_hpp
#define CubicPolynomial_hpp

#include <stdio.h>
#include <complex>

template<class f_pt>
class CubicPolynomial final {

// MARK: - Public

public:

	CubicPolynomial(std::complex<f_pt> b,
					std::complex<f_pt> c,
					std::complex<f_pt> d);

	std::vector<std::complex<f_pt>> calculateRoots();

// MARK: - Private

private:

	// MARK: - Coefficients

	void calculeteAllCoiffecents(std::complex<f_pt> b,
								 std::complex<f_pt> c,
								 std::complex<f_pt> d);

	// B
	std::complex<f_pt> b;
	std::complex<f_pt> squareOfB;
	std::complex<f_pt> cubeOfB;
	std::complex<f_pt> fourthDegreeOfB;

	// C
	std::complex<f_pt> c;
	std::complex<f_pt> squareOfC;
	std::complex<f_pt> cubeOfC;

	// D
	std::complex<f_pt> d;
	std::complex<f_pt> squareOfD;
	std::complex<f_pt> cubeOfD;

	// Helpers
	std::complex<f_pt> multiplicationsOfBAndC;
	std::complex<f_pt> multiplicationsOfSquaresBAndC;

	std::complex<f_pt> sqrtOfDelta0;
	std::complex<f_pt> d_0;
	std::complex<f_pt> coefficient;
	std::complex<f_pt> smallDeltaL;
	std::complex<f_pt> A1;
	std::complex<f_pt> A2;
	std::complex<f_pt> R1;
	std::complex<f_pt> R2;
	std::complex<f_pt> alphaCoefficient;
	std::complex<f_pt> alpha1;
	std::complex<f_pt> alpha2;

	// MARK: - Definition 3.1

	// ∆o
	std::complex<f_pt> calculateDelta0();

	// d_0
	std::complex<f_pt> calculateD0();

	// ∆l
	std::complex<f_pt> calculateDeltaL();

	// MARK: - Definition 3.2

	// δl
	std::complex<f_pt> calculateSmallDeltaL();

	// A1
	std::complex<f_pt> calculateA1();

	// A2
	std::complex<f_pt> calculateA2();

	// MARK: - Theorem 3.3

	// R1 && R2
	std::complex<f_pt> calculateR(bool isR1);

	// α1
	std::complex<f_pt> calculateAlpha1();

	// α2
	std::complex<f_pt> calculateAlpha2();
};

#include "CubicPolynomial.tpp"
#endif
