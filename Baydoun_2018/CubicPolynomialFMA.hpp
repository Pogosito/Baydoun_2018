//
//  CubicPolynomialFMA.hpp
//  Baydoun_2018
//
//  Created by Анесян Погос Артурович on 16.11.2022.
//

#ifndef CubicPolynomialFMA_hpp
#define CubicPolynomialFMA_hpp

#include <stdio.h>
#include <complex>

class CubicPolynomialFMA final {

// MARK: - Public

public:

	CubicPolynomialFMA(long double b,
					   long double c,
					   long double d);

	std::vector<std::complex<long double>> calculateRoots();

// MARK: - Private

private:

	// MARK: - Coefficients

	void calculateAllCoefficients(long double b,
								 long double c,
								 long double d);

	// B
	long double b;
	long double squareOfb;
	long double cubeOfb;

	// C
	long double c;
	long double squareOfC;
	long double cubeOfC;

	// D
	long double d;
	long double squareOfD;
	long double cubeOfd;

	// Helpers
	std::complex<long double> sqrtOfDelta0;
	long double d_0;
	std::complex<long double> coefficient;
	long double smallDeltaL;
	long double A1;
	long double A2;
	std::complex<long double> R1;
	std::complex<long double> R2;
	std::complex<long double> alphaCoefficient;
	std::complex<long double> alpha1;
	std::complex<long double> alpha2;

	// MARK: - Definition 3.1

	// ∆o
	long double calculateDelta0();
	
	// d_0
	long double calculateD0();

	// ∆l
	long double calculateDeltaL();

	// MARK: - Definition 3.2

	// δl
	long double calculateSmallDeltaL();

	// A2
	long double calculateA1();

	// A2
	long double calculateA2();

	// MARK: - Theorem 3.3

	std::complex<long double> calculateR(bool isR1);

	// α1
	std::complex<long double> calculateAlpha1();

	// α2
	std::complex<long double> calculateAlpha2();
};

#endif /* CubicPolynomialFMA_hpp */
