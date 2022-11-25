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

template<class f_pt>
class CubicPolynomialFMA final {

// MARK: - Public

public:

	CubicPolynomialFMA(f_pt b, f_pt c, f_pt d);

	std::vector<std::complex<f_pt>> calculateRoots();

// MARK: - Private

private:

	// MARK: - Coefficients

	void calculateAllCoefficients(f_pt b, f_pt c, f_pt d);

	// B
	f_pt b;
	f_pt squareOfb;
	f_pt cubeOfb;

	// C
	f_pt c;
	f_pt squareOfC;
	f_pt cubeOfC;

	// D
	f_pt d;
	f_pt squareOfD;
	f_pt cubeOfd;

	// Helpers
	std::complex<f_pt> sqrtOfDelta0;
	f_pt d_0;
	std::complex<f_pt> coefficient;
	std::complex<f_pt> smallDeltaL;
	f_pt A1;
	f_pt A2;
	std::complex<f_pt> R1;
	std::complex<f_pt> R2;
	std::complex<f_pt> alphaCoefficient;
	std::complex<f_pt> alpha1;
	std::complex<f_pt> alpha2;

	// MARK: - Definition 3.1

	// ∆o
	f_pt calculateDelta0();
	
	// d_0
	f_pt calculateD0();

	// ∆l
	f_pt calculateDeltaL();

	// MARK: - Definition 3.2

	// δl
	std::complex<f_pt> calculateSmallDeltaL();

	// A2
	f_pt calculateA1();

	// A2
	f_pt calculateA2();

	// MARK: - Theorem 3.3

	std::complex<f_pt> calculateR(bool isR1);

	// α1
	std::complex<f_pt> calculateAlpha1();

	// α2
	std::complex<f_pt> calculateAlpha2();
};

#include "CubicPolynomialFMA.tpp"
#endif /* CubicPolynomialFMA_hpp */
