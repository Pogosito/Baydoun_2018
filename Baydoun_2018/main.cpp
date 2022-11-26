//
//  main.cpp
//  Baydoun_2018
//
//  Created by Pogos Anesyan on 27.10.2022.
//

#include <iostream>
#include <vector>
#include "CubicPolynomial.hpp"
#include "CubicPolynomialFMA.hpp"
#include "excerpt.h"

#define MAX_DISTANCE 10e-5

void smallTest() {

	//	Exmp: 1
	//	coief: 0 -0.0160999
	//	coief: 1 0.19128
	//	coief: 2 -0.757523
	//	coief: 3 1

	//	Exmp: 2
	//	coief: 0 -0.0124015
	//	coief: 1 0.160732
	//	coief: 2 -0.694404
	//	coief: 3 1

	//	Exmp: 3
	//	coief: 0 0.473331
	//	coief: 1 1.82207
	//	coief: 2 2.33799
	//	coief: 3 1

	std::vector<std::vector<long double>> coefficients = {
		{-0.0160999, 0.19128, -0.757523},
		{-0.694404, 0.160732, -0.0124015},
		{ 2.33799,  1.82207, 0.473331 },
	};

	for (int i = 0; i < coefficients.size(); ++i) {
		std::cout << std::endl;
		std::cout << " -----------------" << std::endl;
		std::cout << "Example # " << i << std::endl;
		std::cout << " -----------------" << std::endl;
		
		std::complex<long double> complexnumber(coefficients[i][0], 0.0);
		std::complex<long double> complexnumber2(coefficients[i][1], 0.0);
		std::complex<long double> complexnumber3(coefficients[i][2], 0.0);

		CubicPolynomial complexPolynom = CubicPolynomial(complexnumber, complexnumber2, complexnumber3);

		std::vector<std::complex<long double>> arr = complexPolynom.calculateRoots();
		for (int i = 0; i < arr.size(); ++i) {
			std::cout << "root " << arr[i] << std::endl;
		}

		std::cout << " -----------------" << std::endl;

		CubicPolynomialFMA <long double>helper(coefficients[i][0], coefficients[i][1], coefficients[i][2]);

		std::vector<std::complex<long double>> arrFMA = helper.calculateRoots();
		for (int i = 0; i < arrFMA.size(); ++i) {
			std::cout << "root fma " << arrFMA[i] << std::endl;
		}

		std::cout << " -----------------" << std::endl;
	}
}

template<typename fp_t>
auto testPolynomial(unsigned int roots_count) {
	fp_t max_absolute_error, max_relative_error;
	std::vector<fp_t> roots(roots_count), coefficients(roots_count + 1);
	generate_polynomial<fp_t>(roots_count, 0, roots_count, 0, MAX_DISTANCE, -1, 1, roots, coefficients);

	CubicPolynomialFMA <fp_t>helper(coefficients[0], coefficients[1], coefficients[2]);
	std::vector<std::complex<fp_t>> myRoots = helper.calculateRoots();

	std::vector<fp_t> imags = { abs(myRoots[0].imag()), abs(myRoots[1].imag()), abs(myRoots[2].imag()) };
	int indexOfRealRoot = std::distance(std::begin(imags), std::min_element(std::begin(imags), std::end(imags)));

	std::vector<std::complex<fp_t>> helperRoots = { std::complex<fp_t>(myRoots[indexOfRealRoot].real(), 0) };

//	for (int i = 0; i < myRoots.size(); ++i) {
//		std::cout << myRoots[i] << std::endl;
//	}
//	std::cout << "Real root " << myRoots[indexOfRealRoot] << std::endl;
//	std::cout << " -----------------" << std::endl;

	compare_roots_complex<fp_t>(1, roots.size(), helperRoots, roots,
								max_absolute_error, max_relative_error);
	return max_absolute_error;
}

int main(int argc, const char * argv[]) {

	std::vector<long double> absoluteErrors = {};

	for (int i = 0; i < 1'000'000; ++i) {
		absoluteErrors.push_back(testPolynomial<long double>(3));
	}

// 	Experiments:

//	1. max absolute error 0.598953
//	2. max absolute error 0.598954
//  3. max absolute error 0.598953
//  4. max absolute error 0.59895
//  5. max absolute error 2.46281 ??? 
//  6. max absolute error 0.598952
//  7. max absolute error 0.598951
//  8. max absolute error 0.598952
//  9. max absolute error 0.598952
//  10. max absolute error 0.598953

	std::vector<long double>::iterator it = std::max_element(std::begin(absoluteErrors), std::end(absoluteErrors));
	std::cout << "max absolute error " << (*it) << std::endl;

	return 0;
}
