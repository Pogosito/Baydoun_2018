//
//  main.cpp
//  Baydoun_2018
//
//  Created by Pogos Anesyan on 27.10.2022.
//

#include <iostream>
#include <vector>
#include "CubicPolynomialFMA.hpp"
#include "CubicPolynomial.hpp"
#include "excerpt.h"
#include <limits>
#define MAX_DISTANCE 1e-5

using namespace  std;

/*

 CALGO
 ACM

 rosetta.org

 954 алгоритм
 1010 алгоритм

*/

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
	
	// Exmp: 4
	//	coief: 0 0.263884
	//	coief: 1 1.23423
	//	coief: 2 1.92424
	//	coief: 3 1
	
	std::vector<std::vector<long double>> coefficients = {
		{-0.0160999, 0.19128, -0.757523},
		{-0.535955, 0.0957494, -0.00570194},
		{-0.694404, 0.160732, -0.0124015},
		{ 2.33799,  1.82207, 0.473331},
		{ 1.92424 , 1.23423, 0.263884},
		{ 1.00391, 0.335945, 0.0374732}
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
	fp_t max_absolute_error;
	fp_t max_relative_error;

	vector<fp_t> roots(roots_count), coefficients(roots_count + 1);

	generate_polynomial<fp_t>(roots_count, 0, roots_count, 0, 1e-5, -1, 1, roots, coefficients);

	vector<complex<fp_t>> roots_computed = CubicPolynomialFMA<fp_t>(coefficients[2], coefficients[1], coefficients[0]).calculateRoots();

	auto result = compare_roots_complex<fp_t>(roots_computed.size(), roots.size(), roots_computed, roots, max_absolute_error, max_relative_error);

	if (abs(max_relative_error) > 1)
		cout << "Error rel = " << max_relative_error << " for coef " << endl << coefficients[0] <<
		endl << coefficients[1] << endl << coefficients[2] << endl << coefficients[3] << endl << coefficients[4];

	return pair<fp_t,fp_t>(max_absolute_error,max_relative_error);
}

int main(int argc, const char * argv[]) {
	float max_absolute_error = 0;
	float max_relative_error = 0;
	float absolute_error = 0;
	float relative_error = 0;
	int index_max_abs_err = 0;
	int index_max_rel_err = 0;

	std::cout.precision(11);

	for (auto i = 0; i < 1'000'000; ++i) {

		auto anw = testPolynomial<float>(3);
		absolute_error = anw.first ;
		relative_error = anw.second;

		if (absolute_error > max_absolute_error){
			max_absolute_error = absolute_error;
			index_max_abs_err = i;
		}

		if (relative_error > max_relative_error){
			max_relative_error = relative_error;
			index_max_rel_err = i;
		}
	}

	std::cout << "***max_absolute_error = " << max_absolute_error << " in test " << index_max_abs_err << std::endl;
	std::cout << "***max_relative_error = " << max_relative_error << " in test " << index_max_rel_err << std::endl;

	return 0;
}
