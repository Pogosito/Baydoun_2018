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
		{-0.694404, 0.160732, -0.0124015},
		{ 2.33799,  1.82207, 0.473331 },
		{1.92424 , 1.23423, 0.263884},
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
	
//		std::vector<std::complex<long double>> arr = complexPolynom.calculateRoots();
//		for (int i = 0; i < arr.size(); ++i) {
//			std::cout << "root " << arr[i] << std::endl;
//		}
//
//		std::cout << " -----------------" << std::endl;
		
		CubicPolynomialFMA <long double>helper(coefficients[i][0], coefficients[i][1], coefficients[i][2]);
		
//		std::vector<std::complex<long double>> arrFMA = helper.calculateRoots();
//		for (int i = 0; i < arrFMA.size(); ++i) {
//			std::cout << "root fma " << arrFMA[i] << std::endl;
//		}
//
//		std::cout << " -----------------" << std::endl;
	}
}

// Flexible suppression imaginary part of complex number
template<typename fp_t>
inline complex<fp_t> epsilonComplex(const complex<fp_t> &x)
{
	return abs(x) * numeric_limits<fp_t>::epsilon() > abs(x.imag()) ? complex<fp_t>(x.real(), 0) : x;
}

// Function to test cubic solution:
// testCount - total count of tests
// dist - maximum distance between roots
template <typename fp_t>
void testCubicAdv(const int testCount, const fp_t dist){

	int P = 3; // power, total number of tests
	fp_t low=-1, high=1; // [low, high], max distance between clustered roots
	fp_t absMaxError, relMaxError; // variables for each test Errors
	int cantFind = 0;
	fp_t maxAbsAllofTest = -1, maxRelAllofTest = -1; // maximum from maxAbsoluteError and maxRelError from all [testCount] tests

	long double absErrors = 0;
	long double relError = 0;

	std::vector<fp_t> coefficients(P+1);
	std::vector<fp_t> trueRoots(P);


	for(int i=0; i < testCount; ++i){
		std::vector<fp_t> realRoots;

		generate_polynomial<fp_t>(P, 0, P, 0, dist,
								  low, high, trueRoots, coefficients);

		CubicPolynomialFMA <fp_t>helper(coefficients[2], coefficients[1], coefficients[0]);
		std::vector<std::complex<fp_t>> myRoots = helper.calculateRoots();

		std::complex <fp_t> r1 = epsilonComplex(myRoots[0]);
		std::complex <fp_t> r2 = epsilonComplex(myRoots[1]);
		std::complex <fp_t> r3 = epsilonComplex(myRoots[2]);

		if (!r1.imag()) { realRoots.push_back(r1.real()); }
		if (!r2.imag()) { realRoots.push_back(r2.real()); }
		if (!r3.imag()) { realRoots.push_back(r3.real()); }

		if (realRoots.size() > 0) {
			auto res = compare_roots<fp_t>(realRoots.size(), 3, realRoots,
										   trueRoots, absMaxError, relMaxError);
			if (res == PR_AT_LEAST_ONE_ROOT_IS_NAN) {
				std::cout << "PR_AT_LEAST_ONE_ROOT_IS_NAN" << std::endl;
				continue;
			}
		}

		if(isinf(absMaxError))
			cantFind += 1;
		else{
			maxAbsAllofTest = absMaxError > maxAbsAllofTest? absMaxError : maxAbsAllofTest;
			absErrors += absMaxError;
			maxRelAllofTest = relMaxError > maxRelAllofTest ? relMaxError : maxRelAllofTest;
			relError += relMaxError;
		}
	}
	std::cout<<"Max distance: "<< dist << std::endl;
	std::cout<<"Total count of tests: "<<testCount<<std::endl;
	std::cout<<"Couldn't find roots: " << cantFind <<" times "<<std::endl;
	std::cout<< "----------------------------------------------------" << std::endl;
	std::cout<<"Mean absMaxError = "<< absErrors / (testCount - cantFind) << std::endl;
	std::cout<<"Max {absMaxError_i | i = 0, ..., 1e6} from all of the tests: "<<maxAbsAllofTest<< std::endl;
	std::cout<< "----------------------------------------------------" << std::endl;
	std::cout<<"Mean RelMaxError = "<< relError / (testCount - cantFind) << std::endl;
	std::cout<<"Max {RelMaxError_i | i = 0, ..., 1e6} all of the tests: "<<maxRelAllofTest << std::endl;
	std::cout<< "----------------------------------------------------" << std::endl;
}

int main(int argc, const char * argv[]) {
	testCubicAdv<float>(1'000'000, 1e-5);
	return 0;
}
