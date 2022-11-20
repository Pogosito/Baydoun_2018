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

int main(int argc, const char * argv[]) {

	std::complex<long double> complexnumber(-2.0, 0.0);
	std::complex<long double> complexnumber2(1.0, 0.0);
	std::complex<long double> complexnumber3(2.0, 0.0);

	CubicPolynomial complexPolynom = CubicPolynomial(complexnumber, complexnumber2, complexnumber3);

	std::vector<std::complex<long double>> arr = complexPolynom.calculateRoots();
	for (int i = 0; i < arr.size(); ++i) {
		std::cout << arr[i] << std::endl;
	}

	std::cout << "-----------------" << std::endl;

	CubicPolynomialFMA complexPolynomFMA = CubicPolynomialFMA(-2.0L, 1.0L, 2.0L);

	std::vector<std::complex<long double>> arrFMA = complexPolynom.calculateRoots();
	for (int i = 0; i < arrFMA.size(); ++i) {
		std::cout << arrFMA[i] << std::endl;
	}

	return 0;
}
