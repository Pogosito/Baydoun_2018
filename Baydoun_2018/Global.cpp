//
//  Global.cpp
//  Baydoun_2018
//
//  Created by Анесян Погос Артурович on 21.11.2022.
//

#include "Global.hpp"

std::complex<long double> multyplyComplexNumbersFMA(std::complex<long double> a, std::complex<long double> b) {
	long double aReal = a.real();
	long double aImag = a.imag();
	long double bReal = b.real();
	long double bImag = b.imag();
	return std::complex<long double>(std::fma(aReal, bReal, -aImag * bImag), std::fma(aReal, bImag, bReal * aImag));
}
