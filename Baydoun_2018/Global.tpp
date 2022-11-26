//
//  Global.tpp
//  Baydoun_2018
//
//  Created by Анесян Погос Артурович on 21.11.2022.
//

#include "Global.hpp"

template<typename fp_t>
std::complex<fp_t> multiplyComplexNumbersFMA(std::complex<fp_t> a, std::complex<fp_t> b) {
	fp_t aReal = a.real();
	fp_t aImag = a.imag();
	fp_t bReal = b.real();
	fp_t bImag = b.imag();
	return std::complex<fp_t>(std::fma(aReal, bReal, -aImag * bImag), std::fma(aReal, bImag, bReal * aImag));
}
