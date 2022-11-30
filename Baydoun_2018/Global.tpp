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

	fp_t resultReal = pr_product_difference(aReal, bReal, aImag, bImag);
	fp_t resultImag = pr_product_difference(aReal, bImag, -bReal, aImag);

	return std::complex<fp_t>(resultReal, resultImag);
}

template<typename fp_t>
inline fp_t pr_product_difference(fp_t a, fp_t b, fp_t c, fp_t d) {
	auto tmp = d * c;
	return std::fma(a, b, -tmp) + std::fma(-d, c, tmp);
}
