//
//  Global.tpp
//  Baydoun_2018
//
//  Created by Анесян Погос Артурович on 21.11.2022.
//

#include "Global.hpp"
#include <numbers>

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
	auto tmp = -d * c;
	return std::fma(a, b, tmp) - std::fma(d, c, tmp);
}

template<typename fp_t>
inline std::complex<fp_t> complex_pr_product_difference(std::complex<fp_t> a, std::complex<fp_t> b, std::complex<fp_t> c, std::complex<fp_t> d) {
	auto tmp = -d * c;
	return cfma(a, b, tmp) - cfma(d, c, tmp);
}

template<typename fp_t>
std::complex<fp_t> cfma(std::complex<fp_t> a, std::complex<fp_t> b, std::complex<fp_t> c) {
	fp_t ar, ai, br, bi, cr, ci, p11, p1, p21, p2;

	ar = a.real(); ai = a.imag();
	br = b.real(); bi = b.imag();
	cr = c.real(); ci = c.imag();

	p11 = fma(ar, br, cr);
	p1 = fma(-ai, bi, p11);
	p21 = fma(ai, br, ci);
	p2 = fma(ar, bi, p21);

	std::complex<fp_t> num(p1, p2);
	return num;
}

template<typename fp_t>
fp_t argp(std::complex<fp_t> inp) {
	const long double _PI = std::numbers::pi_v<long double>;
	fp_t x = std::real(inp);
	fp_t y = std::imag(inp);
	if (x > 0)
		return std::atan2(y, x);
	else{
		fp_t _pi = y < 0 ? -_PI : _PI;
		return x == 0 ? _pi / static_cast<fp_t>(2) : std::atan2(y, x) + _pi;
	}
}
