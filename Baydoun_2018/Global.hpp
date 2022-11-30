//
//  Global.hpp
//  Baydoun_2018
//
//  Created by Анесян Погос Артурович on 21.11.2022.
//

#ifndef Global_hpp
#define Global_hpp

#include <complex>

template<typename f_pt>
std::complex<f_pt> multiplyComplexNumbersFMA(std::complex<f_pt> a, std::complex<f_pt> b);

template<typename fp_t>
inline fp_t pr_product_difference(fp_t a, fp_t b, fp_t c, fp_t d);

#include "Global.tpp"
#endif /* Global_hpp */
