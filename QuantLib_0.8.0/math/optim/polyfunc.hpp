#ifndef POLYFUNC_HPP
#define POLYFUNC_HPP

#include "../../defines.hpp"
#include "../../exception.hpp"
#include "../matrix.hpp"

class Polynomial {
public:
    Polynomial(const Vector& coeffs): coeffs_(coeffs) {
        ASSERT(coeffs_.size() > 0,
               "Polynomial: empty vector of coefficients not allowed");
    }

    double operator() (double x) const;
    
private:
    Vector coeffs_;  
};


////////////////////////////////////////////////////////////////////////////.//
// Inline implementations

inline double Polynomial::operator() (double x) const {
    double val = coeffs_(0);
    size_t n = coeffs_.size();
    double pval = x;  // current pow(x,i)
    for (size_t i = 1; i < n; ++i) {
        val += pval * coeffs_(i);
        pval *= x;
    }
    return val;
}

#endif // POLYFUNC_HPP
