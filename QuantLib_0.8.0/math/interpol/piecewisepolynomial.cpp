#include "piecewisepolynomial.hpp"

double PiecewisePolynomial::operator()(double x) const {
    size_t n(size());
    double val;

    if (x < x_(0))
        val = c_(0, 0);
    else if (x >= x_(n-1))
        val = c_(0, n-1);
    else {
        size_t idx = index(x);
        val = derivative(idx, x - x_(idx), 0);
    }
    return val;
}

double PiecewisePolynomial::eval(double x, size_t k) const {
    size_t n(size());
    double val;

    if (x < x_(0))
        val = (k == 0) ? c_(0, 0) : 0.0;
    else if (x >= x_(n-1))
        val = (k == 0) ? c_(0, n-1): 0.0;
    else {
        size_t idx = index(x);
        val = derivative(idx, x - x_(idx), k);
    }
    return val;
}


/** Integrate between a and b */
double PiecewisePolynomial::integral(double a, double b) const {
    int isign = 1;   // the sign of the integral
    if (a == b)
        return 0.0;
    else if (a > b) {  // swap them around, and set isgin = -1, 
        double tmp = a;
        a = b;
        b = tmp;
        isign = -1;
    }
    // now a <= b; find the indices of the breakpoints to the left of a and b
    auto idxBkpt0 = index(a);
    auto idxBkpt1 = index(b);
    if (idxBkpt0 == idxBkpt1) {
        // a, b are between same two breakpoints or outside the bkpt range, just do c_(0,i)*(b-a)
        size_t i = idxBkpt0 < 0 ? 0 : idxBkpt0;  // make sure idxBkpt0 is valid index
        return isign * c_(0, i) * (b - a);
    }
    
    // now idxBkpt0 < idxBkpt1
    double val(0.0);
    // first compute the stub piece between a and the first breakpoint to the right of a
    if (idxBkpt0 == -1) {
        val += c_(0, 0) * (x_(0) - a);
    }
    else {
        val += primitive(idxBkpt0, x_(idxBkpt0 + 1) - x_(idxBkpt0), 1)
            - primitive(idxBkpt0, a - x_(idxBkpt0), 1);
    }
    // iterate over the bkpts in the range [x_(idxBkpt0+1), x_(idxBkpt1) ) and accumulate the sum
    idxBkpt0++;
    for (auto idx = idxBkpt0; idx < idxBkpt1; ++idx)
        val +=  primitive(idx, x_(idx + 1) - x_(idx), 1);
    // finally add the stub piece between x_(idxBkpt1) and b
    val += primitive(idxBkpt1, b - x_(idxBkpt1), 1);

    return val;  
}


PiecewisePolynomial PiecewisePolynomial::operator+(const PiecewisePolynomial& p) const {
    size_t n = size() + p.size();              // the sum has at most n breakpoints
    size_t ord = max(order(), p.order());      // the sum has at most ord order
    Vector bkpts(n);  
    
    // merge the breakpoints
    auto bkend = set_union(x_.begin(), x_.end(),
                           p.breakPoints().begin(), p.breakPoints().end(),
                           bkpts.begin());
    auto bkbeg = bkpts.begin();   
    PiecewisePolynomial psum(bkbeg, bkend, ord);
    size_t nbks = psum.size();                 // number of unique breakpoints

     // polynomial summation; compute derivatives of all orders and add them
    Vector tval(nbks), pval(nbks);
    size_t fct = 1;
    for (size_t i = 0; i <= ord; ++i) {
        eval(bkbeg, bkend, tval.begin(), i);
        p.eval(bkbeg, bkend, pval.begin(), i);
        auto rowbeg = psum.coeff_begin(i);
        for (size_t j = 0; j < nbks; ++rowbeg, ++j)
            *rowbeg = (tval(j) + pval(j))/fct;
        fct *= (i+1);
    }
    return psum;   
}



PiecewisePolynomial PiecewisePolynomial::operator*(const PiecewisePolynomial& p) const {
    size_t n = size() + p.size();              // the sum has at most n breakpoints
    size_t ord = max(order(), p.order());      // the sum has at most ord order
    Vector bkpts(n);
    
    // merge the breakpoints
    auto bkend = set_union(x_.begin(), x_.end(),
                           p.breakPoints().begin(), p.breakPoints().end(),
                           bkpts.begin());
    auto bkbeg = bkpts.begin();
    PiecewisePolynomial pprod(bkbeg, bkend, ord);
    size_t nbks = pprod.size();              // number of unique breakpoints

    // polynomial multiplication: using convolution of coefficients
    Vector tval(nbks), pval(nbks);
    for (size_t i = 0; i <= ord; ++i) {
        Matrix::row_iterator rit = pprod.coeff_begin(i);
        for (size_t k = 0; k <= i; ++k) {
            eval(bkbeg, bkend, tval.begin(), k);
            p.eval(bkbeg, bkend, pval.begin(), i - k);
            size_t fact1 = factorial(i);
            size_t fact2 = factorial(i - k);
            for (size_t j = 0; j < nbks; ++j) {
                double tmp = tval[j] / fact1;
                tmp *= pval[j] / fact2;
                *rit += tmp;
                ++rit;
            }
        }
    }
    return pprod;   
}
