#ifndef PIECEWISEPOLYNOMIAL_HPP
#define PIECEWISEPOLYNOMIAL_HPP

#include "../../defines.hpp"
#include "../../exception.hpp"
#include "../matrix.hpp"

#include <cmath>
#include <algorithm>
#include <functional>

class PiecewisePolynomial {
public:
    /** Default ctor */
    PiecewisePolynomial() {}

    /** Ctor from breakpoints; order = 0: constant, order = 1: linear, ... 
        All polynomial coefficients are set to zero. */
    template<typename ITER>
    PiecewisePolynomial(ITER xFirst, ITER xLast, size_t order);

     /** Special ctor from breakpoints and values;
         order = 0 piecewise constant, order = 1 linear continuous */
    template<typename XITER, typename YITER>
    PiecewisePolynomial(XITER xFirst, XITER xLast, YITER yFirst, size_t order);
    
    // Dtor
    virtual ~PiecewisePolynomial() {}

    // Properties

    /** number of breakpoints */
    size_t size() const { return x_.size(); }

    /** order of polynomial pieces */
    size_t order() const { return c_.n_rows - 1; }

    // Access

    /** Read access by index i */
    double breakPoint(size_t i) const {
        ASSERT(i < x_.n_elem, "PiecewisePolynomial: access breakpoint out of range");
        return x_(i);
    }
    /** Write access by index i */
    void setBreakPoint(size_t i, double val) {
        ASSERT(i < x_.n_elem, "PiecewisePolynomial: access breakpoint out of range");
        x_(i) = val;
    }
    /** Read-only access to all breakpoints */
    Vector const& breakPoints() const { return x_; }
    /** Read-write access to all breakpoints */
    template<typename ITER>
    void setBreakPoints(ITER xFirst, ITER xLast) {
        x_.resize(xLast - xFirst);
        copy(xFirst, xLast, x_.begin());
    }

    /** Read access by index i */
    double coefficient(size_t i, size_t j) const {
        ASSERT(i < c_.n_rows and j < c_.n_cols, "PiecewisePolynomial: access coefficient out of range");
        return c_(i, j);
    }
    /** Write access by index i */
    void setCoefficient(size_t i, size_t j, double val) {
        ASSERT(i < c_.n_rows and j < c_.n_cols, "PiecewisePolynomial: access coefficient out of range");
        c_(i, j) = val;
    }
    /** Read-only access to coefficient matrix */
    Matrix const& coefficients() const { return c_; }
    /** Read-write access to coefficients at row i, i.e. all coefficients for the ith order */
    Matrix::row_iterator coeff_begin(size_t i) { return c_.begin_row(i); }
    Matrix::row_iterator coeff_end(size_t i) { return c_.end_row(i); }

    // Evaluation
    
    /** Evaluation y(x) */
    double operator()(double x) const;

    /** Value or derivative at one point x
        k = 0 : y(x)
        k > 0 : k-th left derivative at x
    */
    double eval(double x, size_t k = 0) const;

    /** Evaluate at each x in [xFirst, xLast)
        The results will be written by advancing yFirst;
        It assumes that [yFirst, yFirst + (xLast - xFirst)) is a valid range.
        k = 0 : y(x)
        k > 0 : k-th left derivative at x
    */
    template<typename XITER, typename YITER>
    void eval(XITER xFirst, XITER xLast, YITER yFirst, size_t k = 0) const;

    /** Integrate between a and b */
    double integral(double a, double b) const;

    /** Integrate from xStart to each x in [xFirst, xLast)
        If stepwise = true then integration
        The results will be written by advancing yFirst;
        It assumes that [yFirst, yFirst + (xLast - xFirst)) is a valid range.
    */
    template<typename XITER, typename YITER>
    void integral(double xStart, XITER xFirst, XITER xLast, YITER yFirst, bool stepwise = false) const;

    // Computed assignments

    /** Add a constant value to this */
    PiecewisePolynomial& operator+=(double a) { c_.row(0) += a; return *this; }
    /** Subtract a constant value from this */
    PiecewisePolynomial& operator-=(double a) { c_.row(0) -= a; return *this; }
    /** Multiply this with a constant value */
    PiecewisePolynomial& operator*=(double a) { c_ *= a; return *this; }
    /** Divide this by a constant value */
    PiecewisePolynomial& operator/=(double a) { c_ /= a; return *this; }

    // Polynomial algebra

    /** Add p to this */
    PiecewisePolynomial operator+(const PiecewisePolynomial& p) const;
    /** Multiply p with this */
    PiecewisePolynomial operator*(const PiecewisePolynomial& p) const;

protected:
    void assertBreakpointOrder() const {
        auto it = adjacent_find(x_.cbegin(), x_.cend(), greater_equal<double>());
        ASSERT(it == x_.cend(), "PiecewisePolynomial: breakpoints must be in strict increasing order");
    }

    // Returns the greatest index in the vector x_ such that x_[idx] <= x;
    // It returns -1 if x < x[0] 
    ptrdiff_t index(double x) const {
        return upper_bound(x_.begin(), x_.end(), x) -  x_.begin() - 1;
    }

    // Helper function for computing factorials
    inline size_t factorial(size_t n) const {
        return n == 0 ? 1 : n * factorial(n - 1);
    }

    // Helper function for computing derivatives. It returns k-th derivative of p at x_[xIdx] + h
    // Assumes that k <= order()
    double derivative(size_t idx, double h, size_t k = 0) const;

    // Helper function for computing primitives. It returns the integral of p at x_[xIdx] + h
    // Assumes that xIdx >= 0 and xIdx <= size() - 1
    double primitive(size_t xIdx, double h, size_t k) const;
    
    Vector x_;  // breakpoints
    Matrix c_;  // polynomial coefficients
};


inline double PiecewisePolynomial::derivative(size_t idx, double h, size_t k) const {
    double val(0.0);
    size_t ord = order();
    size_t factnp1 = factorial(k);
    double tmp = double(factnp1);
    for (size_t i = k; i <= ord; ++i) {
        val +=  c_(i, idx) * tmp;
        tmp *= h * (i + 1)/(i + 1 - k);
    }
    
    /**
    size_t factnp1 = factorial(order() + 1);
    for (size_t i = order(); i >= k; --i) {
        factnp1 = factnp1 / (i + 1);
        val = c_(i, idx) * factnp1 + val * h / (i - k + 1);
    }
    */

    
    return val;
}

inline double PiecewisePolynomial::primitive(size_t idx, double h, size_t k) const {
    double val(0.0);
    size_t ord = order();

    if ((idx == 0 and h < 0) or (idx == size()-1 and h > 0)) {
        // the integration range is outside the breakpoint domain
        // then compute c*h^k/k! (flat extrapolation)
        val = c_(0, idx) * pow(abs(h), k) / factorial(k);
    }
    else {
        //  the integration range is inside the breakpoint domain
        size_t factnp1 = factorial(ord);
        double tmp = pow(h, k) / factnp1;
        for (size_t i = 0; i <= ord; ++i) {
            val += tmp * c_(i, idx);
            tmp *= h * (i + 1) / (i + k + 1);
        }
    }
    return val;
}


template<typename ITER>
inline PiecewisePolynomial::PiecewisePolynomial(ITER xFirst, ITER xLast, size_t order):
    x_(xLast - xFirst), c_(order+1, xLast - xFirst, arma::fill::zeros) {
    copy(xFirst, xLast, x_.begin());
    assertBreakpointOrder();
}


template<typename XITER, typename YITER>
inline PiecewisePolynomial::PiecewisePolynomial(XITER xFirst, XITER xLast, YITER yFirst, size_t order):
    PiecewisePolynomial(xFirst, xLast, order) {
    ASSERT(order < 2, "PiecewisePolynomial: only 0th and 1st order polynomials can be constructed from values");
    size_t n = xLast - xFirst;
    for (size_t j = 0; j < n; ++j)   // c(0,j) = x(j)
        c_(0, j) = *(yFirst + j);

    if (order == 1) {
        // c(1,j) = (y(j+1) - y(j)) / (x(j+1) - x(j))
        for (size_t j = 0; j < n; ++j) {
            if (j < n-1)
                c_(1, j) = (*(yFirst + j + 1) - *(yFirst + j)) / (*(xFirst + j + 1) - *(xFirst + j));
            else
                c_(1, j) = c_(1, j-1);
        }
    }
}


template<typename XITER, typename YITER>
inline void PiecewisePolynomial::eval(XITER xFirst, XITER xLast, YITER yFirst, size_t k) const {
    if (xFirst == xLast)
        return;  // nothing to do

    XITER xit = xFirst;
    YITER yit = yFirst;
    for (; xit < xLast; ++xit, ++yit)
        *yit = eval(*xit, k);
}

template<typename XITER, typename YITER>
inline void PiecewisePolynomial::integral(double xStart, XITER xFirst,
                                          XITER xLast, YITER yFirst, bool stepwise) const {
    if (xFirst == xLast)
        return;  // nothing to do

    XITER xit = xFirst;
    YITER yit = yFirst;

    for (; xit < xLast; ++xit, ++yit)
            *yit = integral(xStart, *xit);
    if (stepwise) {
        --yit; --xit;
        for (; xit != xFirst; --xit, --yit)
            *yit -= *(yit - 1);
    }
    return;
}


#endif // PIECEWISEPOLYNOMIAL_HPP
