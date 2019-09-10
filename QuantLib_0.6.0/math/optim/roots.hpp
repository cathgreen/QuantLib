#ifndef ROOTS_HPP
#define ROOTS_HPP

#include "../../defines.hpp"
#include "../../exception.hpp"
#include "../matrix.hpp"

#include <algorithm>

/**
   Given a function pointer or functor fx defined on the interval [x1, x2], subdivide the interval 
   into n equally spaced segments, and search for zero crossings of the function.
   nroot will be set to the number of bracketing pairs found.
   If nroot is positive, the arrays xb1[0..nroot - 1] and xb2[0..nroot - 1] will be filled 
   sequentially with any bracketing pairs that are found.
   On input, these vectors may have any size, including zero; they will be resized to nroot.
*/
template <typename T>  // T is the type of function
void zbrak(const T& fx, double x1, double x2, int n,
           Vector& xb1, Vector& xb2, int& nroot) {
    int nb = 20;         // initial guess of the number of bracketing subintervals
    xb1.resize(nb);
    xb2.resize(nb);
    
    nroot = 0;
    double dx = (x2 - x1) / n;  // the size of the subintervals
    double x = x1;
    double fp = fx(x);

    for (int i = 0; i < n; ++i) {
        // evaluate the interval betwen x(i) and x(i+1)
        double fc = fx(x + dx);
        if (fp * fc <= 0) {
            // if a sign change occurs, then record this interval
            xb1(nroot) = x;
            xb2(nroot) = x + dx;
            nroot++;
            if (nroot == int(xb1.size())) {
                // Vector tempvec1(xb1), tempvec2(xb2);
                // resize() preserve the original data, and fill in 0 for other elements
                xb1.resize(2 * xb1.size());   
                xb2.resize(2 * xb2.size());
            }    
        }
        fp = fc;
        x += dx;
    }   
}


/**
   Using the secant method, return the root of a function or functor func thought to lie between
   x1 and x2. The root is refined until its accuracy is xacc.
*/
template <typename T>
double rtsec(const T& func, double x1, double x2, double xacc) {
    const int MAXIT = 30;    // maximum allowed number of iterations
    
    double rts = x1, xl = x2;
    double f = func(x1), fl = func(x2);
    if (abs(fl) < abs(f)) {
        rts = x2;
        xl = x1;
        swap(fl, f);
    }
    for (int i = 0; i < MAXIT; ++i) {
        double dx = (xl - rts) * f / (f - fl);  // increment from rts
        xl = rts;
        fl = f;
        rts += dx;
        f = func(rts);
        if (abs(dx) < xacc or f == 0)
            return rts;
    }
    ASSERT(0, "Maximum number of iterations exceeded in rtsec");
}

#endif // ROOTS_HPP
