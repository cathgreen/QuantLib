#ifndef VOLATILITYTERMSTRUCTURE_HPP
#define VOLATILITYTERMSTRUCTURE_HPP

#include "../defines.hpp"
#include "../exception.hpp"
#include "../sptr.hpp"
#include "../math/interpol/piecewisepolynomial.hpp"

#include <algorithm>

class VolatilityTermStructure {
public:
    /** Used to qualify the type of quantities used for building the curve */
    enum class VolType {
        SPOTVOL,
        FWDVOL
    };

    /** Ctor from times to Maturity and corresponding volatility */
    template <typename XITER, typename YITER>
    VolatilityTermStructure(XITER tMatBegin, XITER tMatEnd, YITER volBegin, YITER volEnd,
                            VolType vtype = VolType::SPOTVOL);

    /** Returns the spot vol at time tMat */
    double spotVol(double tMat) const;

    /** Returns the fwd vol between time tMat1 and tMat2 */
    double fwdVol(double tMat1, double tMat2) const;

private:
    // helper functions
    void initFromSpotVol();
    void initFromFwdVol();
    
    PiecewisePolynomial fwdvars_;    // the piecewise constant forward variances
};

using SPtrVolatilityTermStructure = shared_ptr<VolatilityTermStructure>;


////////////////////////////////////////////////////////////////////////////.//
// Inline implementations

template <typename XITER, typename YITER>
inline VolatilityTermStructure::VolatilityTermStructure(XITER tMatBegin, XITER tMatEnd,
                                                 YITER volBegin, YITER volEnd,
                                                 VolType vtype)
    : fwdvars_(tMatBegin, tMatEnd, volBegin, 0)  // create piecewise const (ord = 0), for later change
{
    ASSERT((tMatEnd - tMatBegin) == (volEnd - volBegin),
            "VolatilityTermStructure: different number of maturities and rates");
    auto it = find_if_not(tMatBegin, tMatEnd, [](double x){return x>= 0.0;});
    ASSERT(it == tMatEnd, "VolatilityTermStructure: maturities must be non-negative");
    it = find_if_not(volBegin, volEnd, [](double x){return x>= 0.0;});
    ASSERT(it == volEnd, "VolatilityTermStructure: volatilities must be non-negative");

    switch (vtype) {
    case VolType::SPOTVOL:
        initFromSpotVol();
        break;
    case VolType::FWDVOL:
        initFromFwdVol();
         break;
    default:
        ASSERT(0, "VolatilityTermStructure: unknown volatility term structure input type"); 
    }
}



#endif // VOLATILITYTERMSTRUCTURE_HPP
