#ifndef YIELDCURVE_HPP
#define YIELDCURVE_HPP

#include "../defines.hpp"
#include "../exception.hpp"
#include "../sptr.hpp"
#include "../math/interpol/piecewisepolynomial.hpp"

#include <algorithm>

class YieldCurve {
public:

    /** Used to qualify the type of quantities used for building the curve */
    enum class InputType {
        SPOTRATE,
        FWDRATE,
        ZEROBOND
    };

    /** The swap frequency */
    enum class SwapFreq {
        ANNUAL,              // 1/year
        SEMIANNUAL,          // 2/year
        QUARTERLY,           // 4/year
        MONTHLY,             // 12/year
        WEEKLY               // 52/year
    };

    /** The interest rate compounding frequency */
    enum class RateCmpd {
        CONTINUOUS,
        SIMPLE
    };

    /** Ctor from times to Maturity and corresponding continuous compounded rates */
    template<typename XITER, typename YITER>
    YieldCurve(XITER tMatBegin, XITER tMatEnd, YITER rateBegin, YITER rateEnd,
               InputType intype = InputType::SPOTRATE);

    /** Returns the curve currency */
    const string& ccy() const { return ccy_; }
    
    /** Returns the discount factor from observation date to tMat */
    double discount(double tMat) const;

    /** Returns the forward discount factor from tMat1 to tMat2 */
    double fwdDiscount(double tMat1, double tMat2) const;

    /** Returns the spot rate at time tMat */
    double spotRate(double tMat) const;

    /** Returns the forward rate between times tMat1 and tMat2 */
    double fwdRate(double tMat1, double tMat2) const;

    /** Returns the swap rate at time tMat */
    // TODO Not implemented yet, requires frequency arg
    // double swapRate(double tMat1) const;

    /** Returns the forward swap rate times tMat1 and tMat2 */
    // TODO NOt implemented yet, requires frequency arg
    // double fwdSwapRate(double tMat1, double tMat2) const;

private:
    // helper functions
    void initFromZeroBonds();
    void initFromSpotRates();
    void initFromFwdRates();
    
    string ccy_;   // the curve's currency, default constructor set it to "USD"
    PiecewisePolynomial fwdrates_;   // the piecewise constant forward rates
};

using SPtrYieldCurve = shared_ptr<YieldCurve>;

////////////////////////////////////////////////////////////////////////////.//
// Inline implementations


template<typename XITER, typename YITER>
inline YieldCurve::YieldCurve(XITER tMatBegin, XITER tMatEnd, YITER rateBegin, YITER rateEnd,
                              InputType intype)
    : ccy_("USD"), fwdrates_(tMatBegin, tMatEnd, rateBegin, 0)
{
    ASSERT((tMatEnd - tMatBegin) == (rateEnd - rateBegin),
            "YieldCurve: different number of maturities and rates");
    auto it = find_if_not(tMatBegin, tMatEnd, [](double x){return x>= 0.0;});
    ASSERT(it == tMatEnd, "YieldCurve: maturities must be non-negative");
    it = find_if_not(rateBegin, rateEnd, [](double x){return x>= 0.0;});
    ASSERT(it == rateEnd, "YieldCurve: rates must be non-negative");

    switch (intype) {
        // case YieldCurve::InputType::ZEROBOND:
    case InputType::ZEROBOND:
        initFromZeroBonds();
        break;
        // case YieldCurve::InputType::FWDRATE:
    case InputType::FWDRATE:
        initFromFwdRates();
        break;
        //  case YieldCurve::InputType::SPOTRATE:
    case InputType::SPOTRATE:
        initFromSpotRates();
        break;
    default:
        ASSERT(0, "YieldCurve: unknown yield curve input type");   
    }
    
}

#endif // YIELDCURVE_HPP
