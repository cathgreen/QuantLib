#include "yieldcurve.hpp"

void YieldCurve::initFromZeroBonds() {
    auto cit = fwdrates_.coeff_begin(0);
    size_t n = fwdrates_.size();

    double T1 = 0;       // the observation time is t = 0
    double p1 = 1;       // zero bond price maturing at T1
    for (size_t i = 0; i < n; ++i, ++cit) {
        double T2 = fwdrates_.breakPoint(i);
        double p2 = *cit;
        ASSERT(p2 <= 1.0 && p2 > 0, "YieldCurve: zero bond prices must in (0,1]");
        
        double fwdrate = log(p1 / p2) / (T2 - T1);
        ASSERT(fwdrate >= 0.0,
               "YieldCurve: negative fwd rate between T1 = " + to_string(T1) + " and T2 = " + to_string(T2));
        
        fwdrates_.setBreakPoint(i, T1); // breakpoint(i) = T1, breakpoint(i+1) = T2
        *cit = fwdrate;                 // coeff(i) = fwdrate from T1 to T2
        p1 = p2;
        T1 = T2;
    }
}


void YieldCurve::initFromSpotRates() {
    auto cit = fwdrates_.coeff_begin(0);
    size_t n = fwdrates_.size();

    double T1 = fwdrates_.breakPoint(0);       // the first time > 0
    double R1 = *cit;                          // the first spot rate, no need to change
    for (size_t i = 1; i < n; ++i, ++cit) {
        double T2 = fwdrates_.breakPoint(i);
        double R2 = *cit;
        double fwdrate = (R2 * T2 - R1 * T1) / (T2 - T1);
        ASSERT(fwdrate >= 0.0,
               "YieldCurve: negative fwd rate between T1 = " + to_string(T1) + " and T2 = " + to_string(T2));
        
        fwdrates_.setBreakPoint(i, T1); // breakpoint(i) = T1, breakpoint(i+1) = T2
        *cit = fwdrate;                 // coeff(i) = fwdrate from T1 to T2
        R1 = R2;
        T1 = T2;
    }
}


void YieldCurve::initFromFwdRates() {
    // just validate the fwd rates
    auto cit = fwdrates_.coeff_begin(0);
    size_t n = fwdrates_.size();

    double T1 = 0;       // the first time > 0                    
    for (size_t i = 0; i < n; ++i, ++cit) {
        double T2 = fwdrates_.breakPoint(i);
        double fwdrate = *cit;
        ASSERT(fwdrate >= 0.0,
               "YieldCurve: negative fwd rate between T1 = " + to_string(T1) + " and T2 = " + to_string(T2));
        fwdrates_.setBreakPoint(i, T1); // breakpoint(i) = T1, breakpoint(i+1) = T2
        *cit = fwdrate;                 // coeff(i) = fwdrate from T1 to T2
        T1 = T2;
    }
}


double YieldCurve::discount(double tMat) const {
    ASSERT(tMat >= 0.0, "YieldCurve: discount factors for negative times not allowed");
    double ldf = - fwdrates_.integral(0.0, tMat);
    return exp(ldf);
}

double YieldCurve::fwdDiscount(double tMat1, double tMat2) const {
    ASSERT(tMat1 >= 0.0, "YieldCurve: discount factors for negative times not allowed");
    ASSERT(tMat2 >= tMat1, "YieldCurve: maturities are out of order");
    double ldf =  - fwdrates_.integral(tMat1, tMat2);
    return exp(ldf);
}

double YieldCurve::spotRate(double tMat) const {
    ASSERT(tMat >= 0.0, "YieldCurve: spot rates for negative times not allowed");
    double srate = fwdrates_.integral(0.0, tMat);
    return srate / tMat;
}

double YieldCurve::fwdRate(double tMat1, double tMat2) const {
    ASSERT(tMat1 >= 0.0, "YieldCurve: forward rates for negative times not allowed");
    ASSERT(tMat2 >= tMat1, "YieldCurve: maturities are out of order");
    double frate = fwdrates_.integral(tMat1, tMat2);
    return frate / (tMat2 - tMat1);
}
