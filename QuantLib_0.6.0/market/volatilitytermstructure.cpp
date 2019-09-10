#include "volatilitytermstructure.hpp"

void VolatilityTermStructure::initFromSpotVol() {
    auto vit = fwdvars_.coeff_begin(0);  // read-write access to coef
    size_t n = fwdvars_.size();

    double T1 = fwdvars_.breakPoint(0);  // the first time > 0
    double V1 = *vit;                    // the first volatility
    *vit = V1 * V1;
    double V1squareT1 = V1 * V1 * T1;

    fwdvars_.setBreakPoint(0, 0);  // right coninuous, i.e., var(i) is var between bkp(i) and bkp(i+1)
    ++vit;
    
    for (size_t i = 1; i < n; ++i, ++vit) {
        double T2 = fwdvars_.breakPoint(i);
        double V2 = *vit;
        double V2squareT2 = V2 * V2 * T2;
        double fwdvar = V2squareT2 - V1squareT1;     
        ASSERT(fwdvar >= 0,
               "VolatilityTermStructure: negative variance between T1 = " + to_string(T1) + " and T2 = " + to_string(T2));
        ASSERT( T2 > T1,
               "VolatilityTermStructure: inverse time order between T1 = " + to_string(T1) + " and T2 = " + to_string(T2));
        fwdvar /= (T2 - T1);     // fwd var betwen T1 and T2
        fwdvars_.setBreakPoint(i, T1);   // right coninuous, i.e., var(i) is var between bkp(i) and bkp(i+1)
        *vit = fwdvar;
        T1 = T2;
        V1 = V2;
        V1squareT1 = V2squareT2;
    }
}


void VolatilityTermStructure::initFromFwdVol() {
    auto vit = fwdvars_.coeff_begin(0);
    size_t n = fwdvars_.size();

    double T1 = 0;
    for (size_t i = 0; i < n; ++i, ++vit) {
        double V = *vit;
        *vit = V * V;
        double T2 = fwdvars_.breakPoint(i);
        ASSERT( T2 > T1,
                "VolatilityTermStructure: inverse time order between T1 = " + to_string(T1) + " and T2 = " + to_string(T2));
        fwdvars_.setBreakPoint(i, T1);
        T1 = T2;
    }
}
    

double VolatilityTermStructure::spotVol(double tMat) const {
    ASSERT(tMat >= 0, "VolatilityTermStructure: spot vol for negative times are not allowed");
    if (tMat == 0.0)
        tMat = 1.0e-16;  // handle division by zero
    double svar = fwdvars_.integral(0, tMat);
    double svol = sqrt(svar / tMat);
    return svol;
}
    

double VolatilityTermStructure::fwdVol(double tMat1, double tMat2) const {
    ASSERT(tMat1 >= 0, "VolatilityTermStructure: forward vol for negative times are not allowed");
    ASSERT(tMat2 >= tMat1, "VolatilityTermStructure: maturities for forward vol are out of order");
    if (tMat1 == tMat2) 
        return fwdvars_.eval(tMat1, 0);
    double svar = fwdvars_.integral(tMat1, tMat2);
    double svol = sqrt(svar / (tMat2 - tMat1));
    return svol;
}


