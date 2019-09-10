#include "simplepricers.hpp"
#include "../math/stats/normaldistribution.hpp"
#include <cmath>

double fwdPrice(double spot, double timeToExp, double intRate, double divYield) {
    ASSERT(spot >= 0.0, "spot must be non-negative");
    ASSERT(timeToExp >= 0.0, "time to expiration must be non-negative");
    ASSERT(intRate >= 0.0, "interest rate must be non-negative");
    ASSERT(divYield >= 0.0, "dividend yield must be non-negative");

    return spot * exp((intRate - divYield) * timeToExp);
}


double digitalOptionBS(int payoffType, double spot, double strike, double timeToExp,
                       double intRate, double divYield, double volatility) {
    ASSERT(payoffType == 1 or payoffType == -1, "payoffType must be 1 or -1");
    ASSERT(strike >= 0.0, "strike must be non-negative");
    ASSERT(volatility >= 0.0, "volatility must be non-negative");

    double phi = payoffType;
    double fwd = fwdPrice(spot, timeToExp, intRate, divYield);
    double sigT = volatility * sqrt(timeToExp);
    double d2 = log(fwd / strike) / sigT - 0.5 * sigT;

    NormalDistribution normal;
    double df = exp(-intRate * timeToExp);
    double price = df * normal.cdf(phi * d2);

    return price;
}



double europeanOptionBS(int payoffType, double spot, double strike, double timeToExp,
                        double intRate, double divYield, double volatility) {
    ASSERT(payoffType == 1 or payoffType == -1, "payoffType must be 1 or -1");
    ASSERT(strike >= 0.0, "strike must be non-negative");
    ASSERT(volatility >= 0.0, "volatility must be non-negative");

    double phi = payoffType;
    double fwd = fwdPrice(spot, timeToExp, intRate, divYield);
    double sigT = volatility * sqrt(timeToExp);
    double d1 = log(fwd / strike) / sigT + 0.5 * sigT;
    double d2 = d1 - sigT;

    NormalDistribution normal;
    double df = exp(-intRate * timeToExp);
    double price = fwd * normal.cdf(phi * d1) - strike * normal.cdf(phi * d2);
    price *= phi * df;

    return price;    
}
