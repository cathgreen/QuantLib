#ifndef SIMPLEPRICERS_HPP
#define SIMPLEPRICERS_HPP

#include "../exception.hpp"
#include "../defines.hpp"

double fwdPrice(double spot, double timeToExp, double intRate, double divYield);

double digitalOptionBS(int payoffType, double spot, double strike, double timeToExp,
                       double intRate, double divYield, double volatility);


double europeanOptionBS(int payoffType, double spot, double strike, double timeToExp,
                       double intRate, double divYield, double volatility);

// double capFloorBS();

#endif // SIMPLEPRICERS_HP
