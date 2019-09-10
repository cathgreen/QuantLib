#ifndef UNIVD_HPP
#define UNIVD_HPP

#include "../../defines.hpp"
#include "../../exception.hpp"

/** Abstract class for UnivariateDistribution 
 */

class UnivariateDistribution {
public:
    UnivariateDistribution() { }

    virtual double pdf(double x) const = 0;

    virtual double cdf(double x) const = 0;

    virtual double invcdf(double p) const = 0;
};


#endif // UNIVD_HPP
