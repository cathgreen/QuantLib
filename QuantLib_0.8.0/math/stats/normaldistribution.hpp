#ifndef NORMALD_HPP
#define NORMALD_HPP

#include "univariatedistribution.hpp"
#include "errorfunction.hpp"
#include "../../exception.hpp"
#include "../../defines.hpp"
#include <cmath>

using namespace std;

class NormalDistribution: public UnivariateDistribution {
public:
    NormalDistribution(double mu = 0.0, double sigma = 1.0):
        mu_(mu), sig_(sigma)
    {
        ASSERT(sig_>0, "error: sigma must be positive");
    }

    double pdf(double x) const override {
        return (M_1_SQRT2PI / sig_) * exp(-0.5 * pow((x-mu_)/sig_, 2));
    }

    double cdf(double x) const override {
        ErrorFunction f;
        return 0.5 * f.erfc(-M_SQRT1_2 * (x - mu_)/sig_);
    }

    double invcdf(double p) const override {
        ASSERT(p>0 and p<1, "error: prob. must be in (0,1)");
        ErrorFunction f;
        return -M_SQRT1_2 * sig_ * f.inverfc(2.0 * p) + mu_;
    }


protected:
    double mu_, sig_;
};


#endif // NORMALD_HPP
