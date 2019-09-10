#ifndef ERRORFUNCTION_HPP
#define ERRORFUNCTION_HPP

#include "../../exception.hpp"
#include "../../defines.hpp"

class ErrorFunction {
public:
    double erf(double x) const;
    double erfc(double x) const;
    double inverf(double p) const;
    double inverfc(double p) const;

protected:
    double erfccheb(double z) const;

    // state
    static const int ncof = 28;
    static const double cof[28];
};

///////////////////////////////////////////////////////////////////////////////
// Inline definitions

inline double ErrorFunction::erf(double x) const
{
  // Return erf.x/ for any x.
  if (x >= 0.)
    return 1.0 - erfccheb(x);
  else
    return erfccheb(-x) - 1.0;
}

inline double ErrorFunction::erfc(double x) const
{
  // Return erfc.x/ for any x.
  if (x >= 0.)
    return erfccheb(x);
  else
    return 2.0 - erfccheb(-x);
}

inline double ErrorFunction::inverf(double p) const
{
  return inverfc(1. - p);
}




#endif // ERRORFUNCTION_HP
