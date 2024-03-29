#ifndef INTERPOLATION1D_HPP
#define INTERPOLATION1D_HPP

#include "../../defines.hpp"
#include "../../exception.hpp"
#include "../matrix.hpp"

#include <cmath>
#include <algorithm>
#include <functional>


/** Helper function for finding the bracketing indices of a value in an ordered vector */
template <typename ARRAY>
void findIndices(const ARRAY& v, double y, size_t& i1, size_t& i2)
{
  if (y <= v[1]) {
    i1 = 0; i2 = 1;
    return;
  }
  if (y >= v[v.size() - 2]) {
    i1 = v.size() - 2; i2 = v.size() - 1;
    return;
  }

  for (size_t i = 1; i < v.size(); ++i) {
    if (y == v[i]) {
      i1 = i2 = i;
      return;
    }
    if (y < v[i]) {
      i1 = i - 1;
      i2 = i;
      return;
    }
  }
}


/** The linear interpolator class */
template <typename ARRAY>
class LinearInterpolation1D {
public:
  /** Ctor */
  LinearInterpolation1D(const Vector& x, const ARRAY& y)
    : xvals_(x), yvals_(y) {
    ASSERT(xvals_.size() == yvals_.size(),
           "LinearInterpolation1D: unequal vector sizes");
  }

  int size() const { return xvals_.size();}

  const Vector& xValues() const { return xvals_;}

  const ARRAY& yValues() const { return yvals_;}

  double getValue(size_t i) const { return yvals_[i];}

  /** Returns the y values by linearly interpolating between neighboring values
  */
  double getValue(double x) const {
    size_t i1, i2;
    findIndices(xvals_, x, i1, i2);
    double y1, y2;
    y1 = getValue(i1);
    if (i1 == i2)
      return y1;
    y2 = getValue(i2);
    double x1 = xvals_[i1], x2 = xvals_[i2];
    return y1 + (y2 - y1) * (x - x1) / (x2 - x1);
  }

protected:

  const Vector xvals_;
  const ARRAY yvals_;
};


#endif // INTERPOLATION1D_HPP
