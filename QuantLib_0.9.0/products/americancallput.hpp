#ifndef AMERICANCALLPUT_HPP
#define AMERICANCALLPUT_HPP

#include "europeancallput.hpp"

#include <algorithm>

class AmericanCallPut: public EuropeanCallPut {
public:
  /** Ctor */
  AmericanCallPut(int payoffType, double strike, double timeToExp);

  /** Evaluates the product at fixing time index idx
  */
  virtual void eval(size_t idx, const Vector& spots, double contValue) override;
};

///////////////////////////////////////////////////////////////////////////////
// Inline definitions

inline
AmericanCallPut::AmericanCallPut(int payoffType, double strike, double timeToExp)
  :  EuropeanCallPut(payoffType, strike, timeToExp) {
  // set everyday as fixtimes and paytimes
  // count the number of days between 0 and timeToExp
  size_t nfixings = static_cast<size_t>(timeToExp * DAYS_PER_YEAR) + 1;
  ASSERT(nfixings > 0, "AmericanCallPut: the option has expired");
  fixTimes_.resize(nfixings);
  for (size_t i = 0; i < nfixings - 1; ++i)
    fixTimes_[i] = i / DAYS_PER_YEAR;
  fixTimes_[nfixings - 1] = timeToExp_;

  payTimes_ = fixTimes_;
  payAmounts_.resize(payTimes_.size());
}

inline
void AmericanCallPut::eval(size_t idx, const Vector& spots, double contValue) {
  double spot = spots[0];   // American option only has one variable to eval

  if (idx == payAmounts_.size() - 1) { // this is the last index
    double payoff = (spot - strike_) * payoffType_;
    payAmounts_[idx] = max(payoff, 0.0);
  }
  else { // this is not he last index, check the exercise condition
    double intrinsicValue = (spot - strike_) * payoffType_;
    intrinsicValue = max(intrinsicValue, 0.0);
    payAmounts_[idx] = max(intrinsicValue, contValue);
    // zero out the pay amounts after this index, since American only pays once
    for (size_t j = idx + 1; j < payAmounts_.size(); ++j)
      payAmounts_[j] = 0;
  }
}

#endif // AMERICANCALLPUT_HPP
