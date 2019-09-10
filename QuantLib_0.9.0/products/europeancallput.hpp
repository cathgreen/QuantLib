#ifndef EUROPEANCALLPUT_HPP
#define EUROPEANCALLPUT_HPP

#include "product.hpp"

class EuropeanCallPut: public Product {
public:
  EuropeanCallPut(int payoffType, double strike, double timeToExp);

  /** The number of assets this product depends on */
  virtual size_t nAssets() const override {
    return 1;
  };

  /** Evaluates the product given the passed-in path
      The "pricePath" matrix must have nrows = the number of fixing times
      have ncols = the number of products
      write the results onto payAmounts_
  */
  virtual void eval(const Matrix& pricePath) override;

  virtual void eval(size_t idx, const Vector& spots, double contValue) override;
  
protected:
    int payoffType_;     // 1: call; -1: put
    double strike_;
    double timeToExp_;
};

///////////////////////////////////////////////////////////////////////////////
// Inline definitions

inline
EuropeanCallPut::EuropeanCallPut(int payoffType, double strike, double timeToExp)
    : payoffType_(payoffType), strike_(strike), timeToExp_(timeToExp)
{
    ASSERT(payoffType == -1 or payoffType == 1, "EuropeanCallPut: payoffType must be -1 or 1");
    ASSERT(strike >= 0, "EuropeanCallPut: strike must be non-negative");
    ASSERT(timeToExp >= 0, "EuropeanCallPut: time to expiration must be non-negative");

    // only one fixing time, the expiration
    fixTimes_.resize(1);
    fixTimes_(0) = timeToExp_;

    // only one pay time, the expiration
    payTimes_.resize(1);
    payTimes_(0) = timeToExp_;
    payAmounts_.resize(1);
}

inline
void EuropeanCallPut::eval(const Matrix& pricePath) {
    ASSERT(pricePath.n_rows == fixTimes_.size(), "EuropeanCallPut: pricePath and fixTimes size does not match");
    double S_T = pricePath(0, 0);  // final price
    if (payoffType_ == 1)
        payAmounts_(0) = ( S_T >= strike_ ? (S_T - strike_) : 0 );
    else
        payAmounts_(0) = ( S_T >= strike_ ? 0 : (strike_ -  S_T) );
}

inline
void EuropeanCallPut::eval(size_t idx, const Vector& spots, double contValue)
{
  // the continuation value is not used
  ASSERT(idx == 0, "EuropeanCallPut: wrong fixing time index!");
  double S_T = spots[idx];
  if (payoffType_ == 1)
    payAmounts_[idx] = S_T >= strike_ ? S_T - strike_ : 0.0;
  else
    payAmounts_[idx] = S_T >= strike_ ? 0.0 : strike_ - S_T;
}


#endif // EUROPEANCALLPUT_HPP
