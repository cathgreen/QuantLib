#ifndef ASIANBASKETCALLPUT_HPP
#define ASIANBASKETCALLPUT_HPP

#include "product.hpp"

#include <algorithm>
#include <functional>

class AsianBasketCallPut: public Product {
public:
    /** Initializing ctor */
    AsianBasketCallPut(int payoffType,
                       double strike,
                       const Vector& fixingTimes,
                       const Vector& assetQuantities);

    /** The number of assets this product depends on */
    virtual size_t nAssets() const override {
        return assetQuantities_.size();
    }

    virtual void eval(const Matrix& pricePath) override;

private:
    int payoffType_;          // 1: call; -1: put
    double strike_;
    Vector assetQuantities_;  // number of units of each asset in the basket
};

///////////////////////////////////////////////////////////////////////////////
// Inline definitions

inline
AsianBasketCallPut::AsianBasketCallPut(int payoffType,
                                       double strike,
                                       const Vector& fixingTimes,
                                       const Vector& assetQuantities)
    : payoffType_(payoffType), strike_(strike)
{
    ASSERT(payoffType == -1 or payoffType == 1,
           "AsianBasketCallPut: payoffType must be -1 or 1");
    ASSERT(strike >= 0, "AsianBasketCallPut: strike must be non-negative");
    
    // auto it = find_if_not(assetQuantities.begin(), assetQuantities.end(), [](double i){return i>= 0;});
    // ASSERT(it == assetQuantities.end(), "AsianBasketCallPut: asset quantities must be non-negative");
    ASSERT(assetQuantities.size() > 0, "AsianBasketCallPut: asset quantities cannot be empty");
    ASSERT(!any_of(assetQuantities.begin(), assetQuantities.end(), [](double i){return i < 0;}),
            "AsianBasketCallPut: asset quantities must be non-negative");
    
    ASSERT(fixingTimes.size() > 0, "AsianBasketCallPut: fixing times cannot be empty");
    ASSERT(fixingTimes(0) >= 0, "AsianBasketCallPut: the first fixing time must be non-negative");
    auto it = adjacent_find(fixingTimes.begin(), fixingTimes.end(), greater_equal<double>());
                                            // return the first iterator that *(it) >= *(it+1)
    ASSERT(it == fixingTimes.end(),
           "AsianBasketCallPut: the fixing times must be in strict increasing order");

    assetQuantities_ = assetQuantities;
    fixTimes_ = fixingTimes;

    // only the last fixing time is the paytime
    payTimes_.resize(1);
    payTimes_(0) = fixTimes_(fixTimes_.size() - 1);
    payAmounts_.resize(1);
}


inline
void AsianBasketCallPut::eval(const Matrix& pricePath) {
  /**
  size_t nfixings = fixTimes_.size();
  size_t nassets = assetQuantities_.size();
  pricePath.resize(nfixings, nassets);
  */
  
    size_t nfixings = pricePath.n_rows;
    size_t nassets = pricePath.n_cols;
    ASSERT(nfixings == fixTimes_.size(),
           "AsianBasketCallPut: number of rows of pricePath not equal to number of fixing times");
    ASSERT(nassets == assetQuantities_.size(),
           "AsianBasketCallPut: number of cols of pricePath  not equal to number of assets");
    
    double bsktAvg = 0;   // avg price
    for (size_t t = 0; t < nfixings; ++t) {
        double bsktval = 0;
        for (size_t i = 0; i < nassets; ++i)
            bsktval += assetQuantities_(i) * pricePath(t, i);
        bsktAvg += bsktval;
    }
    bsktAvg /= nfixings;

    if (payoffType_ == 1)
        payAmounts_(0) = bsktAvg >= strike_ ? bsktAvg - strike_ : 0.0;
    else
        payAmounts_(0) = bsktAvg >= strike_ ? 0.0 : strike_ - bsktAvg;
}

#endif // ASIANBASKETCALLPUT_HPP
