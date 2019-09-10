#ifndef ANTITHETICPATHGENERATOR_HPP
#define ANTITHETICPATHGENERATOR_HPP

#include "pathgenerator.hpp"
#include "../../math/random/rng.hpp"
#include "../../math/matrix.hpp"

class AntitheticPathGenerator: public PathGenerator {
public:
  /** Dtor: make any SPtrPathGenerator to antithetic */ 
  AntitheticPathGenerator(SPtrPathGenerator pg);
  
  virtual void next(Matrix& pricePath) override;

protected:
  AntitheticPathGenerator() {}        // default ctor
  
  // state
  SPtrPathGenerator innerpathgen_;
  bool dogenerate_;                   // when dogenerate_ is true a new path must be generated
                                      // iteratively change between true and false
                                      // when false, take the negative of the previous pricePath_
  Matrix pricePath_;                  // the cached generated path  
};


///////////////////////////////////////////////////////////////////////////////
// Inline definitions

inline
AntitheticPathGenerator::AntitheticPathGenerator(SPtrPathGenerator pg)
  : PathGenerator(*pg), innerpathgen_(pg), dogenerate_(true) { }

inline
void AntitheticPathGenerator::next(Matrix& pricePath) {
  if (dogenerate_) {
    innerpathgen_->next(pricePath_);
    dogenerate_ = false;
  }
  else {
    pricePath_ *= -1;
    dogenerate_ = true;
  }
  pricePath = pricePath_;
}




#endif // ANTITHETICPATHGENERATOR_HPP
