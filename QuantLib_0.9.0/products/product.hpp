#ifndef PRODUCT_HPP
#define PRODUCT_HPP

#include "../defines.hpp"
#include "../exception.hpp"
#include "../sptr.hpp"
#include "../math/matrix.hpp"

/** The abstract base class for all financial products.
    It must be inherited by specific product payoffs.
*/
class Product {
public:
  /** Default Ctor */
  explicit Product(const string& payccy = "USD")  
    : payccy_(payccy) {}
    
  virtual ~Product() {}

  /** Returns the fixing (observation) times */
  Vector const& fixTimes() const {
    return fixTimes_;
  }

  /** Returns the payment times */
  Vector const& payTimes() const {
    return payTimes_;
  }

  /** Returns the payment amounts */
  Vector const& payAmounts() const {
    return payAmounts_;
  }

  /** The number of assets this product depends on */
  virtual size_t nAssets() const = 0;
    
  /** Evaluates the product given the passed-in path
      The "pricePath" matrix must have nrows = the number of fixing times
      have ncols = the number of products
      write the results onto payAmounts_
  */
  virtual void eval(const Matrix& pricePath) = 0;

  /** Useful for PDE pricing of products with early exercise features
      Evalutes the product at fixing time index idx, for a vector of current spots,
      and a given continuation value.
  */
  virtual void eval(size_t idx, const Vector& spot, double contValue) = 0;

  /** Sets up the time steps, to be used in a numerical method.
      The timesteps are returned in the std::vector<double> timesteps,
      and for each timestep, the corresponding index in the fixingTimes() array
      is in the std::vector<long> stepindex array.
      If stepindex[i] returns -1, this means that timesteps[i] does not correspond
      to a fixing time.
      Note: timesteps are not necessarily = fixtimes, if the DT between two fixtimes
      are too long, then insert more timesteps in between
  */
  virtual void timeSteps(size_t nsteps,
                         vector<double>& timesteps,
                         vector<ptrdiff_t>& stepindex) const;

protected:
  string payccy_;
  Vector fixTimes_;        // the fixing (observation) times
  Vector payTimes_;        // the payment times
    
  /** the payment amounts, same size as payTimes_
      note, one product only have one vector of payAmounts
      if a basket of stock, then add all the payAmounts at one time
  */
  Vector payAmounts_;      
    
};

using SPtrProduct = shared_ptr<Product>;


///////////////////////////////////////////////////////////////////////////////
// Inline definitions

inline
void Product::timeSteps(size_t nsteps,
                        vector<double>& timesteps,
                        vector<ptrdiff_t>& stepindex) const
{
   // first put all the fixing times into a temp array, starting with t = 0
  vector<double> tstemp(1, 0.0);
  // put the indices also in a temp array
  vector<ptrdiff_t> idxtemp(1, -1);
  for (size_t i = 0; i < fixTimes().size(); ++i) {
    tstemp.push_back(fixTimes()[i]);
    idxtemp.push_back(i);
  }
  // if t = 0 is a fixing time, remove the extra t = 0 point
  if (tstemp[0] == tstemp[1]) {
    tstemp.erase(tstemp.begin());
    idxtemp.erase(idxtemp.begin());
  }

  // NOTE: here we can put other time steps such as ex-div dates for discrete divs.

  // compute the timestep size
  double maxTime = tstemp.back();
  double maxdt = maxTime / max(nsteps, size_t(1));  // the maximum dt allowed

  for (size_t i = 0; i < tstemp.size() - 1; ++i) {
    timesteps.push_back(tstemp[i]);           // add the time step
    stepindex.push_back(idxtemp[i]);          // add the index
    double dt = tstemp[i + 1] - tstemp[i];
    if (dt - maxdt > 1.0e-8) {                // insert more steps in between
      size_t n = size_t(dt / maxdt);
      dt /= n;
      double T1 = tstemp[i];
      for (size_t j = 1; j < n; ++j) {
        timesteps.push_back(T1 + j * dt);
        stepindex.push_back(-1);  // this is not a product event
      }
    }
  }
  timesteps.push_back(tstemp.back());
  stepindex.push_back(idxtemp.back());
}

#endif // PRODUCT_HPP
