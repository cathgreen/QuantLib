#ifndef PDEBASE_HPP
#define PDEBASE_HPP

#include "pdegrid.hpp"
#include "pdeparams.hpp"
#include "../../products/product.hpp"
#include "../../market/yieldcurve.hpp"
#include "../../market/volatilitytermstructure.hpp"

#include <vector>

/**
Abstract base class for all PDE solvers;
It implements the non-virual method solve() that runs the solver.
Various steps inside the solver correspond to pure virtual methods that must be
overriden in the derived classes.
The class also provides implementation for the initGrid() and updateGrid() methods,
which may be overriden in derived classes.
*/

class PdeBase {
public:

  /** Ctor from product and market data */
  PdeBase(SPtrProduct product,
          SPtrYieldCurve discountYieldCurve,
          vector<double> spots,
          vector<SPtrYieldCurve> const & accrualYieldCurves,
          vector<double> divyields,
          vector<SPtrVolatilityTermStructure> vols)
    : spprod_(product), spdiscyc_(discountYieldCurve), spots_(spots),
      spaccrycs_(accrualYieldCurves), divyields_(divyields), vols_(vols)
  {
    nAssets_ = spprod_->nAssets();
    ASSERT(spots_.size() == nAssets_, "PdeBase: missing or redundant spots!");
    ASSERT(spaccrycs_.size() == nAssets_, "PdeBase: missing or redundant accrual yield curves!");
    ASSERT(divyields_.size() == nAssets_, "PdeBase: missing or redundant dividend yields!");
    ASSERT(vols_.size() == nAssets_, "PdeBase: missing or redundant volatilities!");
  }

  /** Dtor */
  virtual ~PdeBase() {}

  /** Returns the number of factors, i.e. number of dynamic variables.
      This is equal to the number of grid spatial (asset) dimensions. 
      Should equal to nAssets_;
   */
  size_t nFactors() const { return gridAxes_.size(); }

  /** The entry point for the solver; this is the method that the client needs to call */
  void solve(const PdeParams& params);

  /** Initializes the grid axes, sets up the nodes and the bounds */
  virtual void initGrid(double T, const PdeParams& params);

  /** Updates the drift and variance coefficients for the current time step */
  virtual void updateGrid(const PdeParams& params,
                          const Matrix& fwdFactors,
                          const Matrix& fwdVols,
                          size_t stepIdx);

  /** Solves backwards from one time step to the previous */
  virtual void solveFromStepToStep(ptrdiff_t step, double DT) = 0;

  /** Initializes the value layers (grid functions) */
  virtual void initValLayers() = 0;

  /** Evaluates the product at the passed-in time step index */
  virtual void evalProduct(size_t stepIdx) = 0;

  /** Stores the solver results */
  virtual void storeResults() = 0;

  /** Discounts the grid functions on the current time step, by applying
      the passed-in one-step discount factor. */
  virtual void discountFromStepToStep(double df) = 0;
  
protected:

  /** Default ctor */
  PdeBase() {}

  /** Ctor from product; inherited classes must set the other market data */
  PdeBase(SPtrProduct product) : spprod_(product) {}

  // state
  size_t nSteps_;       // num. times steps
  size_t nAssets_;      // num. assets to diffuse
  size_t nLayers_;      // num. PDE variables being solved on the same grid
  double theta_;        // variable control implicit and explicit scheme

  SPtrProduct spprod_;                     // the product being priced
  SPtrYieldCurve spdiscyc_;                // the discounting yield curve 
  vector<double> spots_;                   // the asset spots, of size nAssets
  vector<SPtrYieldCurve> spaccrycs_;       // the accrual yield curve (used for forward calculation), of size nAssets
  vector<double> divyields_;               // the dividend yields, of size nAssets
  vector<SPtrVolatilityTermStructure> vols_;  // the volatility term structures, of size nAssets

  vector<GridAxis> gridAxes_;     // the grid axes, of size nAsset
  vector<double> spotAxis_;       // what is this? I think it is sam as spots_    = =b
  vector<double> alignments_;     // one value per axis at which a grid node must pass through
  vector<double> timesteps_;      // the time axis, real time, size nTimeSteps
  vector<ptrdiff_t> stepindex_;   // the index of timesteps_[i] in product's fixing times; of >= 0, product must be evaluated
                                  // need to to initialized by spprod_->timeSteps((params.nTimeSteps, timesteps_, stepindex_)

  vector<SPtrVolatilityTermStructure> fxvols_;  // the volatility term structure for fx
  double correl_;                                    // the correlation between asset and fx
  bool quantoFlag_;                                  // flag for quanto payoff
  
};

#endif // PDEBASE_HPP


