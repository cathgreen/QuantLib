#ifndef PDE1DSOLVER_HPP
#define PDE1DSOLVER_HPP

#include "pdebase.hpp"
#include "tridiagonalops1d.hpp"
#include "pderesults.hpp"

#include <memory>

/** The 1-d pde solver class, nAsset_ = nLayers_ = 1 */
class Pde1DSolver: public PdeBase {
public:
  /** 1st Ctor, for regular product */
  Pde1DSolver(SPtrProduct product,
              SPtrYieldCurve discountYieldCurve,
              double spot,
              double divyield,
              SPtrVolatilityTermStructure vol,
              Pde1DResults& results,
              bool storeAllResults = false)
    : PdeBase(product), results_(results), storeAllResults_(storeAllResults)
  {
    ASSERT(product->nAssets() == 1, "Pde1DSolver: can only handle one asset only!")
    nAssets_ = product->nAssets();  // should = 1
    nLayers_ = 1; // currently we solve for only one variable
    spdiscyc_ = discountYieldCurve;
    spots_.push_back(spot);
    spaccrycs_.push_back(discountYieldCurve);
    divyields_.push_back(divyield);
    vols_.push_back(vol);
    quantoFlag_ = false;
  }

  /** 2nd Ctor, for quanto product */
  Pde1DSolver(SPtrProduct product, 
              SPtrYieldCurve discountYieldCurve, 
              SPtrYieldCurve growthYieldCurve,
              double spot,
              double divyield,
              SPtrVolatilityTermStructure assetVol,
              SPtrVolatilityTermStructure fxVol,
              double correl,
              Pde1DResults& results,
              bool storeAllResults = false)
    : PdeBase(product), results_(results), storeAllResults_(storeAllResults)
  {
    ASSERT(product->nAssets() == 1, "Pde1DSolver: can only handle one asset only!")
    nAssets_ = product->nAssets();
    nLayers_ = 1;
    spdiscyc_ = discountYieldCurve;
    spots_.push_back(spot);
    spaccrycs_.push_back(growthYieldCurve);
    divyields_.push_back(divyield);
    vols_.push_back(assetVol);
    fxvols_.push_back(fxVol);
    correl_ = correl;
    quantoFlag_ = true;
  }

  /** Solves backwards from one time step to the previous */
  virtual void solveFromStepToStep(ptrdiff_t step, double DT) override;

  /** Initializes the layers, i.e., preValues, currValues
      Each layer (grid function) corresponds to variable solved on the grid.
  */
  virtual void initValLayers() override;

  /** Evaluates the product at the passed-in time step index */
  virtual void evalProduct(size_t stepIdx) override;

  /** Stores the solver results */
  virtual void storeResults() override;

  /** Discounts the grid functions on the current time step, by applying
      the passed-in one-step discount factor. */
  virtual void discountFromStepToStep(double df) override;
  
protected:
  
  // state
  DeltaOp1D<Vector> deltaOpExplicit_, deltaOpImplicit_;
  GammaOp1D<Vector> gammaOpExplicit_, gammaOpImplicit_;
  TridiagonalOp1D<Vector> opExplicit_, opImplicit_;

  bool storeAllResults_;
  Pde1DResults& results_;

  // Matrix values1, values2;  // row: spot node, col: variable, should have only one variable
  //  Matrix* prevValues, * currValues;  // pointer to values1, values2
  shared_ptr<Matrix> prevValues, currValues;
  
};

#endif // PDE1DSOLVER_HPP
