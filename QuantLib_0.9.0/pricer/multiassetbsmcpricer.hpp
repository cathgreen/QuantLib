#ifndef MULTIASSETBSMCPRICER_HPP
#define MULTIASSETBSMCPRICER_HPP

#include "../market/yieldcurve.hpp"
#include "../market/volatilitytermstructure.hpp"
#include "../products/product.hpp"
#include "../methods/montecarlo/mcparams.hpp"
#include "../methods/montecarlo/pathgenerator.hpp"
#include "../math/stats/statisticscalculator.hpp"

/** Multiasset Monte Carlo pricer in the Black-Scholes model (deterministic rates and vols).
    it is used for single product whose value depends on multiple assets, such as AsianBasketCallPuts
    Current constraint: all assets must be in the same economy, i.e. share the same yield curve.
    All assets must have constant volatilities. Different assets have different spots and volatilities.
*/

class MultiAssetBsMcPricer {
public:
  /** Initializing ctor */
  MultiAssetBsMcPricer(SPtrProduct prod,
                       SPtrYieldCurve discyc,
                       const Vector& vols,
                       const Vector& spots,
                       const Vector& divylds,
                       const Matrix& corrMat,
                       const McParams& mcparams);

  /** Returns the number of variables that can be tracked for stats */
  size_t nVariables() const { return 1; }

  /** Runs the simulation and collects statistics */
  template <typename ITER>
  void simulate(StatisticsCalculator<ITER>& statsCalc, unsigned long npaths);

protected:
  /** Creates and processes one price path.
      It returns the PV of the product
  */
  double processOnePath(Matrix& pricePath);

private:
  // state
  SPtrProduct prod_;               // pointer to product
  SPtrYieldCurve discyc_;          // pointer to yield curve
  Vector vols_;                    // volatilities for all assets
  Vector spots_;                   // spots for all assets
  Vector divylds_;                 // divylds for all assets
  Matrix corrMat_;                 // correlation matrix between assets

  // state to be precomputed
  Vector discfactors_;   // caches the pre-computed discount factors
  Matrix drifts_;        // caches the pre-computed asset drifts, ntimesteps x nassets
  Matrix stdevs_;        // caches the pre-computed asset stdevs, ntimesteps x nassets

  // the final results
  Vector payamts_;       // scratch array for writing the payments after current simulation
  //  Vector currspots_;     // scratch array with the current spots, one per asset
  
  // state of the MC
  McParams mcparams_;   // the Monte Carlo parameters
  SPtrPathGenerator pathgen_; //  pointer to the path generator
};

///////////////////////////////////////////////////////////////////////////////
// Inline definitions

template <typename ITER>
void MultiAssetBsMcPricer::simulate(StatisticsCalculator<ITER>& statsCalc, unsigned long npaths) {
  ASSERT(statsCalc.nVariables() == 1, "MultiAssetBsMcPricer: statsCalc can only take one variable");
  Matrix pricePath;
  double pv = 0;
  for (unsigned long i = 0; i < npaths; ++i) {
    pv = processOnePath(pricePath);
    statsCalc.addSample(&pv, &pv + 1);
  }
}


#endif // MULTIASSETBSMCPRICER_HPP
