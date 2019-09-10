#ifndef BSMCPRICER_HPP
#define BSMCPRICER_HPP

#include "../market/yieldcurve.hpp"
#include "../market/volatilitytermstructure.hpp"
#include "../products/product.hpp"
#include "../methods/montecarlo/mcparams.hpp"
#include "../methods/montecarlo/pathgenerator.hpp"
#include "../math/stats/statisticscalculator.hpp"

class BsMcPricer {
public:
    /** pathgen_.nFactors() must equal to prod.nAssets()
        For this class, only prod.nAssets() = 1 is allowed
     */
    BsMcPricer(SPtrProduct prod,
               SPtrYieldCurve discountYieldCurve,
               double divYield,
               SPtrVolatilityTermStructure volTermStruct,
               double spot,
               McParams mcparams);

    /** Returns the number of variables that can be tracked for stats 
        should be 1, since only track the value, mean and var are both calculated
     */
    size_t nVariables() const;

    /** Runs the simulation and collects statistics,
        Results can be retrieved by statsCalc.results()
        ITER must be double*, StatsCalc.nVariables() must = 1
     */
    template <typename ITER>
    void simulate(StatisticsCalculator<ITER>& statsCalc, unsigned long npaths);
    
    // protected:

    /** Creates and processes one price path.
        It returns the PV of the product
    */
    double processOnePath(Matrix& pricePath);
    
    // state of the product
    SPtrProduct prod_;   // pointer to product
    double divyld_;      // the constant dividend yield

    // the market
    double spot_;        // the spot price
    SPtrYieldCurve discyc_;                        // pointer to yield curve
    SPtrVolatilityTermStructure volterm_;          // pointer to vol term structure

    // state to be precomputed
    Vector discfactors_;  // caches the pre-computed discount factors
    Vector drifts_;       // caches the pre-computed asset drifts
    Vector stdevs_;       // caches the pre-computed standard deviations

    // the final results
    Vector payamts_;      // scratch array for writing the payments after current simulation

    // state of the MC
    McParams mcparams_;   // the Monte Carlo parameters
    SPtrPathGenerator pathgen_; //  pointer to the path generator
};


///////////////////////////////////////////////////////////////////////////////
// Inline definitions

inline
size_t BsMcPricer::nVariables() const {
    return 1;
}


template <typename ITER>
void BsMcPricer::simulate(StatisticsCalculator<ITER>& statsCalc, unsigned long npaths) {
    Matrix pricePath(pathgen_->nTimeSteps(), pathgen_->nFactors());
    ASSERT(statsCalc.nVariables() == 1, "BsMcPricer: statsCalc must have only one variable");
    
    for (unsigned long i = 0; i < npaths; ++i) {
        double pv = processOnePath(pricePath);
        statsCalc.addSample(&pv, &pv + 1);
    }        
};




#endif // BSMCPRICER_HPP
