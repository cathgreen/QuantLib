#include "bsmcpricer.hpp"
#include "../methods/montecarlo/eulerpathgenerator.hpp"

#include <cmath>

BsMcPricer::BsMcPricer(SPtrProduct prod,
                       SPtrYieldCurve discountYieldCurve,
                       double divYield,
                       SPtrVolatilityTermStructure volTermStruct,
                       double spot,
                       McParams mcparams)
    : prod_(prod), divyld_(divYield), spot_(spot),
      discyc_(discountYieldCurve), volterm_(volTermStruct),
      mcparams_(mcparams)
{
    // Check if input are legal
    ASSERT(spot_ >= 0, "BsMcPricer: spot must be non-negative");
    ASSERT(divyld_ >= 0, "BsMcPricer: divYield must be non-negative");
    ASSERT(prod_->nAssets() == 1, "BsMcPricer: only single asset product is allowed")
 
    // initialize pathgen_
    const Vector& fixtimes = prod_->fixTimes();
    size_t ntimesteps = fixtimes.size();
    
    if (mcparams.pathGenType == McParams::PathGenType::EULER) {
        if (mcparams.urngType == McParams::UrngType::MINSTDRAND)
            pathgen_ = SPtrPathGenerator(new EulerPathGenerator<NormalRngMinStdRand>(ntimesteps, 1));
        else if (mcparams.urngType == McParams::UrngType::MT19937)
            pathgen_ = SPtrPathGenerator(new EulerPathGenerator<NormalRngMt19937>(ntimesteps, 1));
        else if (mcparams.urngType == McParams::UrngType::RANLUX3)
            pathgen_ = SPtrPathGenerator(new EulerPathGenerator<NormalRngRanLux3>(ntimesteps, 1));
        else if (mcparams.urngType == McParams::UrngType::RANLUX4)
            pathgen_ = SPtrPathGenerator(new EulerPathGenerator<NormalRngRanLux4>(ntimesteps, 1));
        else
            ASSERT(0, "BsMcPricer: mcparams has invalid urngType!");
        }
    else
        ASSERT(0, "BsMcPricer: mcparams has invalid pathGenType");

    // Pre-compute the discount factor at all payTimes of prod
    const Vector& paytimes = prod_->payTimes();
    discfactors_.resize(paytimes.size());
    for (size_t i = 0; i < paytimes.size(); ++i) 
        discfactors_(i) = discyc_->discount(paytimes(i));

    // Pre-compute the stdevs and drifts at all fixTimes of prod
    drifts_.resize(ntimesteps);
    stdevs_.resize(ntimesteps);
    double t1 = 0;
    for (size_t i = 0; i < ntimesteps; ++i) {
        double t2 = fixtimes(i);
        double fwdrate = discyc_->fwdRate(t1, t2);
        double fwdvol = volterm_->fwdVol(t1, t2);
        drifts_(i) = (fwdrate - divyld_ - 0.5 * fwdvol * fwdvol) * (t2 - t1) ;
        stdevs_(i) = fwdvol * sqrt(t2 - t1);
        t1 = t2;
    }
    
    // Resize payamts_
    payamts_.resize(prod_->payAmounts().size());
    payamts_.fill(0);
}


double BsMcPricer::processOnePath(Matrix& pricePath) {
    size_t ntimesteps = pathgen_->nTimeSteps();
    pricePath.resize(pathgen_->nTimeSteps(), 1);
    pathgen_->next(pricePath);  // normal distribution with mean 0 and stdev 1 as seed
    double St = spot_;          // current price
    for (size_t i = 0; i < ntimesteps; ++i) {
        double rd = pricePath(i, 0);
        St *= exp(drifts_(i) + stdevs_(i) * rd);
        pricePath(i, 0) = St;
    }
    prod_->eval(pricePath);
    payamts_ = prod_->payAmounts();

    double pv = 0;
    for (size_t i = 0; i < payamts_.size(); ++i)
        pv += payamts_(i) * discfactors_(i);
    
    return pv;
}
