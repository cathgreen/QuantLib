#include "multiassetbsmcpricer.hpp"
#include "../methods/montecarlo/eulerpathgenerator.hpp"
#include "../methods/montecarlo/antitheticpathgenerator.hpp"

MultiAssetBsMcPricer::MultiAssetBsMcPricer(SPtrProduct prod,
                                           SPtrYieldCurve discyc,
                                           const Vector& vols,
                                           const Vector& spots,
                                           const Vector& divylds,
                                           const Matrix& corrMat,
                                           const McParams& mcparams)
  : prod_(prod), discyc_(discyc), vols_(vols), spots_(spots),
    divylds_(divylds), corrMat_(corrMat), mcparams_(mcparams)
{
  size_t nassets = prod_->nAssets();
  ASSERT(nassets == vols_.size(),
         "MultiAssetBsMcPricer: size dismatch between prod.nassets and vols");
  ASSERT(nassets == spots_.size(),
         "MultiAssetBsMcPricer: size dismatch between prod.nassets and spots");
  ASSERT(nassets == divylds_.size(),
         "MultiAssetBsMcPricer: size dismatch between prod.nassets and divylds");
  ASSERT(corrMat_.is_square(),
          "MultiAssetBsMcPricer: corrMat is not symmetric");
  ASSERT(nassets == corrMat_.n_rows,
          "MultiAssetBsMcPricer: corrMat size must equal to nasset x nasset");

  // Get the simulation times and number of assets
  Vector fixtimes = prod_->fixTimes();
  size_t ntimesteps = fixtimes.size();

  // initialize pathgen_
  if (mcparams.pathGenType == McParams::PathGenType::EULER) {
    if (mcparams.urngType == McParams::UrngType::MINSTDRAND)
      pathgen_ = make_shared<EulerPathGenerator<NormalRngMinStdRand>>
        (fixtimes.begin(), fixtimes.end(), nassets, corrMat_);
    else if (mcparams.urngType == McParams::UrngType::MT19937)
      pathgen_ = make_shared<EulerPathGenerator<NormalRngMt19937>>
        (fixtimes.begin(), fixtimes.end(), nassets, corrMat_);
    else if (mcparams.urngType == McParams::UrngType::RANLUX3)
      pathgen_ = make_shared<EulerPathGenerator<NormalRngRanLux3>>
        (fixtimes.begin(), fixtimes.end(), nassets, corrMat_);
    else if (mcparams.urngType == McParams::UrngType::RANLUX4)
      pathgen_ = make_shared<EulerPathGenerator<NormalRngRanLux4>>
        (fixtimes.begin(), fixtimes.end(), nassets, corrMat_);
    else if (mcparams.urngType == McParams::UrngType::SOBOL)
      pathgen_ = make_shared<EulerPathGenerator<NormalRngSobol>>
        (fixtimes.begin(), fixtimes.end(), nassets, corrMat_);
    else
       ASSERT(0, "MultiAssetBsMcPricer: unknown UrngType");
  }
  else
    ASSERT(0, "MultiAssetBsMcPricer: unknown PathGenType");

  if (mcparams.controlVarType == McParams::ControlVarType::ANTITHETIC)
    pathgen_ = make_shared<AntitheticPathGenerator>(pathgen_);

  // precompute the discfactors_
  Vector paytimes = prod->payTimes();
  size_t npaytimes = paytimes.size();
  discfactors_.resize(npaytimes);
  for (size_t i = 0; i < npaytimes; ++i)
    discfactors_(i) = discyc_->discount(paytimes(i));

  // precompute the drifts_ and stdevs_ for each asset
  drifts_.resize(ntimesteps, nassets);
  stdevs_.resize(ntimesteps, nassets);
  double T1 = 0;
  for (size_t i = 0; i < ntimesteps; ++i) {
    double T2 = fixtimes(i);
    double dT = T2 - T1;
    for (size_t j = 0; j < nassets; ++j) {
      double fwdrate = discyc_->fwdRate(T1, T2);
      drifts_(i, j) = (fwdrate - divylds_(j) - 0.5 * vols_(j)) * dT;
      stdevs_(i, j) = vols_(j) * sqrt(dT);
    }
    T1 = T2;
  }

  // intialize payamnts_
  payamts_.resize(prod->payTimes().size());
  payamts_.fill(0);
}

double MultiAssetBsMcPricer::processOnePath(Matrix& pricePath) {
  pathgen_->next(pricePath);
  // Get the simulation times and number of assets
  Vector fixtimes = prod_->fixTimes();
  size_t ntimesteps = fixtimes.size();
  size_t nassets = prod_->nAssets();
  ASSERT(pricePath.n_rows == ntimesteps,
         "MultiAssetBsMcPricer: pricePath row number must equal to prod fix time steps");
  ASSERT(pricePath.n_cols == nassets,
         "MultiAssetBsMcPricer: pricePath col number must equal to nassets");

  Vector spots = spots_;
  for (size_t i = 0; i < ntimesteps; ++i) {
    for (size_t j = 0; j < nassets; ++j) {
      double rndnormal = pricePath(i, j);
      pricePath(i, j) = spots(j) * exp(drifts_(i,j) + stdevs_(i,j) * rndnormal);
      spots(j) = pricePath(i, j);
    }
  }
  
  prod_->eval(pricePath);
  Vector payAmounts = prod_->payAmounts();   // payment of length npaytimes, before discounted
  
  double value = 0;
  size_t npaytimes = discfactors_.size();
  ASSERT(payAmounts.size() == npaytimes,
         "MultiAssetBsMcPricer: payAmounts must equal to npaytimes of prod");
  for (size_t i = 0; i < npaytimes; ++i)
    value += payAmounts(i) * discfactors_(i);

  return value;
}
