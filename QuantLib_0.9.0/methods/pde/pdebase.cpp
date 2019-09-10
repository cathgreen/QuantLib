#include "pdebase.hpp"


void PdeBase::initGrid(double T, const PdeParams& params) {
  ASSERT(nAssets_ == params.nSpotNodes.size(),
         "PdeBase: unequal nAssets and pde parameter axes specs");
  gridAxes_.resize(nAssets_);

  // loop over assets
  for (size_t i = 0; i < nAssets_; ++i) {
    GridAxis& grax = gridAxes_[i];
    grax.NX = params.nSpotNodes[i];
    double S0 = spots_[i];
    spotAxis_.push_back(S0);

    // compute forward to maturity
    double rate = spaccrycs_[i]->spotRate(T);
    double forward = S0 * exp((rate - divyields_[i]) * T);
    double vol = vols_[i]->spotVol(T);

    // initialize the coordinate transform for this axis
    grax.coordinateChange->init(params); // set theta
    double X0 = grax.coordinateChange->fromRealToDiffused(S0);
    double forwardX = forward, volX = vol;
    grax.coordinateChange->forwardAndVariance(forwardX, volX, T);
    grax.coordinateChange->bounds(X0, forwardX, volX, T, params.nStdDevs[i],
                                  grax.Xmin, grax.Xmax);
    grax.DX = (grax.Xmax - grax.Xmin) / (grax.NX + 1);

    // align the grid axis so that a node passes through the alignment value
    double alignValue = grax.coordinateChange->fromRealToDiffused(alignments_[i]);
    int X0NodeIdx = int(0.5 + (alignValue - grax.Xmin) / grax.DX);
    double closestX = grax.Xmin + X0NodeIdx * grax.DX;
    grax.Xmin -= closestX - alignValue;
    grax.Xmax -= closestX - alignValue;

    // fill in the original and the transformed spot nodes
    grax.Xlevels.resize(grax.NX + 2);  // add 2 for the boundary nodes
    grax.Slevels.resize(grax.NX + 2);
    double Xnow = grax.Xmin;
    for (size_t j = 0; j < grax.NX + 2; ++j) {
      grax.Xlevels[j] = Xnow;
      grax.Slevels[j] = grax.coordinateChange->fromDiffusedToReal(grax.Xlevels[j]);
      Xnow += grax.DX;
    }

    // resize the drift, variance and vol vectors
    // no need to add boundary points here
    grax.drifts.resize(grax.NX);
    grax.variances.resize(grax.NX);
    grax.vols.resize(grax.NX);
  } 
}

/** timesteps_[stepIdx] is the current time 
    need to update the drifts, variance, vols
    no need to update Xlevels, Slevels
 */
void PdeBase::updateGrid(const PdeParams& params,
                         const Matrix& fwdFactors, // fwdFactors and fwdVols are constant, caculated once
                         const Matrix& fwdVols,
                         size_t stepIdx) {

  double T1 = timesteps_[stepIdx];
  double T2 = timesteps_[stepIdx + 1];
  double DT = T2 - T1;

  // loop over assets
  for (size_t i = 0; i < nAssets_; ++i) {
    for (size_t j = 1; j <= params.nSpotNodes[i]; ++j) {
       GridAxis& grax = gridAxes_[i];
       double RealS = grax.Slevels[j];
       double aCoeff = fwdFactors(stepIdx, i);
       double RealF = RealS * aCoeff;
       //changed for step localvol
       double RealLNvol = fwdVols(stepIdx, i);
       // set the drift, variance and vol values for this time step
       grax.coordinateChange->driftAndVariance(RealS, RealF, theta_, DT,
                                               RealLNvol, aCoeff, grax.DX,
                                               grax.drifts[j - 1],
                                               grax.variances[j - 1],
                                               grax.vols[j - 1]);
    }
  } 
}

/** Initializes the grid axes gridAxes_, sets up the nodes and the bounds
*/
void PdeBase::solve(const PdeParams& params) {
  // store the Theta
  theta_ = params.theta;

  // get the time steps
  // the resulted timesteps_ size not necessarily = params.nTimeSteps
  spprod_->timeSteps(params.nTimeSteps, timesteps_, stepindex_); 
  nSteps_ = timesteps_.size();

  // set the alignment values to the corresponding spots
  alignments_ = spots_;

  // initialize the grid
  double T = timesteps_.back();
  initGrid(T, params);

  // Pre-compute fwdFactors, used for calculate drifts
  Matrix fwdFactors(nSteps_, nAssets_);
  for (size_t j = 0; j < nAssets_; ++j) {
    SPtrYieldCurve spyc = spaccrycs_[j];
    double divyld = divyields_[j];
    for (size_t i = 0; i < nSteps_ - 1; ++i) {
      double T1 = timesteps_[i];
      double T2 = timesteps_[i + 1];
      double fwdRate = spyc->fwdRate(T1, T2);
      if (quantoFlag_) {
        double fwdAssetVol = vols_[j]->fwdVol(T1, T2);
        double fwdFXVol = fxvols_[j]->fwdVol(T1, T2);
        fwdFactors(i, j) = exp((fwdRate - divyld + correl_ * fwdAssetVol * fwdFXVol) * (T2 - T1));
      }
      else {
        fwdFactors(i, j) = exp((fwdRate - divyld) * (T2 - T1));
      }
    }
  }

  // Pre-compute fwdVols, used for calculate drifts
   Matrix fwdVols(nSteps_, nAssets_);
   for (size_t j = 0; j < nAssets_; ++j) {
     for (size_t i = 0; i < nSteps_ - 1; ++i) {
       double T1 = timesteps_[i];
       double T2 = timesteps_[i + 1];
       fwdVols(i, j) = vols_[j]->fwdVol(T1, T2);
     } 
   }

   // initialize the value layers (grid functions, one per variable to solve)
   initValLayers();

   // evaluate the product at maturity
   evalProduct(nSteps_ - 1);

   // the main loop
   for (ptrdiff_t stepIdx = nSteps_ - 2; stepIdx >= 0; --stepIdx) {
     updateGrid(params, fwdFactors, fwdVols, stepIdx);

     // solve
     double dT = timesteps_[stepIdx + 1] - timesteps_[stepIdx];
     solveFromStepToStep(stepIdx, dT);

     // discount
     double df = spdiscyc_->fwdDiscount(timesteps_[stepIdx], timesteps_[stepIdx + 1]);
     discountFromStepToStep(df);

     // eval product for next iteration
     evalProduct(stepIdx);
   }
   storeResults();
}
