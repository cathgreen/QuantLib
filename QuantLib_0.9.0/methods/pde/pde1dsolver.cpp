#include "pde1dsolver.hpp"
#include "../../math/interpol/interpolation1d.hpp"

/** Solves backwards from one time step to the previous */
void Pde1DSolver::solveFromStepToStep(ptrdiff_t step, double DT) {

  // initialise operators
  GridAxis& grax = gridAxes_[0];
  deltaOpExplicit_.init(grax.drifts, DT, grax.DX, 1.0 - theta_);
  deltaOpImplicit_.init(grax.drifts, DT, grax.DX, theta_);

  gammaOpExplicit_.init(grax.variances, DT, grax.DX, 1.0 - theta_);
  gammaOpImplicit_.init(grax.variances, DT, grax.DX, theta_);

  // build the explicit and implicit operators
  opExplicit_.init(grax.NX, 0.0, 1.0, 0.0); // initialize to identity matrix
  opExplicit_ += deltaOpExplicit_;
  opExplicit_ += gammaOpExplicit_;
  opImplicit_.init(grax.NX, 0.0, 1.0, 0.0); // initialize to identity matrix
  opImplicit_ -= deltaOpImplicit_;
  opImplicit_ -= gammaOpImplicit_;

  // adjust the operators for boundary conditions
  adjustOpsForBoundaryConditions(opExplicit_, opImplicit_, grax.DX);

   // Main loop over the layers
  for (size_t j = 0; j < nLayers_; ++j) {
    // NOTE: v1 and v2 are read-write views into the corresponding columns
    // They are not independent copies, so we are modifying in place prevValues and currValues
    auto v1 = prevValues->col(j);
    auto v2 = currValues->col(j);
    opExplicit_.apply(v1, v2);
    opImplicit_.applyInverse(v2, v1);
  }

  // apply boundary coditions to solution
  applyBoundaryConditions(*prevValues);
}


/** Initializes the layers, i.e., preValues, currValues
    Each layer (grid function) corresponds to variable solved on the grid.
*/
void Pde1DSolver::initValLayers() {
  ASSERT(gridAxes_.size() == 1,
         "Pde1DSolver: can only handle one asset only!");
  
  // initiliaze all values to zero
  prevValues = make_shared<Matrix>(gridAxes_[0].NX + 2, nLayers_, arma::fill::zeros);
  currValues = make_shared<Matrix>(gridAxes_[0].NX + 2, nLayers_, arma::fill::zeros);
  
  // prepare the results
  results_.times.resize(nSteps_);
  results_.values.resize(nSteps_);
  results_.prices.resize(nLayers_);
}

/** Evaluates the product at the passed-in time step index */
void Pde1DSolver::evalProduct(size_t stepIdx) {
  ptrdiff_t eventIdx = stepindex_[stepIdx];
  if (eventIdx >= 0) {   // product event, must evaluate;
                         // otherwise, just equal value at eventIdx = 0 discounted back
    const Vector& payTms = spprod_->payTimes();      // the payment times
    double fixtime = spprod_->fixTimes()[eventIdx];  // the current fixing time
    Vector spots(1);     // one underlying spot
    for (size_t node = 0; node <= gridAxes_[0].NX + 1; ++node) {
      spots[0] = gridAxes_[0].Slevels[node];
      spprod_->eval(eventIdx, spots, (*prevValues)(node, 0));   // need eventIdx do determine if the current price should be used for calculating price

      // read out the amounts
      const Vector& payAms = spprod_->payAmounts();
      (*prevValues)(node, 0) = payAms[eventIdx]; // assume payments and fixtimes of prod are same
    }
  }
  results_.times[stepIdx] = timesteps_[stepIdx];
  if (storeAllResults_)
    results_.values[stepIdx] = *prevValues;
}


/** Stores the solver results */
void Pde1DSolver::storeResults() {
  results_.gridAxes = gridAxes_;
   for (size_t j = 0; j < nLayers_; ++j) {
     double X0 = gridAxes_[0].coordinateChange->fromRealToDiffused(spots_[0]);
     Vector temp(prevValues->col(j));
     LinearInterpolation1D<Vector> interp(gridAxes_[0].Xlevels, temp);
     results_.prices[j] = interp.getValue(X0);
  }
}


/** Discounts the grid functions on the current time step, by applying
    the passed-in one-step discount factor. */
void Pde1DSolver::discountFromStepToStep(double df) {
  *prevValues *= df;
}
