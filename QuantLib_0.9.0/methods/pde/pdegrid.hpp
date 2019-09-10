#ifndef PDEGRID_HPP
#define PDEGRID_HPP


#include "../../exception.hpp"
#include "../../math/matrix.hpp"

#include <algorithm>
#include <memory>

/** Abstract coordinate change interface. 
    Note: only member variable is theta
*/
class CoordinateChangeBase {
public:

  double theta;  // calculate drift needs theta

  virtual ~CoordinateChangeBase() {}

  template <typename PARAMS>
  void init(const PARAMS& params) {
    theta = params.theta;
  }

  virtual double fromRealToDiffused(double S1) = 0;

  virtual double fromDiffusedToReal(double S1) = 0;

  /** Change the fwd price and vol in new coordinate */
  virtual void forwardAndVariance(double& fwd, double& vol, double T) = 0;

  /** Calculate drift and variance in new coordinate */
  virtual void driftAndVariance(double realS,
                                double realF,   // real are in original coord
                                double theta,
                                double DT,
                                double realLNVol,
                                double aCoeff,
                                double DX,
                                double& drift,
                                double& variance,
                                double& finalVol) = 0; 

  /** Computes the grid bounds Xmin and Xmax */
  virtual void bounds(double S0,
                      double fwd,
                      double vol,
                      double T,    // total time, (X0 - Xmin) ~ nstds * vol * sqrt(T)
                      double nstds,
                      double& Xmin,
                      double& Xmax) = 0;
};


/** No coordinate change, i.e. Diffused = Real */
class NoCoordinateChange: public CoordinateChangeBase {
public:
  
  virtual double fromRealToDiffused(double S) {
    return S;
  }

  virtual double fromDiffusedToReal(double X) {
    return X;
  }

  /** Change the fwd price and vol in new coordinate */
  virtual void forwardAndVariance(double& fwd, double& vol, double T) {
    // nothing to do 
    return;
  }

  /** Calculate drift and variance in new coordinate */
  virtual void driftAndVariance(double realS,
                                double realF,   
                                double theta,
                                double DT,
                                double realLNVol,
                                double aCoeff,
                                double DX,
                                double& drift,
                                double& variance,
                                double& finalVol) {

    double vol = realLNVol * realS;     // vol = S * LNvol
    double corr = (theta * aCoeff + 1 - theta);
    drift = (realF - realS) / corr / DT;
    variance = vol * vol;           
    finalVol = sqrt(variance);
  }

  /** Computes the grid bounds Smin and Smax */
  virtual void bounds(double S0,
                      double F,     // fwd in log
                      double vol,
                      double T,    
                      double nstds,
                      double& Smin,
                      double& Smax) {
    Smin =  min(S0, F) * exp(-0.5 * vol * vol * T - nstds * vol * sqrt(T));
    Smax = max(S0, F) * exp(-0.5 * vol * vol * T + nstds * vol * sqrt(T));
  }
};



/** Logarithmic coordinate change, i.e. Diffused = log(Real) */
class LogCoordinateChange: public CoordinateChangeBase {
public:
  
  virtual double fromRealToDiffused(double S) {
    return log(S);
  }

  virtual double fromDiffusedToReal(double X) {
    return exp(X);
  }

  /** Change the fwd price and vol in new coordinate */
  virtual void forwardAndVariance(double& fwd, double& vol, double T) {
    fwd = log(fwd) - 0.5 * vol * vol * T;
  }

  /** Calculate drift and variance in new coordinate */
  virtual void driftAndVariance(double realS,
                                double realF,   
                                double theta,
                                double DT,
                                double realLNVol,
                                double aCoeff,
                                double DX,
                                double& drift,
                                double& variance,
                                double& finalVol) {

    double Xi = fromRealToDiffused(realS);
    double Deltaip1 = (fromDiffusedToReal(Xi + DX) - fromDiffusedToReal(Xi - DX)) / (2.0 * DX);
    double Gammaip1 = (fromDiffusedToReal(Xi + DX) - 2 * fromDiffusedToReal(Xi)
                       + fromDiffusedToReal(Xi - DX));
    Gammaip1 /= (DX * DX);
    double corr = (theta * aCoeff + 1 - theta);
    drift = (realF - realS) / corr / DT / Deltaip1
      - 0.5 * realLNVol * realLNVol * Gammaip1 / Deltaip1;
    variance = realLNVol * realLNVol;  // variance and vol in log is same as in linear
    finalVol = realLNVol;
  }

  /** Computes the grid bounds Xmin and Xmax */
  virtual void bounds(double X0,
                      double F,     // fwd in log
                      double vol,
                      double T,    
                      double nstds,
                      double& Xmin,
                      double& Xmax) {
    Xmin = min(X0, F) - nstds * vol * sqrt(T);
    Xmax = max(X0, F) + nstds * vol * sqrt(T);
  }
};



/** Describes the discretization of a grid coordinate axis
    with a mbmer a pointer to CoordinateChangeBase, to do coord change
    Note: only store the X grid and corresponding drifts, variances, 
    vols on the grid at ONE timestep
 */
class GridAxis {
public:
  double Xmin, Xmax, DX;     // max, min and distance between nodes
  size_t NX;                 // number of interior nodes
  Vector Xlevels, Slevels;
  Vector drifts, variances, vols;
  shared_ptr<CoordinateChangeBase> coordinateChange;
  // the coordinate change rules for this axis

  /** Default ctor uses logarithmic coordinate changes */
  GridAxis()
    : coordinateChange(new LogCoordinateChange()) {}

  /** Sets the coordinate changes */
  void setCoordinateChange(shared_ptr<CoordinateChangeBase> c) {
    coordinateChange = c;
  }
};

#endif // PDEGRID_HPP
