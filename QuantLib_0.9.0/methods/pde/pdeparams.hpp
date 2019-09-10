#ifndef PDEPARAMS_HPP
#define PDEPARAMS_HPP

#include <vector>
using namespace std;

/** nSpotNodes and nStdDevs are used for X = log S axis. 
    nStdDevs is used to determine Xmin, Xmax:
    (X0 - Xmin) ~ nStdDevs * volX * sqrt(T) 
    nSpotNodes is used to calculate DX
*/

struct PdeParams {
  size_t nTimeSteps;
  vector<size_t> nSpotNodes;   // spot nodes for each asset
  vector<double> nStdDevs;     // num. standard deviations for each asset
  double theta;                // used for adjust between explicit / implicit scheme

  /** Default ctor, n assets, explicit */
  PdeParams(size_t n = 1) : nTimeSteps(1), nSpotNodes(n, 10), nStdDevs(n, 4), theta(0) {};

  PdeParams(size_t nT, size_t Nspot, double stdDev, double tt)
    : nTimeSteps(nT), nSpotNodes(1, Nspot), nStdDevs(1, stdDev), theta(tt) {}
};

#endif // PDEPARAMS_HPP
