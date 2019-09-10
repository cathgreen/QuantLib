#ifndef PDERESULTS_HPP
#define PDERESULTS_HPP

#include "pdegrid.hpp"
#include <vector>


class PdeResults {
public:
  Vector prices;   // vector of size nLayers, with the prices at the current spots
  Vector times;    // vector of time nodes
  vector<GridAxis> gridAxes; // vector of size nAssets with the grid axes

  /** Computes the spot axis (linear coord) for the asset with index assetIdx */
  void getSpotAxis(size_t assetIdx, Vector& axis) {
    ASSERT(assetIdx < gridAxes.size(), "PdeResults: assetIdx out of range");
    GridAxis& gridAxis = gridAxes[assetIdx];
    axis.resize(gridAxis.NX + 2);   // axis includes boundary nodes
    double xnow = gridAxis.Xmin, xmin = gridAxis.Xmin;
    double dx = gridAxis.DX;
    for (size_t i = 0; i < axis.size(); ++i) {
      axis[i] = gridAxis.coordinateChange->fromDiffusedToReal(xnow);
      xnow += dx;
    }
  }
 
  size_t nAssets() const {
    return gridAxes.size();
  }
  
protected:
  /** Computes the spot axes for all assets */
  void computeSpotAxes() {
    size_t nAssets = gridAxes.size();
    spotAxes.resize(nAssets);
    for (size_t i = 0; i < nAssets; ++i)
      getSpotAxis(i, spotAxes[i]);
  }
  
  // state
  vector<Vector> spotAxes;  // a vector of spot axes vectors, could be jagged, get from gridAxes
};


/** Pde1D nLayers = 1, nAssets = 1 
 */
class Pde1DResults : public PdeResults {
public:
  vector<Matrix> values;   // a vector of ntimesteps, each matrix is of nSpots x nLayers
                           // is the values of nLayers variables for a given time
  
   /** Returns the vector of times, the vector of spots and the matrix of values for
      a variable with index varIdx
      Note: zValues = ntimes x nLayers, sAxis of Nspots for the one asset
  */
  void getValues(size_t varIdx, Vector& timeAxis, Vector& sAxis, Matrix& zValues) {
    ASSERT(varIdx == 0, "Pde1DResults: can only check one variable");
    timeAxis = times;
    getSpotAxis(0, sAxis);
    zValues.resize(times.size(), sAxis.size());

    for (size_t i = 0; i < times.size(); ++i) {
      for (size_t j = 0; j < sAxis.size(); ++j) {
        zValues(i, j) = values[i](j, varIdx);
      }
    }
  }
};

#endif // PDERESULTS_HPP
