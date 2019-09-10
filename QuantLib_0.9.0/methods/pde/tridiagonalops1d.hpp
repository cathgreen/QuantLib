#ifndef TRIDIAGONALOPS1D_HPP
#define TRIDIAGONALOPS1D_HPP

#include "../../exception.hpp"
#include "../../math/matrix.hpp"

#include <algorithm>


/** Utility function for inverting a tridiagonal matrix T, i.e. solving the system T*x=y.
    The matrix is input as three equal size vectors (lower, diag, upper) of size N + 2.
    The right-hand side is input as the vector y of size N + 2.
    The solution vector x is modified in place.
    CAUTION: all vectors must be of size N + 2. This is not checked by the function.

    The first and last elements of lower, diag, upper and y are ignored.
    Only the elements x[1] ... x[N] are modified.
*/
template <typename ARRAY1, typename ARRAY2>
void solveTridiagonal(ARRAY2& x,
                      const ARRAY1& lower,
                      const ARRAY1& diag,
                      const ARRAY1& upper,
                      const ARRAY2& y);

/** Utility function that adjusts the explicit and implicit operators for boundary conditions.
    The adjustment implements constant first derivative in spot space at the edge nodes
    (zero second derivative in spot space)
*/
template <typename EXPOP, typename IMPOP>
void adjustOpsForBoundaryConditions(EXPOP& opExplicit,
                                    IMPOP& opImplicit,
                                    double DX)
{
  // do implicit first
  double lowAdjustmentValueImp = opImplicit.adjustForLowerBoundaryCondition(3, 0.0, DX, 0.0, 0.0);
  double highAdjustmentValueImp = opImplicit.adjustForHigherBoundaryCondition(3, 0.0, DX, 0.0, 0.0);

  double lowAdjustmentValueExp = opExplicit.adjustForLowerBoundaryCondition(3, 0.0, DX, 0.0, 0.0);
  double highAdjustmentValueExp = opExplicit.adjustForHigherBoundaryCondition(3, 0.0, DX, 0.0, 0.0);

  opExplicit.addToLowerVal(-lowAdjustmentValueImp + lowAdjustmentValueExp);
  opExplicit.addToUpperVal(-highAdjustmentValueImp + highAdjustmentValueExp);
}


/** Adjusts the solution at the edge notes */
inline
void applyBoundaryConditions(Matrix& solution)  // solution = ntimesteps x nLayers
{
  size_t n = solution.n_rows - 2;     // n is the num. of interior nodes
  size_t nLayers = solution.n_cols;   // nLayers is the num. of variables
  for (size_t j = 0; j < nLayers; ++j) {
    solution(0, j) = 2 * solution(1, j) - solution(2, j);
    solution(n + 1, j) = 2 * solution(n, j) - solution(n - 1, j);
  }
}


/** Base class representing a tridiagonal operator arising in discretization of
    1-dimensional PDEs.
*/
template <class ARRAY = Vector>
class TridiagonalOp1D
{
public:
  /** default Ctor */
  TridiagonalOp1D(): N_(0), LowerVal_(0), UpperVal_(0) {}

  /** default Dtor */
  virtual ~TridiagonalOp1D() {}

  /** initializing ctor from the three diagonal vectors*/
  TridiagonalOp1D(const ARRAY& lower, const ARRAY& diag, const ARRAY& upper)
    : lower_(lower), diag_(diag), upper_(upper), LowerVal_(0), UpperVal_(0)
  {
    N_ = lower_.size() - 2;
  }

  /** initializing ctor from size and constant values for the three diagonal vectors */
  TridiagonalOp1D(size_t N, double lowerConst, double diagConst, double upperConst)
    : N_(N), LowerVal_(0), UpperVal_(0)
  {
    init(N, lowerConst, diagConst, upperConst);
  }

  void init(size_t N, double lowerConst, double diagConst, double upperConst)
  {
    N_ = N;
    LowerVal_ = UpperVal_ = 0;
    lower_.resize(N + 2);
    lower_.fill(lowerConst);
    diag_.resize(N + 2);
    diag_.fill(diagConst);
    upper_.resize(N + 2);
    upper_.fill(upperConst);
    
  }

  /** Adds to the lower value */
  void addToLowerVal(double lowerVal) { LowerVal_ += lowerVal; }

  /** Adds to the upper value */
  void addToUpperVal(double upperVal) { UpperVal_ += upperVal; }

  /////////////////////////////
  // Adjust Boundary conditions

  /** Adjust for the standard (log-linear interpolation) boundary conditions */
  void adjustStandardBoundaryConditions(double DX) {
    adjustForLowerBoundaryCondition(3, 0.0, DX, 0.0, 0.0);
    adjustForHigherBoundaryCondition(3, 0.0, DX, 0.0, 0.0);
  }

  /** Adjust lower boundary conditions, return extra LowerVal_ */
  double adjustForLowerBoundaryCondition(int degree,
                                         double value,
                                         double DX,
                                         double diagAdjust,
                                         double upAdjust);
  
  /** Adjust upper boundary conditions, return extra UpperVal_ */
  double adjustForHigherBoundaryCondition(int degree,
                                          double value,
                                          double DX,
                                          double diagAdjust,
                                          double lowAdjust);
  
  //////////////////////////////////////
  // Apply or inverse apply the operator

  /** vals, results must of size N + 2, but only results[1] ... results[N] are written */
  template <typename ARRAY1, typename ARRAY2>
  void apply(const ARRAY1& vals, ARRAY2& results) const {
    /** 
    ASSERT(vals.size() == N_ + 2,
           "TridiagonalOp1D: input vector size must equal num. interior pts + 2");
    ASSERT(results.size() == N_ + 2,
           "TridiagonalOp1D: output vector size must equal num. interior pts + 2");
    */
    
    results[1] = LowerVal_ + diag_[1] * vals[1] + upper_[1] * vals[2];
    for (size_t i = 2; i < N_; ++i) 
      results[i] = lower_[i] * vals[i - 1] + diag_[i] * vals[i] + upper_[i] * vals[i + 1];
    results[N_] = lower_[N_] * vals[N_ - 1] + diag_[N_] * vals[N_] + UpperVal_;
  }
  

  template <typename ARRAY1, typename ARRAY2>
  void applyPlus(const ARRAY1& vals, ARRAY2& results) const {
    /**
    ASSERT(vals.size() == N_ + 2,
           "TridiagonalOp1D: input vector size must equal num. interior pts + 2");
    ASSERT(results.size() == N_ + 2,
           "TridiagonalOp1D: output vector size must equal num. interior pts + 2");
    */
    results[1] += LowerVal_ + diag_[1] * vals[1] + upper_[1] * vals[2];
    for (size_t i = 2; i < N_; ++i) 
      results[i] += lower_[i] * vals[i - 1] + diag_[i] * vals[i] + upper_[i] * vals[i + 1];
    results[N_] += lower_[N_] * vals[N_ - 1] + diag_[N_] * vals[N_] + UpperVal_;
  }

  /** vals, results must of size N + 2, but only results[1] ... results[N] are written */
  template <typename ARRAY1, typename ARRAY2>
  void applyInverse(const ARRAY1& vals, ARRAY2& results) const {
    /**
    ASSERT(vals.size() == N_ + 2,
           "TridiagonalOp1D: input vector size must equal num. interior pts + 2");
    ASSERT(results.size() == N_ + 2,
           "TridiagonalOp1D: output vector size must equal num. interior pts + 2");    
    */
    solveTridiagonal(results, lower_, diag_, upper_, vals);
  }


  ///////////
  // Algebra

  /** Adds the rhs in place */
  template <typename ARRAY1>
  TridiagonalOp1D & operator+=(const TridiagonalOp1D<ARRAY1> & rhs);

  /** Subtracts the rhs in place */
  template <typename ARRAY1>
  TridiagonalOp1D & operator-=(const TridiagonalOp1D<ARRAY1> & rhs);

  /** Multiplies this operator by the rhs scalar */
  TridiagonalOp1D & operator*=(double rhs);

  template <typename ARRAY1>
  friend TridiagonalOp1D<ARRAY1> operator+(const TridiagonalOp1D<ARRAY1> & a,
                                          const TridiagonalOp1D<ARRAY1> & b);
  template <typename ARRAY1>
  friend TridiagonalOp1D<ARRAY1> operator-(const TridiagonalOp1D<ARRAY1> & a,
                                          const TridiagonalOp1D<ARRAY1> & b);

  template <typename ARRAY1>
  friend TridiagonalOp1D<ARRAY1> operator*(double coeff,
                                          const TridiagonalOp1D<ARRAY1> & a);
  
protected:
  size_t N_;     // num. of interior points
  ARRAY lower_, diag_, upper_;  // all of them have size N_+2

private:
  double LowerVal_, UpperVal_;
};


/** The identity operator */
template <typename ARRAY = Vector >
class IdentityOp1D : public TridiagonalOp1D<ARRAY> {
public:
  /** Ctor, has N interior pts */
  IdentityOp1D(size_t N): TridiagonalOp1D<ARRAY>(N, 0, 1, 0) {}  
};


/** The Delta operator */
template <typename ARRAY = Vector >
class DeltaOp1D : public TridiagonalOp1D<ARRAY> {
public:
  /** Default Ctor */
  DeltaOp1D() = default;

  /** Ctor for Delta Op x theta, drifts only has interior pts, [0, N - 1] */
  template <typename ARRAY2>
  DeltaOp1D(const ARRAY2& drifts, double DT, double DX, double theta) {  
    init(drifts, DT, DX, theta);
  }

  template <typename ARRAY2>
  void init(const ARRAY2& drifts, double DT, double DX, double theta) {
    size_t N = drifts.size();
    TridiagonalOp1D<ARRAY>::init(N, 0, 0, 0);
    double tmp, f1 = theta * DT / (2 * DX);
    for (size_t i = 1; i <= N; ++i) {
      tmp = f1 * drifts[i - 1];
      TridiagonalOp1D<ARRAY>::lower_[i] = -tmp;
      TridiagonalOp1D<ARRAY>::upper_[i] = tmp;
    }
  }
};


/** The Gamma operator */
template <typename ARRAY = Vector >
class GammaOp1D : public TridiagonalOp1D < ARRAY > {
public:
  /** Default Ctor */
  GammaOp1D() = default;

  /** Ctor for Delta Op x theta, drifts only has interior pts, [0, N - 1] */
  template <typename ARRAY2>
  GammaOp1D(const ARRAY2& variances, double DT, double DX, double theta) {
    init(variances, DT, DX, theta);
  }

  template <typename ARRAY2>
  void init(const ARRAY2& variances, double DT, double DX, double theta) {
     size_t N = variances.size();
     TridiagonalOp1D<ARRAY>::init(N, 0, 0, 0);
     double tmp, f1 = theta * DT / (2 * DX * DX);
     for (size_t i = 1; i <= N; ++i) {
       tmp = f1 * variances[i - 1];
       TridiagonalOp1D<ARRAY>::lower_[i] = tmp;
       TridiagonalOp1D<ARRAY>::diag_[i] = -2 * tmp;
       TridiagonalOp1D<ARRAY>::upper_[i] = tmp;
     }
  }
  
};


///////////////////////////////////////////////////////////////////////////////
// Inline definitions

template <typename ARRAY>
template <typename ARRAY1>
inline
TridiagonalOp1D<ARRAY> &
TridiagonalOp1D<ARRAY>::operator+=(const TridiagonalOp1D<ARRAY1> & rhs) {
  ASSERT(N_ == rhs.N_,
         "TridiagonalOperator1D: cannot add two operators of different sizes");
  for (size_t i = 0; i < lower_.size(); ++i) {
    lower_[i] += rhs.lower_[i];
    diag_[i] += rhs.diag_[i];
    upper_[i] += rhs.upper_[i];
  }
  return *this;
}



template<typename ARRAY>
template<typename ARRAY1>
inline
TridiagonalOp1D<ARRAY> &
TridiagonalOp1D<ARRAY>::operator-=(TridiagonalOp1D<ARRAY1> const& rhs)
{
  ASSERT(N_ == rhs.N_,
         "TridiagonalOperator1D: Cannot subtract two operators of different sizes");
  for (size_t i = 0; i < lower_.size(); ++i) {
    lower_[i] -= rhs.lower_[i];
    diag_[i] -= rhs.diag_[i];
    upper_[i] -= rhs.upper_[i];
  }
  return *this;
}



template<typename ARRAY>
inline
TridiagonalOp1D<ARRAY> &
TridiagonalOp1D<ARRAY>::operator*=(double rhs)
{
  for (size_t i = 0; i < lower_.size(); ++i) {
    lower_[i] *= rhs;
    diag_[i] *= rhs;
    upper_[i] *= rhs;
  }
  return *this;
}

template <typename ARRAY>
inline
TridiagonalOp1D<ARRAY> operator+(const TridiagonalOp1D<ARRAY> & a,
                                 const TridiagonalOp1D<ARRAY> & b) {
  ASSERT(a.N_ == b.N_,
         "TridiagonalOperator1D: Cannot add two operators of different sizes");
  size_t N = a.N_;
  TridiagonalOp1D<ARRAY> res(N, 0, 0, 0);
  for (size_t i = 1; i <= N; ++i) {
    res.lower_[i] = a.lower_[i] + b.lower_[i];
    res.upper_[i] = a.upper_[i] + b.upper_[i];
    res.diag_[i] = a.diag_[i] + b.diag_[i];
  }
  return res;
}

template <typename ARRAY>
inline
TridiagonalOp1D<ARRAY> operator-(const TridiagonalOp1D<ARRAY> & a,
                                 const TridiagonalOp1D<ARRAY> & b) {
  ASSERT(a.N_ == b.N_,
         "TridiagonalOperator1D: Cannot subtract two operators of different sizes");
  size_t N = a.N_;
  TridiagonalOp1D<ARRAY> res(N, 0, 0, 0);
  for (size_t i = 1; i <= N; ++i) {
    res.lower_[i] = a.lower_[i] - b.lower_[i];
    res.upper_[i] = a.upper_[i] - b.upper_[i];
    res.diag_[i] = a.diag_[i] - b.diag_[i];
  }
  return res;
}




template <typename ARRAY>
inline
TridiagonalOp1D<ARRAY> operator*(double coeff,
                                 const TridiagonalOp1D<ARRAY> & a) {
  TridiagonalOp1D<ARRAY> res(a);
  res *= coeff;
  return res;
}



template <typename ARRAY1, typename ARRAY2>
void solveTridiagonal(ARRAY2& x,
                      const ARRAY1& lower,
                      const ARRAY1& diag,
                      const ARRAY1& upper,
                      const ARRAY2& y) {
  // all arguments should of size n + 2
  size_t n = diag.size() - 2;

  static Vector D(1);    // initialize once, then store for future function call
  static Vector Y(1); 
  
  if (D.size() != n + 1)
    D = Vector(n + 1);
  if (Y.size() != n + 1)
    Y = Vector(n + 1);
  D[n] = diag[n];
  Y[n] = y[n];

  for (size_t i = n - 1; i >= 1; --i) {
    D[i] = diag[i] - upper[i] * lower[i + 1] / D[i + 1];
    Y[i] = y[i] - upper[i] * Y[i + 1] / D[i + 1];
  }

  x[1] = Y[1] / D[1];

  for (size_t i = 2; i <= n; ++i)
    x[i] = (Y[i] - lower[i - 1] * x[i - 1]) / D[i];
}



template<typename ARRAY>
inline
double TridiagonalOp1D<ARRAY>::adjustForLowerBoundaryCondition(int degree,
                                                               double value,
                                                               double DX,
                                                               double diagAdjust,
                                                               double upAdjust)
{
  ASSERT(diag_.size() >= 4, "TridiagonalOperator1D: grid is too small!");
  switch (degree) {
  case 0:
    return value;       // we set the actual value
    break;
  case 1:
    upper_[1] += lower_[1];
    return -lower_[1] * value;
    break;
  case 2:
    diag_[1] += 2 * lower_[1];
    upper_[1] -= lower_[1];
    return lower_[1] * value;
    break;
  case 3:
    ASSERT(value == 0.0,
           "TridiagonalOperator1D: cannot do non-zero 2nd derivative boundary condition");
    diag_[1] += 2 / (1.0 + DX / 2.0) * lower_[1];
    upper_[1] -= (1.0 - DX / 2.0) / (1.0 + DX / 2.0) * lower_[1];
    return 0;
    break;
  case 4:
    ASSERT(value == 0.0,
           "TridiagonalOperator1D: cannot do non-zero 2nd derivative boundary condition");
    diag_[1] += diagAdjust * lower_[1];
    upper_[1] += upAdjust * lower_[1];
    return 0;
    break;
  default:
    ASSERT(0, "TridiagonalOperator1D: invalid degree for boundary condition");
  }
}




template<typename ARRAY>
inline
double TridiagonalOp1D<ARRAY>::adjustForHigherBoundaryCondition(int degree,
                                                               double value,
                                                               double DX,
                                                               double diagAdjust,
                                                               double lowAdjust)
{
  ASSERT(diag_.size() >= 4, "TridiagonalOperator1D: grid is too small!");
  switch (degree) {
  case 0:
    return value;  // we set the actual value
  case 1:
    lower_[N_] += upper_[N_];
    return -upper_[N_] * value;
  case 2:
    diag_[N_] += 2 * upper_[N_];
    lower_[N_] -= upper_[N_];
    return upper_[N_] * value;
  case 3:
    ASSERT(value == 0.0,
           "TridiagonalOperator1D: cannot do non-zero 2nd derivative boundary condition");
    diag_[N_] += 2.0 / (1.0 - DX / 2.0) * upper_[N_];
    lower_[N_] -= (1.0 + DX / 2.0) / (1.0 - DX / 2.0) * upper_[N_];
    return 0.0;
  case 4:
    ASSERT(value == 0.0,
           "TridiagonalOperator1D: cannot do non-zero 2nd derivative boundary condition");
    diag_[N_] += diagAdjust * upper_[N_];
    lower_[N_] += lowAdjust  * upper_[N_];
    return 0.0;
  default:
    ASSERT(0, "TridiagonalOperator1D: invalid degree for boundary condition");
  }
}


#endif // TRIDIAGONALOPS1D_HPP
