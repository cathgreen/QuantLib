#ifndef STATISTICSCALCULATOR_HPP
#define STATISTICSCALCULATOR_HPP

#include "../../defines.hpp"
#include "../../exception.hpp"
#include "../matrix.hpp"

/** Statistics calculator base class */
template <typename ITER>
class StatisticsCalculator {
public:
    StatisticsCalculator(size_t nresults, size_t nvars);

    virtual ~StatisticsCalculator() {}

    /** Adds one sample; requires end - big == nVariables() */
    virtual void addSample(ITER begin, ITER end) = 0;

    /** Clears samples and results, but keep nresults, nvars*/
    virtual void reset();

    /** Returns the number of samples addes so far */
    size_t nSample() const {
        return nsamples_;
    }

    /** Returns the number of variables */
    virtual size_t nVariables() const {
        return results_.n_cols;
    }

    /** Returns the results, one column per variable */
    virtual const Matrix& results() = 0;
    
    
protected:
    size_t nsamples_;         // number of observations to calculate results
    // mutable Matrix results_;  // a matrix of nresults x nvars to store statistics
                              // results(i, j) is the statistic i for variable j
    Matrix results_;
};


///////////////////////////////////////////////////////////////////////////////
// Inline definitions

template <typename ITER>
inline StatisticsCalculator<ITER>::StatisticsCalculator(size_t nresults, size_t nvars)
    : nsamples_(0), results_()
{
    ASSERT(nresults > 0, "StatisticalCalculator: nresults must be positive");
    ASSERT(nvars > 0, "StatisticalCalculator: nvars must be positive");
    results_.resize(nresults, nvars);
    results_.fill(0);   // default fill in zeros
}

template <typename ITER>
void StatisticsCalculator<ITER>::reset() {
    nsamples_ = 0;
    results_.fill(0);
}


#endif // STATISTICSCALCULATOR_HPP
