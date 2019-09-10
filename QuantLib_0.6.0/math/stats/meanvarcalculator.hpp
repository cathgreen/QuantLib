#ifndef MEANVARCALCULATOR_HPP
#define MEANVARCALCULATOR_HPP

#include "statisticscalculator.hpp"

/** Statistics calculator for mean and variance */
template <typename ITER>
class MeanVarCalculator: public StatisticsCalculator<ITER> {
public:
    MeanVarCalculator(size_t nvars);

    virtual void addSample(ITER begin, ITER end) override;

    virtual void reset() override;

    virtual const Matrix& results() override;

    // there is a implicit virtual destructor heritate from Base class
    
protected:

    // state
    Vector runningSum_;   // sum of the samples have seen
    Vector runningSum2_;  // sum of the square of samples have seen
    
};


///////////////////////////////////////////////////////////////////////////////
// Inline definitions

template <typename ITER>
MeanVarCalculator<ITER>::MeanVarCalculator(size_t nvars)
    : StatisticsCalculator<ITER>(2, nvars), runningSum_(nvars), runningSum2_(nvars)
{
    runningSum_.fill(0);
    runningSum2_.fill(0);
}


template <typename ITER>
void MeanVarCalculator<ITER>::addSample(ITER begin, ITER end) {
    ASSERT((size_t)(end - begin) == this->nVariables(), "MeanVarCalculator: new sample size not match");
    ITER it = begin;
    for (size_t i = 0; i < this->nVariables(); ++i, ++it) {
        runningSum_(i) += *it;
        runningSum2_(i) += (*it) * (*it);
    }
    ++this->nsamples_;
}

template <typename ITER>
void MeanVarCalculator<ITER>::reset() {
    this->nsamples_ = 0;
    runningSum_.fill(0);
    runningSum2_.fill(0);
    this->results_.fill(0);
}

template <typename ITER>
const Matrix& MeanVarCalculator<ITER>::results() {
    ASSERT(this->nsamples_ > 0, "MeanVarCalculator: no sample to calculate mean & variance");
    for (size_t i = 0; i < this->nVariables(); ++i) {
        this->results_(0, i) = runningSum_(i) / this->nsamples_;
        double mean = this->results_(0, i);
        this->results_(1, i) = (runningSum2_(i) - this->nsamples_ * mean * mean) / (this->nsamples_ - 1);
    }
    return this->results_;
}


    
#endif // MEANVARCALCULATOR_HPP
