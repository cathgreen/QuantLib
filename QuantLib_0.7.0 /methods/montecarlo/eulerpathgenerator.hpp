#ifndef EULERPATHGENERATE_HPP
#define EULERPATHGENERATE_HPP

#include "pathgenerator.hpp"
#include "../../math/random/rng.hpp"
#include "../../math/matrix.hpp"

/** PathGenerator to fill in a ntimestpes * nfactors matrix with 
    N(0, 1) random numbers, for further use to calculate price path
*/
template <typename NRNG>
class EulerPathGenerator: public PathGenerator {
public:
    /** Ctor for EulterPathGenerator for independent factors 
        EulerPathGenerator(size_t ntimesteps, size_t nfactors)
        : PathGenerator(ntimesteps, nfactors), nrng_(ntimesteps, 0, 1) {} 
    */

    
    /** Ctor for generating increments for correlated factors.
        If the correlation matrix is not passed in, it assumes independent factors
    */
    template <typename ITER>
    EulerPathGenerator(ITER timestepsBegin, ITER timestepsEnd, size_t nfactors,
                       const Matrix& correlMat = Matrix());
   

    /** Returns the dimension of random number generator nrng_
        should = ntimesteps_
    */
    size_t dim() const {
        return nrng_.dim();
    }

     /** Write ont he next price path on pricePath
        = ntimesteps_ * nfactors_, if size not match, then resize pricePath
    */
    virtual void next(Matrix& pricePath) override;
    
protected:
    
    // state
    NRNG nrng_;  // rng class deinfed in math/random/rng.hpp
    Vector sqrtDeltaT_;              // sqrt(T1), sqrt(T2-T1), ... I don't think it is necessary here
    Vector normalDevs_;              // scratch array
};

///////////////////////////////////////////////////////////////////////////////
// Inline definitions

template <typename NRNG>
template <typename ITER>
EulerPathGenerator<NRNG>::EulerPathGenerator(ITER timestepsBegin, ITER timestepsEnd,
                                             size_t nfactors, const Matrix& correlMat)
    : PathGenerator((size_t)(timestepsEnd - timestepsBegin), nfactors, correlMat),
      // nrng_((size_t)(timestepsEnd - timestepsBegin) * nfactors, 0, 1)
      nrng_((size_t)(timestepsEnd - timestepsBegin), 0, 1)
{
    ASSERT(ntimesteps_ > 0, "EulerPathGenerator: no time steps");
    normalDevs_.resize(ntimesteps_);
    sqrtDeltaT_.resize(ntimesteps_);
    sqrtDeltaT_(0) = sqrt(*timestepsBegin);
    ITER it = timestepsBegin;
    for (size_t i = 1; i < ntimesteps_; ++i, ++it) {
        double deltaT = *(it+1) - *it;
        ASSERT(deltaT > 0,
               "EulerPathGenerator: time steps are not unique or not in increasing order");
        sqrtDeltaT_(i) = sqrt(deltaT);
    }
}

template <typename NRNG>
void EulerPathGenerator<NRNG>::next(Matrix& pricePath) {
    pricePath.resize(ntimesteps_, nfactors_);
    // iterate over columns; the matrix will be filled column by column
    for (size_t j = 0; j < nfactors_; ++j) {
        nrng_.next(normalDevs_.begin(), normalDevs_.end());
        pricePath.col(j) = normalDevs_;
    }

    
     // finally apply the Cholesky factor if not empty
    if (!sqrtCorrel_.is_empty()) {
        pricePath *= sqrtCorrel_.t();
    }
}


#endif // EULERPATHGENERATE_HPP
