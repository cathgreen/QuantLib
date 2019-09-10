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
    /** Ctor for EulterPathGenerator for independent factors */
    EulerPathGenerator(size_t ntimesteps, size_t nfactors)
        : PathGenerator(ntimesteps, nfactors), nrng_(ntimesteps, 0, 1) {}

    
    /** Ctor for generating increments for correlated factors.
        If the correlation matrix is not passed in, it assumes independent factors
    
    template <typename ITER>
    EulerPathGenerator(ITER timestepsBegin, ITER timestepsEnd, size_t nfactors,
                       const Matrix& correlMat = Matrix());
    */

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
    
};

///////////////////////////////////////////////////////////////////////////////
// Inline definitions

template <typename NRNG>
void EulerPathGenerator<NRNG>::next(Matrix& pricePath) {
    pricePath.resize(ntimesteps_, nfactors_);
    for (size_t j = 0; j < nfactors_; ++j) {
        auto itb = pricePath.begin_col(j);
        auto ite =  pricePath.end_col(j);
        nrng_.next(itb, ite);
    }
}


#endif // EULERPATHGENERATE_HPP
