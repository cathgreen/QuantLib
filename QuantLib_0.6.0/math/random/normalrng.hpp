#ifndef NORMALRNG_HPP
#define NORMALRNG_HPP

#include "../../defines.hpp"
#include "../../exception.hpp"

#include <random>

/** Generator of normal deviates. It is templatized on the underlying uniform RNG
*/
template <typename URNG>
class NormalRng {
public:
    /** Ctor from distribution parameters */
    explicit NormalRng(size_t dimension, double mean = 0, double stdev = 1, const URNG& urng = URNG());

    /** Returns the dimension of the generator */
    size_t dim() const {
        return dim_;
    }

    /** Returns a batch of random deviates using ITER
        CAUTION: it requires end - begin == dimension()
    */
    template <typename ITER>
    void next(ITER begin, ITER end);

    /** Returns the underlying uniform rng. */
    URNG & urng() {
        return urng_;
    }

private:
    // state
    size_t dim_;     // the dimension of the generator
    URNG urng_;      // the uniform random number generator
    normal_distribution<double> normcdf_;   // the normal distribution
    
};


///////////////////////////////////////////////////////////////////////////////
// Inline definitions

template <typename URNG>
NormalRng<URNG>::NormalRng(size_t dimension, double mean, double stdev, const URNG& urng)
    : dim_(dimension),  urng_(urng)
{
    ASSERT(stdev > 0, "NormalRng: the standard deviation must be positive");
    normcdf_ = normal_distribution<double>(mean, stdev);
}


template <typename URNG>
template <typename ITER>
inline void NormalRng<URNG>::next(ITER begin, ITER end) {
    ASSERT((size_t)(end - begin) == dim_, "NormalRng: dimension not match");
    for (ITER it = begin; it != end; ++it)
        *it = normcdf_(urng_);
}



    
#endif // NORMALRNG_HPP
