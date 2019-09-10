#ifndef PATHGENERATE_HPP
#define PATHGENERATE_HPP

#include "../../defines.hpp"
#include "../../exception.hpp"
#include "../../sptr.hpp"

class PathGenerator {
public:
    PathGenerator(size_t ntimesteps, size_t nfactors)
        : ntimesteps_(ntimesteps), nfactors_(nfactors) {}
    
    virtual ~PathGenerator() {}

    size_t nTimeSteps() const {
        return ntimesteps_;
    }

    size_t nFactors() const {
        return nfactors_;
    }

    /** Write ont he next price path on pricePath
        = ntimesteps_ * nfactors_
    */
    virtual void next(Matrix& pricePath) = 0;
    
protected:

    // state
    size_t ntimesteps_;   
    size_t nfactors_;
};

using SPtrPathGenerator = shared_ptr<PathGenerator>;  // used in BsMcPricers

#endif // PATHGENERATE_HPP
