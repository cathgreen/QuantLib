#ifndef PATHGENERATE_HPP
#define PATHGENERATE_HPP

#include "../../defines.hpp"
#include "../../exception.hpp"
#include "../../sptr.hpp"
#include "../../math/matrix.hpp"

class PathGenerator {
public:
    /*
    PathGenerator(size_t ntimesteps, size_t nfactors)
        : ntimesteps_(ntimesteps), nfactors_(nfactors) {}
    */
        
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

    PathGenerator() = default; 
    PathGenerator(size_t ntimesteps, size_t nfactors, const Matrix& correlMatrix);

    void initCorrelation(const Matrix& correlMatrix);

    // state
    size_t ntimesteps_;   
    size_t nfactors_;
    // new added state
    Matrix sqrtCorrel_;    // the Cholesky factor of the correlation matrix
                           // is a lower triangular matrix
};

using SPtrPathGenerator = shared_ptr<PathGenerator>;  // used in BsMcPricers


///////////////////////////////////////////////////////////////////////////////
// Inline definitions

inline
PathGenerator::PathGenerator(size_t ntimesteps, size_t nfactors, const Matrix& correlMatrix)
    : ntimesteps_(ntimesteps), nfactors_(nfactors)
{
    if (!correlMatrix.is_empty()) {
         ASSERT(correlMatrix.is_square(), "PathGenerator: the correlation matrix is not square");
         ASSERT(correlMatrix.n_rows == nfactors,
                "PathGenerator: the correlation matrix number of rows is not equal to nfactors");
    }
    initCorrelation(correlMatrix);    
}



#endif // PATHGENERATE_HPP
