#include "linalg.hpp"
#include "../../exception.hpp"

/** Cholesky decomposition of a positive semi-definite matrix inMat.
    The returned matrix is lower triangular, such that 
    outMat * outMat^T = inMat
*/

void choldcmp(Matrix const& inMat, Matrix& outMat)
{
    // check square
    ASSERT(inMat.is_square(), "choldcmp: input matrix must be square!");
    // check symmetric
    ASSERT(arma::approx_equal(inMat, inMat.t(), "absdiff", 1.0e-16), "choldcmp: input matrix must be symmetric!");
    
    bool  ok = arma::chol(outMat, inMat, "lower");
    ASSERT(ok, "choldcmp: input matrix not positive definite!");
    return;
}
