#include "linalg.hpp"
#include "../../exception.hpp"

/** 
* Eigenvalues and eigenvectors of a real symmetric matrix
*/


void eigensym(Matrix const& inMat, Vector& eigenValues, Matrix& eigenVectors)
{
    // check square
    ASSERT(inMat.is_square(), "eigensym: input matrix must be square!");
    // check symmetric
    ASSERT(arma::approx_equal(inMat, inMat.t(), "absdiff", 1.0e-16), "eigensym: input matrix must be symmetric!");
    try {
        arma::eig_sym(eigenValues, eigenVectors, inMat);
    }
    catch (...) {
        ASSERT(0, "eigensym: failed to diagonalize the correlation matrix!");
    }
}
