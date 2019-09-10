#ifndef LINALG_HPP
#define LINALG_HPP

#include "../matrix.hpp"

/** 
 * Cholesky decomposition of a positive semi-definite matrix inMat.
 * It computes the lower triangular part of outMat such that outMat * trans(outMat) = inMat.
 */
void choldcmp(Matrix const& inMat, Matrix& outMat);

/** 
 * Eigenvalues and eigenvectors of a real symmetric matrix
 */
void eigensym(Matrix const& inputMatrix, Vector& eigenValues, Matrix& eigenVectors);

/** 
 * Spectral truncation of the input correlation matrix.
 * The input matrix must be symmetric with ones along the diagonal.
 * Spectral truncation happens in place and the returned matrix is symmetric, 
 * positive semi-definite and with ones along the diagonal.
 */
void spectrunc(Matrix& corrmat, double tolerance = 1e-8);

#endif // LINALG_HPP
