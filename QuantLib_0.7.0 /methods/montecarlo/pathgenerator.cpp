#include "pathgenerator.hpp"
#include "../../math/linalg/linalg.hpp"

void PathGenerator::initCorrelation(const Matrix& corrMat)
{
    if (corrMat.is_empty())
        return;     // nothing to do, so sqrtCorrel_ is empty
    Matrix fixedCorrel = corrMat;
    spectrunc(fixedCorrel);                 // spectral truncation
    choldcmp(fixedCorrel, sqrtCorrel_);     // Cholesky decomposition   
}
