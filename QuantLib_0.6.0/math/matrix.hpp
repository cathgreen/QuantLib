#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <armadillo>
#include "../defines.hpp"


/** The Vector class is an alias for the armadillo column vector, a sequence of doubles */
using Vector = arma::vec;

using RowVector = arma::rowvec;

/** The Matrix class is an alias for the armadillo matrix, a dense matrix of doubles.
    The storage is column-wise. Access to the underlying data can be gained via the .memptr() method.
*/
using Matrix = arma::mat;


#endif // MATRIX_HPP
