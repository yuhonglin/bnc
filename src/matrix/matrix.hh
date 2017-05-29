/**
 * Currently, the matrix system is based on eigen.
 * 
 */

#ifndef MATRIX_H
#define MATRIX_H

#include <Eigen/Dense>

namespace bnc {
    const int Dynamic = Eigen::Dynamic;

    // matrix
    template<int m=Dynamic, int n=Dynamic, class DataType=double>
    using MatrixTemplate = Eigen::Matrix<DataType, m, n>;

    using Matrix = MatrixTemplate<>;

    // vector
    template<int m=Dynamic, class DataType=double>
    using VectorTemplate = Eigen::Matrix<DataType, m, 1>;

    using Vector = VectorTemplate<>;

    
}  // namespace bnc

#endif /* MATRIX_H */
