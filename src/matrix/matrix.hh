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
    template<class DataType=double, int m=Dynamic, int n=Dynamic>
    using MatrixTemplate = Eigen::Matrix<DataType, m, n>;

    using Matrix = MatrixTemplate<>;
    using dMatrix = MatrixTemplate<>;
    using iMatrix = MatrixTemplate<int>;

    // vector
    template<class DataType=double, int m=Dynamic>
    using VectorTemplate = Eigen::Matrix<DataType, m, 1>;

    using Vector = VectorTemplate<>;
    using dVector = VectorTemplate<>;
    using iVector = VectorTemplate<int>;

    
}  // namespace bnc

#endif /* MATRIX_H */
