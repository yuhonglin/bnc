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


// serialisation of matrix and vectors by cereal
namespace cereal
{
  template <class Archive, class _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols> inline
    typename std::enable_if<traits::is_output_serializable<BinaryData<_Scalar>, Archive>::value, void>::type
    save(Archive & ar, Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> const & m)
    {
      int32_t rows = m.rows();
      int32_t cols = m.cols();
      ar(rows);
      ar(cols);
      ar(binary_data(m.data(), rows * cols * sizeof(_Scalar)));
    }

  template <class Archive, class _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols> inline
    typename std::enable_if<traits::is_input_serializable<BinaryData<_Scalar>, Archive>::value, void>::type
    load(Archive & ar, Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> & m)
    {
      int32_t rows;
      int32_t cols;
      ar(rows);
      ar(cols);

      m.resize(rows, cols);

      ar(binary_data(m.data(), static_cast<std::size_t>(rows * cols * sizeof(_Scalar))));
    }
}


#endif /* MATRIX_H */
