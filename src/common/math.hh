#ifndef COMMON_MATH_H
#define COMMON_MATH_H

#include <cmath>

#include "matrix/matrix.hh"
#include "util/logger.hh"

namespace bnc {

    // matrix decomposition methods
    enum MAT_DECOMP{
	EIGEN_DECOMP,    // eigen (spectral) decomposition
	CHOL_DECOMP,     // Cholesky decomposition
	RCHOL_DECOMP     // robust Cholesky decomposition
    };

    enum MAT_TYPE {
	MAT_SPD,
	MAT_PD
    };

    enum MAT_STRUCTURE {
	UPPER_TRIANGLE,
	LOWER_TRIANGLE,
	MAT_FULL
    };
    
    inline double fmin2(const double& x, const double& y) {
	if (std::isnan(x) || std::isnan(y))
	    return x + y;
	return (x < y) ? x : y;
    }

    inline double fmax2(const double& x, const double& y)
    {
	if (std::isnan(x) || std::isnan(y))
	    return x + y;
	return (x < y) ? y : x;
    }

    inline int imin2(const int &x, const int &y)
    {
	return (x < y) ? x : y;
    }

    inline int imax2(const int &x, const int &y)
    {
	return (x < y) ? y : x;
    }

    inline double fsign(const double &x, const double &y)
    {
	if (std::isnan(x) || std::isnan(y))
	    return x + y;
	return ((y >= 0) ? std::fabs(x) : -std::fabs(x));
    }

    /**
     * Convert a triangular matrix to a symmetric matrix
     * by filling the other half.
     */
    template<MAT_STRUCTURE instructure, class T>
    inline void tri2sym_inplace(T& A) {
	if (instructure==UPPER_TRIANGLE) {
	    int idx = 1;
	    int idx2;
	    for (int j=0; j<A.cols(); j++) {       // col
		idx2 = (j+1)*A.cols() + j;
		for (int i=j+1; i<A.cols(); i++) { // row
		    A(idx) = A(idx2);
		    idx++;
		    idx2 += A.cols();
		}
		idx += j;
		idx ++;
		idx ++;		
	    }
	    return;
	} else if (instructure==LOWER_TRIANGLE) {
	    int idx;
	    int idx2 = 1;
	    for (int j=0; j<A.cols(); j++) {       // col
		idx = (j+1)*A.cols() + j;
		for (int i=j+1; i<A.cols(); i++) { // row
		    A(idx) = A(idx2);
		    idx2++;
		    idx += A.cols();
		}
		idx2 += j;
		idx2 ++;
		idx2 ++;
	    }
	    return;
	}

	LOG_ERROR("Matrix structure not supported.");
    }
    
    /**
     * Compute A'*A (Matrix dot product)
     * 
     * Naive way in Eigen is to 
     * 
     *            Matrix result = A.transpose()*A;
     * 
     * but since the result is always symmetric, it is not
     * necessary to times the other half. 
     * 
     * @return: a matrix with meanfully upper/lower triangular
     *          part. The other half is undefined.
     *
     */
    template <MAT_STRUCTURE ms>
    inline Matrix mdot(const Matrix& A) {
	const double * a = A.data();
	Matrix ret = Matrix(A.cols(),A.cols());
	double * r = ret.data();
	int idx = 0;
	if (ms == LOWER_TRIANGLE) {
	    // lower triangular
	    for (int j=0; j<A.cols(); j++) {     // col
		for (int i=j; i<A.cols(); i++) { // row
		    // update r[idx] (ret(i,j))
		    r[idx] = 0;
		    int k1 = i*A.rows();
		    int k2 = j*A.rows();
		    for (int k=0; k<A.rows(); k++) {
			r[idx] += A(k1)*A(k2);
			k1++;
			k2++;
		    }
		    idx++;
		}
		idx += j;
		idx ++;
	    }
	} else if (ms == UPPER_TRIANGLE) {
            // upper triangular
	    for (int j=0; j<A.cols(); j++) {     // col
		idx -= j;
		for (int i=0; i<=j; i++) { // row
		    // update r[idx] (ret(i,j))
		    r[idx] = 0;
		    int k1 = i*A.rows();
		    int k2 = j*A.rows();
		    for (int k=0; k<A.rows(); k++) {
			r[idx] += A(k1)*A(k2);
			k1++;
			k2++;
		    }
		    idx++;
		}
		idx += A.cols();
	    }
	} else if (ms == MAT_FULL) {
	    ret = mdot<UPPER_TRIANGLE>(A);
	    ret.triangularView<Eigen::Lower>() = ret.transpose(); // this may unnecessarily copy the diagnal
	    //tri2sym_inplace<UPPER_TRIANGLE>(ret);
	} else {
	    LOG_ERROR("Matrix structure not supported.");
	}
	return ret;
    }

    // // imitate the cov2cor function in R
    // Matrix cov2cor(const Matrix &m) {
    // 	Vector is = 1./m.diagonal().array().sqrt();
    // 	return (m.array().rowwise()*is.transpose().array()).array().colwise()*is.array();
    // }

    // imitate R's seq function.
    
} // namespace bnc

#endif /* COMMON_MATH_H */
