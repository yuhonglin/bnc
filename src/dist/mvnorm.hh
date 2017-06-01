/**
  Partly adapted from jags.
 */

#ifndef MVNORM_H
#define MVNORM_H

#include <limits>

#include <Eigen/Eigenvalues>
#include <util/logger.hh>
#include <common/math.hh>
#include <matrix/matrix.hh>
#include <dist/norm.hh>
#include <dist/dist.hh>


namespace bnc {

    enum MVNORM_INPUT { PRECISION, VARIANCE };
    
    /* R function */

    // generate only 1 sample
    template<MVNORM_INPUT prec=VARIANCE,
	     MAT_DECOMP MD=EIGEN_DECOMP, class RNGType>
    Vector rmvnorm(const Vector& mu, const Matrix& sigma,
		   RNGType* rng)
    {
	if (MD==CHOL_DECOMP) {
	    // use cholesky decomposition (fastest, but only for PD sigma)
	    Vector ret = rnorm(mu.size(), 0., 1., rng);
	    Matrix L = sigma.llt().matrixL();
	    
	    if (prec==PRECISION) {
		// sigma is precision matrix
		return L.triangularView<Eigen::Lower>().solve(ret) + mu;
	    } else {
		// sigma is covariance matrix
		return L.triangularView<Eigen::Lower>()*ret + mu;
	    }
	}

	if (MD==RCHOL_DECOMP) {
	    // use robest cholesky decomposition
	    Vector ret = rnorm(mu.size(), 0., 1., rng);
	    Matrix L = sigma.ldlt().matrixL();
	    
	    if (prec==PRECISION) {
		// sigma is precision matrix
		return L.triangularView<Eigen::Lower>().solve(ret) + mu;
	    } else {
		// sigma is covariance matrix
		return L.triangularView<Eigen::Lower>()*ret + mu;
	    }
	}

	if (MD==EIGEN_DECOMP) {
	    // use eigen decomposition
	    Eigen::EigenSolver<Matrix> es(sigma);
	    Vector ret;
	    if (prec==PRECISION) {
		// sigma is precision matrix
		ret = rnorm(mu.size(), 0., es.eigenvalues().real().
			    array().sqrt().inverse(), rng);
	    } else {
		// sigma is covariance matrix
		ret = rnorm(mu.size(), 0., es.eigenvalues().array().
			    real().sqrt(), rng);
	    }
	    return es.eigenvectors().real()*ret + mu;
	}

	LOG_ERROR("Matrix decomposition method not supported, return NaN.");
	return Vector::Constant(mu.size(),
				std::numeric_limits<double>::quiet_NaN());
    }

    // generate one sample but use eigen decomposition as input
    template<MVNORM_INPUT prec=VARIANCE,
	     MAT_DECOMP MD=EIGEN_DECOMP, class RNGType>
    Vector rmvnorm(const Vector& mu, const Eigen::EigenSolver<Matrix>& es,
		   RNGType* rng)
    {
	Vector ret;
	if (prec==PRECISION) {
	    // sigma is precision matrix
	    ret = rnorm(mu.size(), 0., es.eigenvalues().real().
			array().sqrt().inverse(), rng);
	} else {
	    // sigma is covariance matrix
	    ret = rnorm(mu.size(), 0., es.eigenvalues().array().
			real().sqrt(), rng);
	}
	return es.eigenvectors().real()*ret + mu;
    }

    // generate multiple samples
    template<MVNORM_INPUT prec=VARIANCE,
	     MAT_DECOMP MD=EIGEN_DECOMP, class RNGType>
    Matrix rmvnorm(const int& n, const Vector& mu,
		   const Matrix& sigma, RNGType* rng)
    {
	static_assert(MD==EIGEN_DECOMP || CHOL_DECOMP,
		      "Matrix factorisation method not supported");

	Matrix ret(n, mu.size());
	Eigen::EigenSolver<Matrix> es(sigma);
	
	for (int i=0; i<n; i++) {
	    ret.row(i) = rmvnorm<prec,MD,RNGType>(mu, es, rng);
	}

	return ret;
    }


    /* D functions */

    // for row vector
    template <int R, int C>
    double dmvnorm(const Eigen::Matrix<double,R,C>& x,
		   const typename enable_if<R==1,Vector>::type& mu,
		   const Matrix& sigma, const SCALE& s=NORMAL) {

	auto diff = x.transpose() - mu;
	double logden = -diff.transpose()*sigma.inverse()
	    .dot(diff.transpose())/2;

	if (s == LOG) {
	    logden -= log(sigma.determinant())/2 +
		mu.size()* 0.572364942924700087071713675677;
	    return logden;
	} else {
	    return exp(logden)*
		pow(0.398942280401432677939946059934,mu.size())/
		sqrt(sigma.determinant());
	}
    }

    // for column vector
    template <int R, int C>
    double dmvnorm(const Eigen::Matrix<double,R,C>& x,
		   const typename enable_if<C==1,Vector>::type& mu,
		   const Matrix& sigma, const SCALE& s=NORMAL) {

	auto diff = x - mu;
	double logden = -(diff.transpose()*sigma.inverse())
	    .dot(diff.transpose())/2;

	if (s == LOG) {
	    logden -= log(sigma.determinant())/2 +
		mu.size()*0.918938533204672741780329736406;
	    return logden;
	} else {
	    return exp(logden)*
		pow(0.398942280401432677939946059934,mu.size())/
		sqrt(sigma.determinant());
	}
	
    }

    
}  // namespace bnc

#endif /* MVNORM_H */
