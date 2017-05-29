/**
  Partly adapted from jags.
 */

#ifndef MVNORM_H
#define MVNORM_H

#include <Eigen/Eigenvalues>
#include <util/logger.hh>
#include <common/math.hh>
#include <matrix/matrix.hh>
#include <dist/norm.hh>


namespace bnc {

    enum MVNORM_INPUT { PRECISION, VARIANCE };
    
    /* R function */

    // generate only 1 sample
    template<class RNGType, MVNORM_INPUT prec=VARIANCE,
	     MAT_DECOMP MD=EIGEN_DECOMP>
    Vector mvrnorm(const Vector& mu, const Matrix& sigma,
		   RNGType* rng)
    {
	static_assert(MD==EIGEN_DECOMP || CHOL_DECOMP,
		      "Matrix factorisation method not supported");

	if (MD!=EIGEN_DECOMP) {
	    TODO
	}
	
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
	return es.eigenvectors().real()*ret;
    }

    // generate one sample but use eigen decomposition as input
    template<class RNGType, MVNORM_INPUT prec=VARIANCE,
	     MAT_DECOMP MD=EIGEN_DECOMP>
    Vector mvrnorm(const Vector& mu, const Eigen::EigenSolver<Matrix>& es,
		   RNGType* rng)
    {
	static_assert(MD==EIGEN_DECOMP || CHOL_DECOMP,
		      "Matrix factorisation method not supported");

	if (MD!=EIGEN_DECOMP) {
	    TODO
	}
	
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
	return es.eigenvectors().real()*ret;
    }

    // generate multiple samples
    template<class RNGType, MVNORM_INPUT prec=VARIANCE,
	     MAT_DECOMP MD=EIGEN_DECOMP>
    Matrix mvrnorm(const int& n, const Vector& mu,
		   const Matrix& sigma, RNGType* rng)
    {
	static_assert(MD==EIGEN_DECOMP || CHOL_DECOMP,
		      "Matrix factorisation method not supported");

	Matrix ret(n, mu.size());
	Eigen::EigenSolver<Matrix> es(sigma);
	
	for (int i=0; i<n; i++) {
	    ret.row(i) = mvrnorm<RNGType,prec,MD>(mu, es, rng);
	}

	return ret;
    }

    
}  // namespace bnc

#endif /* MVNORM_H */
