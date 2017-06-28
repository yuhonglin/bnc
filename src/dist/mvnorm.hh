/**
  Partly adapted from jags.
 */

#ifndef MVNORM_H
#define MVNORM_H

#include <limits>

#include <Eigen/Eigenvalues>
#include <util/logger.hh>
#include <util/constant.hh>
#include <common/math.hh>
#include <matrix/matrix.hh>
#include <dist/norm.hh>
#include <dist/dist.hh>
#include <dist/mvt.hh>

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
	    Eigen::SelfAdjointEigenSolver<Matrix> es(sigma);
	    Vector ret;
	    if (prec==PRECISION) {
		// sigma is precision matrix
		ret = rnorm(mu.size(), 0., es.eigenvalues().
			    array().sqrt().inverse(), rng);
	    } else {
		// sigma is covariance matrix
		ret = rnorm(mu.size(), 0., es.eigenvalues().array().
			    sqrt(), rng);
	    }
	    return es.eigenvectors()*ret + mu;
	}

	LOG_ERROR("Matrix decomposition method not supported, return NaN.");
	return Vector::Constant(mu.size(),
				NaN);
    }

    // generate one sample but use eigen decomposition as input
    template<MVNORM_INPUT prec=VARIANCE,
	     MAT_DECOMP MD=EIGEN_DECOMP, class RNGType>
    Vector rmvnorm(const Vector& mu, const Eigen::SelfAdjointEigenSolver<Matrix>& es,
		   RNGType* rng)
    {
	Vector ret;
	if (prec==PRECISION) {
	    // sigma is precision matrix
	    ret = rnorm(mu.size(), 0., es.eigenvalues().
			array().sqrt().inverse(), rng);
	} else {
	    // sigma is covariance matrix
	    ret = rnorm(mu.size(), 0., es.eigenvalues().array().
			sqrt(), rng);
	}
	return es.eigenvectors()*ret + mu;
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
	Eigen::SelfAdjointEigenSolver<Matrix> es(sigma);
	
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

    // P functions, based on mvt.hh
    // Notice that the input matrix is correlations matrix
    template <class RNGType>
    double pmvnorm(Vector &lower, Vector &upper, Vector &mean,
		   Matrix &corr, RNGType *rng, double& error, int &inform) {

	if (!(corr.diagonal().array()==1).all()) {
	    // not correlation matrix
	    LOG_WARNING("Input corr matrix' diagonal is not all one");
	}
	
	Vector l = lower - mean;
	Vector u = upper - mean;
	Vector m = Vector::Constant(l.size(), 0.);
	
	int df = 0;
	return mvt::mvt(l, u, df, corr, m, error, inform, rng);
    }
    
}  // namespace bnc

#endif /* MVNORM_H */
