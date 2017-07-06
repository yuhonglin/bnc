#ifndef TNORM_H
#define TNORM_H

#include <cmath>
#include <vector>
#include <memory>
#include <numeric>
#include <algorithm>

#include <util/logger.hh>
#include <util/constant.hh>
#include <common/math.hh>
#include <matrix/matrix.hh>
#include <dist/norm.hh>
#include <dist/mvnorm.hh>
#include <dist/chisq.hh>
#include <dist/mvnorm.hh>

#include <dist/truncatednormal.hh>

namespace bnc {
    // R function
    template <class RNGType>
    inline Matrix rtmvnorm(const int& n, const Vector& mean, const Matrix& Sig,
			   const Vector& l, const Vector& u, RNGType* rng) {
	return truncatednormal::rtmvnorm(n, l, u, Sig, rng).array()
	    .colwise() + mean.array();
    }

    // P function
    // notice the input is correlation matrix
    template <class RNGType>    
    double ptmvnorm(Mvt& mvt, Vector& x, Vector& mean,
		    Matrix& corr, Vector& lower, Vector& upper,
		    const SCALE& s, RNGType* rng) {
	if (!((corr.diagonal().array()-1).abs()<1e-15).all()) {
	    // not correlation matrix
	    LOG_WARNING("Input corr matrix' diagonal is not all one");
	}

	if ((x.array()<lower.array()).any()) {
	    if (s == NORMAL)
		return 0.;
	    else
		return NEGINF;
	}

	double l = pmvnorm(mvt, lower, x, mean, corr, rng);
	double u = pmvnorm(mvt, x, upper, mean, corr, rng);
	double p_normal = l/(l+u);
	if (s == NORMAL) {
	    return p_normal;
	} else {
	    return std::log(p_normal);
	}
    }
    template <class RNGType>    
    double ptmvnorm(Vector& x, Vector& mean,
		    Matrix& corr, Vector& lower, Vector& upper,
		    const SCALE& s, RNGType* rng) {
	Mvt mvt;
	return ptmvnorm(mvt, x, mean, corr, lower, upper, s, rng);
    }

    
    // Q function
    // Not sure if this is useful..., not implemented
    
    // D function
    // Notice that the input matrix is both sigma and correlation matrix
    template <class RNGType>    
    double dtmvnorm(Mvt& mvt, const Vector& x, Vector& mean, Matrix& sigma, Matrix& corr,
		    Vector& lower, Vector& upper, const SCALE& s, RNGType* rng) {
	if (!((corr.diagonal().array()-1).abs()<1e-15).all()) {
	    // not correlation matrix
	    LOG_WARNING("Input corr matrix' diagonal is not all one");
	}

	if ((x.array()<lower.array()).any() ||
	    (x.array()>upper.array()).any()) {
	    if (s == NORMAL)
		return 0.;
	    else
		return NEGINF;
	}
	
	if (s == NORMAL) {
	    return dmvnorm(x, mean, sigma, NORMAL) /
		pmvnorm(mvt, lower, upper, mean, corr, rng);
	} else {
	    return dmvnorm(x, mean, sigma, LOG) -
		std::log(pmvnorm(mvt, lower, upper, mean, corr, rng));
	}
    }
    template <class RNGType>    
    double dtmvnorm(const Vector& x, Vector& mean, Matrix& sigma, Matrix& corr,
		    Vector& lower, Vector& upper, const SCALE& s, RNGType* rng) {
	Mvt mvt;
	return dtmvnorm(mvt, x, mean, sigma, corr, lower, upper, s, rng);
    }

    
}  // namespace bnc

#endif /* TNORM_H */
