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

#define TRUNCATEDNORMAL_IDENTICAL
#include <dist/truncatednormal.hh>

namespace bnc {
    // R function
    template <class RNGType>
    inline Matrix rtmvnorm(const int& n, const Vector& l, const Vector& u,
			   const Matrix& Sig, RNGType* rng) {
	return truncatednormal::rtmvnorm(n, l, u, Sig, ng);
    }

    // D function
    // Notice that the input matrix is both sigma and correlation matrix
    double dtmvnorm(Mvt& mvt, const Vector& x, Vector& mean, Matrix& sigma, Matrix& corr,
		    Vector& lower, Vector& upper, RNGType* rng, const SCALE& s=NORMAL) {
	if (!(corr.diagonal().array()==1).all()) {
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
    
}  // namespace bnc

#endif /* TNORM_H */
