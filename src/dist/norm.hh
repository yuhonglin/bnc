#ifndef NORM_H
#define NORM_H

#include <type_traits>
#include <limits>
#include <cmath>
#include <cfloat>

#include <dist/jrutil.hh>

namespace bnc {

    /* R functions */
    
    template<class RNGType>
    double rnorm(const double mu, const double sigma, RNGType *rng)
    {
	if (isnan(mu) || !isfinite(sigma) || sigma < 0.)
	{
	    LOG_ERROR("invalid input, return NaN");
	    return numeric_limits<double>::quiet_NaN();
	}

	if (sigma == 0. || !isfinite(mu))
	    return mu; /* includes mu = +/- Inf with finite sigma */
	else
	    return mu + sigma * rng->normal();
    }
    
    // Imitate R's rnorm function
    R_RFUNC_INTERFACE_2ARG(rnorm);

}  // namespace bnc

#endif /* NORM_H */
