#ifndef CAUCHY_H
#define CAUCHY_H

#include <dist/jrutil.hh>

namespace bnc {
    template<class RNGType>
    double rcauchy(const double &location, const double &scale, RNGType *rng)
    {
	if (ISNAN(location) || !R_FINITE(scale) || scale < 0)
	    ML_ERR_return_NAN;
	if (scale == 0. || !R_FINITE(location))
	    return location;
	else
	    return location + scale * tan(M_PI * unif_rand(rng));
    }

    // imitate R interfaces
    R_RFUNC_INTERFACE_2ARG(rcauchy);
    
}  // namespace bnc

#endif /* CAUCHY_H */
