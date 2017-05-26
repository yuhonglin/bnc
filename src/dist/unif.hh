#ifndef UNIF_H
#define UNIF_H

#include <dist/jrutil.hh>

namespace bnc {

    /* R functions */
    
    template<class RNGType>
    double runif(const double& a, const double& b, RNGType *rng)
    {
	if (!isfinite(a) || !isfinite(b) || b < a)
	    ML_ERR_return_NAN;

	if (a == b)
	    return a;
	else {
	    double u;
	    /* This is true of all builtin generators, but protect against
	       user-supplied ones */
	    do { u = rng->uniform(); } while (u <= 0 || u >= 1);
	    return a + (b - a) * u;
	}
    }

    // Imitate R's rnorm function
    R_RFUNC_INTERFACE_2ARG(runif);

}  // namespace bnc

#endif /* UNIF_H */
