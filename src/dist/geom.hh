#ifndef GEOM_H
#define GEOM_H

#include <dist/pois.hh>

namespace bnc {

    /* R functions */
    
    template<class RNGType>
    double rgeom(double p, RNGType *rng)
    {
	if (!R_FINITE(p) || p <= 0 || p > 1) ML_ERR_return_NAN;
	
	return rpois(exp_rand(rng) * ((1 - p) / p), rng);
    }

    R_RFUNC_INTERFACE_1ARG(rgeom);
    
}  // namespace bnc

#endif /* GEOM_H */
