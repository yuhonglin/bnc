#ifndef GEOM_H
#define GEOM_H

namespace bnc {
    
    
    template<class RNGType>
    double rgeom(double p, RNGType *rng)
    {
	if (!R_FINITE(p) || p <= 0 || p > 1) ML_ERR_return_NAN;
	
	return rpois(exp_rand(rng) * ((1 - p) / p), rng);
    }

}  // namespace bnc

#endif /* GEOM_H */
