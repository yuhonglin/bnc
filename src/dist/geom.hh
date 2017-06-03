#ifndef GEOM_H
#define GEOM_H

#include <dist/pois.hh>

extern "C" {
    double br_dgeom(double x, double p, int give_log);
    double br_pgeom(double x, double p, int lower_tail, int log_p);
    double br_qgeom(double p, double prob, int lower_tail, int log_p);
}

namespace bnc {

    /* R functions */
    
    template<class RNGType>
    double rgeom(double p, RNGType *rng)
    {
	if (!R_FINITE(p) || p <= 0 || p > 1) ML_ERR_return_NAN;
	
	return rpois(exp_rand(rng) * ((1 - p) / p), rng);
    }

    R_RFUNC_INTERFACE_1ARG(rgeom);

    /* D functions */
    R_DFUNC_INTERFACE_3ARG(dgeom, br_dgeom);

    /* P functions */
    R_PFUNC_INTERFACE_4ARG(pgeom, br_pgeom);

    /* Q functions */
    R_QFUNC_INTERFACE_4ARG(qgeom, br_qgeom);
    
    
}  // namespace bnc

#endif /* GEOM_H */
