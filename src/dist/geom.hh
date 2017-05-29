#ifndef GEOM_H
#define GEOM_H

#include <dist/pois.hh>

extern "C" {
    double bnc_dgeom(double x, double p, int give_log);
    double bnc_pgeom(double x, double p, int lower_tail, int log_p);
    double bnc_qgeom(double p, double prob, int lower_tail, int log_p);
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
    R_DFUNC_INTERFACE_3ARG(dgeom, bnc_dgeom);

    /* P functions */
    R_PFUNC_INTERFACE_4ARG(pgeom, bnc_pgeom);

    /* Q functions */
    R_QFUNC_INTERFACE_4ARG(qgeom, bnc_qgeom);
    
    
}  // namespace bnc

#endif /* GEOM_H */
