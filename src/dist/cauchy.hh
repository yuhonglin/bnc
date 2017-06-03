#ifndef CAUCHY_H
#define CAUCHY_H

#include <dist/brutil.hh>

extern "C" {
    double br_dcauchy(double x, double location, double scale, int give_log);
    double br_pcauchy(double x, double location, double scale,
		       int lower_tail, int log_p);
    double br_qcauchy(double p, double location, double scale,
		       int lower_tail, int log_p);
}

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

    /*
     *  D functions 
     */
    R_DFUNC_INTERFACE_4ARG(dcauchy, br_dcauchy);

    /* 
     *  P functions
     */
    R_PFUNC_INTERFACE_5ARG(pcauchy, br_pcauchy);

    /* 
     *  Q functions
     */
    R_QFUNC_INTERFACE_5ARG(qcauchy, br_qcauchy);
    
    
}  // namespace bnc

#endif /* CAUCHY_H */
