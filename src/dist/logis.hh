#ifndef LOGIS_H
#define LOGIS_H

#include <dist/jrutil.hh>

extern "C" {
    double bnc_dlogis(double x, double location, double scale, int give_log);
    double bnc_plogis(double x, double location, double scale,
		      int lower_tail, int log_p);
    double bnc_qlogis(double p, double location, double scale,
		      int lower_tail, int log_p);
}


namespace bnc {
    
    /* R function */
    
    template<class RNGType>
    double rlogis(double location, double scale, RNGType *rng)
    {
	if (ISNAN(location) || !R_FINITE(scale))
	    ML_ERR_return_NAN;

	if (scale == 0. || !R_FINITE(location))
	    return location;
	else {
	    double u = unif_rand(rng);
	    return location + scale * log(u / (1. - u));
	}
    }

    // Imitate R's rnorm function
    R_RFUNC_INTERFACE_2ARG(rlogis);

    /*
     *  D functions 
     */
    R_DFUNC_INTERFACE_4ARG(dlogis, bnc_dlogis);

    /* 
     *  P functions
     */
    R_PFUNC_INTERFACE_5ARG(plogis, bnc_plogis);

    /* 
     *  Q functions
     */
    R_QFUNC_INTERFACE_5ARG(qlogis, bnc_qlogis);

    
}  // namespace bnc

#endif /* LOGIS_H */
