#ifndef LNORM_H
#define LNORM_H

#include <dist/brutil.hh>

extern "C" {
    double br_dlnorm(double x, double meanlog, double sdlog,
		  int give_log);
    double br_plnorm(double x, double meanlog, double sdlog,
		  int lower_tail, int log_p);
    double br_qlnorm(double p, double meanlog, double sdlog,
		      int lower_tail, int log_p);
}

namespace bnc {

    template<class RNGType>
    double rlnorm(double meanlog, double sdlog, RNGType *rng)
    {
	if(ISNAN(meanlog) || !R_FINITE(sdlog) || sdlog < 0.)
	    ML_ERR_return_NAN;

	return exp(rnorm(meanlog, sdlog, rng));
    }

    // Imitate R's rnorm function
    R_RFUNC_INTERFACE_2ARG(rlnorm);

    /*
     *  D functions 
     */
    R_DFUNC_INTERFACE_4ARG(dlnorm, br_dlnorm);

    /* 
     *  P functions
     */
    R_PFUNC_INTERFACE_5ARG(plnorm, br_plnorm);

    /* 
     *  Q functions
     */
    R_QFUNC_INTERFACE_5ARG(qlnorm, br_qlnorm);
    
    
}  // namespace bnc

#endif /* LNORM_H */
