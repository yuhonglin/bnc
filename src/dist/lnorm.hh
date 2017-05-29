#ifndef LNORM_H
#define LNORM_H

#include <dist/jrutil.hh>

extern "C" {
    double bnc_dlnorm(double x, double meanlog, double sdlog,
		  int give_log);
    double bnc_plnorm(double x, double meanlog, double sdlog,
		  int lower_tail, int log_p);
    double bnc_qlnorm(double p, double meanlog, double sdlog,
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
    R_DFUNC_INTERFACE_4ARG(dlnorm, bnc_dlnorm);

    /* 
     *  P functions
     */
    R_PFUNC_INTERFACE_5ARG(plnorm, bnc_plnorm);

    /* 
     *  Q functions
     */
    R_QFUNC_INTERFACE_5ARG(qlnorm, bnc_qlnorm);
    
    
}  // namespace bnc

#endif /* LNORM_H */
