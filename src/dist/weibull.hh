#ifndef WEIBULL_H
#define WEIBULL_H

#include <dist/brutil.hh>

extern "C" {
    double br_dweibull(double x, double shape,
			double scale, int give_log);
    double br_pweibull(double x, double shape, double scale,
			int lower_tail, int log_p);
    double br_qweibull(double p, double shape, double scale,
			int lower_tail, int log_p);
	
}

namespace bnc {
    
    /* R function */
    
    template<class RNGType>
    double rweibull(double shape, double scale, RNGType *rng)
    {
	if (!R_FINITE(shape) || !R_FINITE(scale) || shape <= 0. || scale <= 0.) {
	    if(scale == 0.) return 0.;
	    /* else */
	    ML_ERR_return_NAN;
	}

	return scale * pow(-log(unif_rand(rng)), 1.0 / shape);
    }

    // Imitate R's rnorm function
    R_RFUNC_INTERFACE_2ARG(rweibull);

    /*
     *  D functions 
     */
    R_DFUNC_INTERFACE_4ARG(dweibull, br_dweibull);

    /* 
     *  P functions
     */
    R_PFUNC_INTERFACE_5ARG(pweibull, br_pweibull);

    /* 
     *  Q functions
     */
    R_QFUNC_INTERFACE_5ARG(qweibull, br_qweibull);

}  // namespace bnc

#endif /* WEIBULL_H */
