#ifndef WEIBULL2_H
#define WEIBULL2_H

#include <dist/brutil.hh>

extern "C" {
    double br_dweibull2(double x, double shape,
			double scale, int give_log);
    double br_pweibull2(double x, double shape, double scale,
			int lower_tail, int log_p);
    double br_qweibull2(double p, double shape, double scale,
			int lower_tail, int log_p);
	
}

namespace bnc {
    
    /* R function */
    
    template<class RNGType>
    double rweibull2(double shape, double scale, RNGType *rng)
    {
	if (!R_FINITE(shape) || !R_FINITE(scale) || shape <= 0. || scale <= 0.) {
	    if(scale == 0.) return 0.;
	    /* else */
	    ML_ERR_return_NAN;
	}

	return scale * pow(-log(unif_rand(rng)), 1.0 / shape);
    }

    // Imitate R's rnorm function
    R_RFUNC_INTERFACE_2ARG(rweibull2);

    /*
     *  D functions 
     */
    R_DFUNC_INTERFACE_4ARG(dweibull2, br_dweibull2);

    /* 
     *  P functions
     */
    R_PFUNC_INTERFACE_5ARG(pweibull2, br_pweibull2);

    /* 
     *  Q functions
     */
    R_QFUNC_INTERFACE_5ARG(qweibull2, br_qweibull2);

}  // namespace bnc

#endif /* WEIBULL2_H */
