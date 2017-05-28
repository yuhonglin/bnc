#ifndef UNIF_H
#define UNIF_H

#include <dist/jrutil.hh>

extern "C" {
    double bnc_dunif(double x, double a, double b, int give_log);
    double bnc_punif(double x, double a, double b, int lower_tail, int log_p);
    double bnc_qunif(double p, double a, double b, int lower_tail, int log_p);
}

namespace bnc {

    /* R functions */
    
    template<class RNGType>
    double runif(const double& a, const double& b, RNGType *rng)
    {
	if (!isfinite(a) || !isfinite(b) || b < a)
	    ML_ERR_return_NAN;

	if (a == b)
	    return a;
	else {
	    double u;
	    /* This is true of all builtin generators, but protect against
	       user-supplied ones */
	    do { u = rng->uniform(); } while (u <= 0 || u >= 1);
	    return a + (b - a) * u;
	}
    }

    // Imitate R's rnorm function
    R_RFUNC_INTERFACE_2ARG(runif);

    /*
     *  D functions 
     */
    R_DFUNC_INTERFACE_4ARG(dunif, bnc_dunif);

    /* 
     *  P functions
     */
    R_PFUNC_INTERFACE_5ARG(punif, bnc_punif);

    /* 
     *  Q functions
     */
    R_QFUNC_INTERFACE_5ARG(qunif, bnc_qunif);
    

}  // namespace bnc

#endif /* UNIF_H */
