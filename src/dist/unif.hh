#ifndef UNIF_H
#define UNIF_H

#include <dist/brutil.hh>

extern "C" {
    double br_dunif(double x, double a, double b, int give_log);
    double br_punif(double x, double a, double b, int lower_tail, int log_p);
    double br_qunif(double p, double a, double b, int lower_tail, int log_p);
}

namespace bnc {

    /* R functions */
    
    template<class RNGType>
    double runif(const double& a, const double& b, RNGType *rng)
    {
	if (!std::isfinite(a) || !std::isfinite(b) || b < a)
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
    R_DFUNC_INTERFACE_4ARG(dunif, br_dunif);

    /* 
     *  P functions
     */
    R_PFUNC_INTERFACE_5ARG(punif, br_punif);

    /* 
     *  Q functions
     */
    R_QFUNC_INTERFACE_5ARG(qunif, br_qunif);
    

}  // namespace bnc

#endif /* UNIF_H */
