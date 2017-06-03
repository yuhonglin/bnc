#ifndef EXP_H
#define EXP_H

#include <dist/brutil.hh>

extern "C" {
    double br_dexp(double x, double scale, int give_log);
    double br_pexp(double x, double scale, int lower_tail, int log_p);
    double br_qexp(double p, double scale, int lower_tail, int log_p);
}

namespace bnc {
    
    /* R functions */
    
    template<class RNGType>
    double rexp(double scale, RNGType *rng)
    {
	if (!R_FINITE(scale) || scale <= 0.0) {
	    if(scale == 0.) return 0.;
	    /* else */
	    ML_ERR_return_NAN;
	}
	return scale * exp_rand(rng); // --> in ./sexp.c
    }

    R_RFUNC_INTERFACE_1ARG(rexp);

    /* D functions */
    R_DFUNC_INTERFACE_3ARG(dexp, br_dexp);

    /* P functions */
    R_PFUNC_INTERFACE_4ARG(pexp, br_pexp);

    /* Q functions */
    R_QFUNC_INTERFACE_4ARG(qexp, br_qexp);
    
    
}  // namespace bnc

#endif /* EXP_H */
