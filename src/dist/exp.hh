#ifndef EXP_H
#define EXP_H

#include <dist/exp.hh>

namespace bnc {
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

}  // namespace bnc

#endif /* EXP_H */
