#ifndef BETA_H
#define BETA_H

#include <limits>
#include <cmath>

#include <rng/rng.hh>
#include <util/logger.hh>
#include <dist/dist.hh>
#include <dist/rmath.hh>
#include <dist/lbeta.hh>
#include <dist/dbinom.hh>

namespace bnc {

    template<SCALE_P scale_p=NOR_P>
    double dbeta(double x, double a, double b)
    {
	/* NaNs propagated correctly */
	if (ISNAN(x) || ISNAN(a) || ISNAN(b)) return x + a + b;

	if (a < 0 || b < 0) ML_ERR_return_NAN;
	if (x < 0 || x > 1) return(R_D__0);

	// limit cases for (a,b), leading to point masses
	if(a == 0 || b == 0 || !R_FINITE(a) || !R_FINITE(b)) {
	    if(a == 0 && b == 0) { // point mass 1/2 at each of {0,1} :
		if (x == 0 || x == 1) return(ML_POSINF); /* else */ return(R_D__0);
	    }
	    if (a == 0 || a/b == 0) { // point mass 1 at 0
		if (x == 0) return(ML_POSINF); /* else */ return(R_D__0);
	    }
	    if (b == 0 || b/a == 0) { // point mass 1 at 1
		if (x == 1) return(ML_POSINF); /* else */ return(R_D__0);
	    }
	    // else, remaining case:  a = b = Inf : point mass 1 at 1/2
	    if (x == 0.5) return(ML_POSINF); /* else */ return(R_D__0);
	}

	if (x == 0) {
	    if(a > 1) return(R_D__0);
	    if(a < 1) return(ML_POSINF);
	    /* a == 1 : */ return(R_D_val(b));
	}
	if (x == 1) {
	    if(b > 1) return(R_D__0);
	    if(b < 1) return(ML_POSINF);
	    /* b == 1 : */ return(R_D_val(a));
	}

	double lval;
	if (a <= 2 || b <= 2)
	    lval = (a-1)*log(x) + (b-1)*log1p(-x) - lbeta(a, b);
	else
	    lval = log(a+b-1) + dbinom_raw<LOG_P>(a-1, a+b-2, x, 1-x);

	return R_D_exp(lval);
    }


}  // namespace bnc

#endif /* BETA_H */
