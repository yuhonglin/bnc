#ifndef DBINOM_H
#define DBINOM_H

#include <limits>
#include <cmath>

#include <rng/rng.hh>
#include <util/logger.hh>
#include <dist/dist.hh>
#include <dist/rmath.hh>
#include <dist/lbeta.hh>

namespace bnc {

    template<SCALE_P scale_p>
    double dbinom_raw(double x, double n, double p, double q)
    {
	double lf, lc;

	if (p == 0) return((x == 0) ? R_D__1 : R_D__0);
	if (q == 0) return((x == n) ? R_D__1 : R_D__0);

	if (x == 0) {
	    if(n == 0) return R_D__1;
	    lc = (p < 0.1) ? -bd0(n,n*q) - n*p : n*log(q);
	    return( R_D_exp(lc) );
	}
	if (x == n) {
	    lc = (q < 0.1) ? -bd0(n,n*p) - n*q : n*log(p);
	    return( R_D_exp(lc) );
	}
	if (x < 0 || x > n) return( R_D__0 );

	/* n*p or n*q can underflow to zero if n and p or q are small.  This
	   used to occur in dbeta, and gives NaN as from R 2.3.0.  */
	lc = stirlerr(n) - stirlerr(x) - stirlerr(n-x) - bd0(x,n*p) - bd0(n-x,n*q);

	/* f = (M_2PI*x*(n-x))/n; could overflow or underflow */
	/* Upto R 2.7.1:
	 * lf = log(M_2PI) + log(x) + log(n-x) - log(n);
	 * -- following is much better for  x << n : */
	lf = M_LN_2PI + log(x) + log1p(- x/n);

	return R_D_exp(lc - 0.5*lf);
    }

    template<SCALE_P scale_p>
    double dbinom(double x, double n, double p)
    {
	/* NaNs propagated correctly */
	if (ISNAN(x) || ISNAN(n) || ISNAN(p)) return x + n + p;

	if (p < 0 || p > 1 || R_D_negInonint(n))
	    ML_ERR_return_NAN;
	R_D_nonint_check(x);
	if (x < 0 || !R_FINITE(x)) return R_D__0;

	n = R_forceint(n);
	x = R_forceint(x);

	return dbinom_raw<scale_p>(x, n, p, 1-p);
    }
    

}  // namespace bnc


#endif /* DBINOM_H */
