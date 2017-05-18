#ifndef LBETA_H
#define LBETA_H

#include <limits>
#include <cmath>

#include <util/logger.hh>
#include <rng/rng.hh>
#include <dist/dist.hh>
#include <dist/rmath.hh>
#include <dist/gammafn.hh>
#include <dist/lgammacor.hh>
#include <dist/lgamma.hh>

namespace bnc {

    double lbeta(double a, double b)
    {
	double corr, p, q;

	if(ISNAN(a) || ISNAN(b))
	    return a + b;

	p = q = a;
	if(b < p) p = b;/* := min(a,b) */
	if(b > q) q = b;/* := max(a,b) */

	/* both arguments must be >= 0 */
	if (p < 0)
	    ML_ERR_return_NAN
	    else if (p == 0) {
		return ML_POSINF;
	    }
	    else if (!R_FINITE(q)) { /* q == +Inf */
		return ML_NEGINF;
	    }

	if (p >= 10) {
	    /* p and q are big. */
	    corr = lgammacor(p) + lgammacor(q) - lgammacor(p + q);
	    return log(q) * -0.5 + M_LN_SQRT_2PI + corr
		+ (p - 0.5) * log(p / (p + q)) + q * log1p(-p / (p + q));
	}
	else if (q >= 10) {
	    /* p is small, but q is big. */
	    corr = lgammacor(q) - lgammacor(p + q);
	    return lgammafn(p) + corr + p - p * log(p + q)
		+ (q - 0.5) * log1p(-p / (p + q));
	}
	else
	    /* p and q are small: p <= q < 10. */
	    /* R change for very small args */
	    if (p < 1e-306) return lgamma(p) + (lgamma(q) - lgamma(p+q));
	return log(gammafn(p) * (gammafn(q) / gammafn(p + q)));
    }
    
}  // namespace bnc

#endif /* LBETA_H */
