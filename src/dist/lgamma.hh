/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 1998 Ross Ihaka
 *  Copyright (C) 2000-2012 The R Core Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  http://www.r-project.org/Licenses/
 *
 *  SYNOPSIS
 *
 *    #include <Rmath.h>
 *    double lgammafn_sign(double x, int *sgn);
 *    double lgammafn(double x);
 *
 *  DESCRIPTION
 *
 *    The function lgammafn computes log|gamma(x)|.  The function
 *    lgammafn_sign in addition assigns the sign of the gamma function
 *    to the address in the second argument if this is not NULL.
 *
 *  NOTES
 *
 *    This routine is a translation into C of a Fortran subroutine
 *    by W. Fullerton of Los Alamos Scientific Laboratory.
 *
 *    The accuracy of this routine compares (very) favourably
 *    with those of the Sun Microsystems portable mathematical
 *    library.
 */

#ifndef LGAMMA_H
#define LGAMMA_H

#include <limits>
#include <cmath>

using namespace std;

#include <rng/rng.hh>
#include <util/logger.hh>

#include <dist/dist.hh>
#include <dist/rmath.hh>
#include <dist/cospi.hh>
#include <dist/lgammacor.hh>

namespace bnc {
    
    double gammafn(double x);

    double lgammafn_sign(const double& x, int *sgn)
    {
	double ans, y, sinpiy;

/* For IEEE double precision DBL_EPSILON = 2^-52 = 2.220446049250313e-16 :
   xmax  = DBL_MAX / log(DBL_MAX) = 2^1024 / (1024 * log(2)) = 2^1014 / log(2)
   dxrel = sqrt(DBL_EPSILON) = 2^-26 = 5^26 * 1e-26 (is *exact* below !)
*/
#undef xmax	
#define xmax  2.5327372760800758e+305
#undef dxrel
#define dxrel 1.490116119384765625e-8

	if (sgn != NULL) *sgn = 1;

	if(isnan(x)) return x;

	if (sgn != NULL && x < 0 && fmod(floor(-x), 2.) == 0)
	    *sgn = -1;

	if (x <= 0 && x == trunc(x)) { /* Negative integer argument */
	    LOG_ERROR("invalid x");
	    return ML_POSINF;/* +Inf, since lgamma(x) = log|gamma(x)| */
	}

	y = fabs(x);

	if (y < 1e-306) return -log(y); // denormalized range, R change
	if (y <= 10) return log(fabs(gammafn(x)));
	/*
	  ELSE  y = |x| > 10 ---------------------- */

	if (y > xmax) {
	    LOG_ERROR("y>xmax, return inf");
	    return ML_POSINF;
	}

	if (x > 0) { /* i.e. y = x > 10 */
	    if(x > 1e17)
		return(x*(log(x) - 1.));
	    else if(x > 4934720.)
		return(M_LN_SQRT_2PI + (x - 0.5) * log(x) - x);
	    else
		return M_LN_SQRT_2PI + (x - 0.5) * log(x) - x + lgammacor(x);
	}
	/* else: x < -10; y = -x */
	sinpiy = fabs(sinpi(y));

	if (sinpiy == 0) { /* Negative integer argument ===
			      Now UNNECESSARY: caught above */
	    LOG_WARNING(" ** should NEVER happen! ** ");
	    ML_ERR_return_NAN;
	}

	ans = M_LN_SQRT_PId2 + (x - 0.5) * log(y) - x - log(sinpiy) - lgammacor(y);

	if(fabs((x - trunc(x - 0.5)) * ans / x) < dxrel) {

	    /* The answer is less than half precision because
	     * the argument is too near a negative integer. */

	    LOG_ERROR("not precise");
	}

	return ans;
    }

    double lgammafn(const double& x)
    {
	return lgammafn_sign(x, NULL);
    }

}

#endif /* LGAMMA_H */
