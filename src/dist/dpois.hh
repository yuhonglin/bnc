#ifndef DPOIS_H
#define DPOIS_H

#include <limits>
#include <cmath>

using namespace std;

#include <rng/rng.hh>
#include <util/logger.hh>
#include <dist/dist.hh>
#include <dist/rmath.hh>
#include <dist/bd0.hh>
#include <dist/stirlerr.hh>
#include <dist/lgamma.hh>

namespace bnc {

    template<SCALE_P scale_p=NOR_P>
    double dpois_raw(double x, double lambda)
    {
	/*       x >= 0 ; integer for dpois(), but not e.g. for pgamma()!
		 lambda >= 0
	*/
	if (lambda == 0) return( (x == 0) ? R_D__1 : R_D__0 );
	if (isinf(lambda)) return R_D__0;
	if (x < 0) return( R_D__0 );
	if (x <= lambda * DBL_MIN) return(R_D_exp(-lambda) );
	if (lambda < x * DBL_MIN) return(R_D_exp(-lambda + x*log(lambda) -lgammafn(x+1)));
	return(R_D_fexp( M_2PI*x, -stirlerr(x)-bd0(x,lambda) ));
    }

}  // bnc


#endif /* DPOIS_H */
