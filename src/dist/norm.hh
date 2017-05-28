#ifndef NORM_H
#define NORM_H

#include <dist/jrutil.hh>
#include <dist/dist.hh>

extern "C" {
    double bnc_dnorm4(double x, double mu, double sigma, int give_log);
}

namespace bnc {

    /*
     *   R functions 
     */
    template<class RNGType>
    double rnorm(const double mu, const double sigma, RNGType *rng)
    {
	if (isnan(mu) || !isfinite(sigma) || sigma < 0.)
	    ML_ERR_return_NAN;

	if (sigma == 0. || !isfinite(mu))
	    return mu; /* includes mu = +/- Inf with finite sigma */
	else
	    return mu + sigma * rng->normal();
    }
    
    // Imitate R's rnorm function
    R_RFUNC_INTERFACE_2ARG(rnorm);

    /*
     *  D functions 
     */
    R_DFUNC_INTERFACE_4ARG(dnorm, bnc_dnorm4);
    

}  // namespace bnc

#endif /* NORM_H */
