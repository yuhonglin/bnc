#ifndef NBINOM_H
#define NBINOM_H

#include <dist/jrmath.hh>

extern "C" {
    double bnc_dnbinom(double x, double size, double prob,
		       int give_log);
    double bnc_pnbinom(double x, double size, double prob,
		       int lower_tail, int log_p);
    double bnc_qnbinom(double p, double size, double prob,
		       int lower_tail, int log_p)
}

namespace bnc {
    /* R function */
    
    template<class RNGType>
    double rnbinom(double size, double prob, RNGType *rng)
    {
	if(!R_FINITE(size) || !R_FINITE(prob) || size <= 0 || prob <= 0 || prob > 1)
	    /* prob = 1 is ok, PR#1218 */
	    ML_ERR_return_NAN;
	return (prob == 1) ? 0 : rpois(rgamma(size, (1 - prob) / prob, rng), rng);
    }

    template<class RNGType>
    double rnbinom_mu(double size, double mu, RNGType *rng)
    {
	if(!R_FINITE(size) || !R_FINITE(mu) || size <= 0 || mu < 0)
	    ML_ERR_return_NAN;
	return (mu == 0) ? 0 : rpois(rgamma(size, mu / size, rng), rng);
    }

    // Imitate R's rnorm function
    R_RFUNC_INTERFACE_2ARG(rnorm);

    /*
     *  D functions 
     */
    R_DFUNC_INTERFACE_4ARG(dnbinom, bnc_dnbinom);

    /* 
     *  P functions
     */
    R_PFUNC_INTERFACE_5ARG(pnbinom, bnc_pnbinom);

    /* 
     *  Q functions
     */
    R_QFUNC_INTERFACE_5ARG(qnbinom, bnc_qnbinom);

}  // namespace bnc

#endif /* NBINOM_H */
