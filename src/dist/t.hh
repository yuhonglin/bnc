#ifndef T_H
#define T_H

#include <dist/jrutil.hh>

extern "C" {
    double bnc_dt(double x, double n, int give_log);
    double bnc_pt(double x, double n, int lower_tail, int log_p);
    double bnc_qt(double p, double ndf, int lower_tail, int log_p);
}

namespace bnc {
    /* R function */
    
    template<class RNGType>
    double rt(double df, RNGType *rng)
    {
	if (ISNAN(df) || df <= 0.0)	ML_ERR_return_NAN;

	if(!R_FINITE(df))
	    return norm_rand(rng);
	else {
/* Some compilers (including MW6) evaluated this from right to left
   return norm_rand(rng) / sqrt(rchisq(df, rng) / df); */
	    double num = norm_rand(rng);
	    return num / sqrt(rchisq(df, rng) / df);
	}
    }

    R_RFUNC_INTERFACE_1ARG(rt);

    /* D functions */
    R_DFUNC_INTERFACE_3ARG(dt, bnc_dt);

    /* P functions */
    R_PFUNC_INTERFACE_4ARG(pt, bnc_pt);

    /* Q functions */
    R_QFUNC_INTERFACE_4ARG(qt, bnc_qt);
    

}  // namespace bnc

#endif /* T_H */
