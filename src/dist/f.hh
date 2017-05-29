#ifndef F_H
#define F_H

extern "C" {
    double bnc_dF(double x, double m, double n, int give_log);
    double bnc_pF(double x, double df1, double df2, int lower_tail, int log_p);
    double bnc_qF(double p, double df1, double df2, int lower_tail, int log_p);
}

namespace bnc {

    /* R functions */
    
    template<class RNGType>
    double rf(double n1, double n2, RNGType *rng)
    {
	double v1, v2;
	if (ISNAN(n1) || ISNAN(n2) || n1 <= 0. || n2 <= 0.)
	    ML_ERR_return_NAN;

	v1 = R_FINITE(n1) ? (rchisq(n1, rng) / n1) : 1;
	v2 = R_FINITE(n2) ? (rchisq(n2, rng) / n2) : 1;
	return v1 / v2;
    }

    R_RFUNC_INTERFACE_2ARG(rexp);

    /*
     *  D functions 
     */
    R_DFUNC_INTERFACE_4ARG(df, bnc_dF);

    /* 
     *  P functions
     */
    R_PFUNC_INTERFACE_5ARG(pf, bnc_pF);

    /* 
     *  Q functions
     */
    R_QFUNC_INTERFACE_5ARG(qf, bnc_qF);
    

}  // namespace bnc

#endif /* F_H */
