// #ifndef PARETO_H
// #define PARETO_H

// #include <cmath>
// #include <limits>

// #include <util/logger.hh>
// #include <dist/brutil.hh>

// namespace bnc {
//     template<class RNGType>
//     double rpareto(const double& m, const double& s, RNGType *rng)
//     {
// 	ASSERT_MSG(std::isnan(m) || !std::isfinite(s) || m<=0 || s<=1,
// 	    "Invalid inputs.");
// 	return s*exp(rng->exponential()/m);
//     }

//     // Imitate R's rnorm function
//     R_RFUNC_INTERFACE_2ARG(rpareto);
    
//     /*
//      *  D functions 
//      */
//     double br_dpareto(const double&x, const double& m, const double& s,
// 		   const int& give_log)
//     {
// 	if (give_log == LOG) {
// 	    return log(m) + m*log(s) - (m+1)*log(s);
// 	} else {
// 	    return m/pow(x,m+1)*pow(s,m);
// 	}
//     }
//     R_DFUNC_INTERFACE_4ARG(dpareto, br_dpareto);

//     /*
//      *  Q functions
//      */
//     double br_qpareto(const double& p, const double& m, const double& s,
// 		      const int& lower_tail, const int& log_p)
//     {
// 	ASSERT_MSG( (log_p==LOG && p>0) || (log_p==NORMAL && (p<0 || p>1)),
// 		    "Invalid inputs." );
// 	double logp;
// 	if (lower_tail==LOWTAIL)
// 	{
// 	    if (log_p==LOG)
// 		logp = p;
// 	    else
// 		logp = log(p);
// 	} else
// 	{
// 	    if (log_p==LOG)
// 		logp = log(1-exp(p));
// 	    else
// 		logp = log(1-p);
// 	}

// 	return exp(log(s)-logp/m);
//     }
//     R_QFUNC_INTERFACE_5ARG(qpareto, br_ppareto);

//     /* 
//      *  P functions
//      */
//     double br_ppareto(const double& x, const double& m, const double& s,
// 		      const int& lower_tail, const int& give_log)
//     {
// 	if (x < s)
// 	{
// 	    if (give_log==LOG)
// 	    {
// 		return -std::numeric_limits<double>::quiet_NaN();
// 	    }
// 	    else
// 	    {
// 		return 0.;
// 	    }
// 	}
	
// 	double logq = m*log(s/x);
// 	if (lower_tail==UPTAIL)
// 	    return give_log==LOG ? logq : exp(logq);
// 	else
// 	    return give_log==LOG ? log(1-exp(logq)) : 1-exp(logq);
//     }
//     R_PFUNC_INTERFACE_5ARG(ppareto, br_ppareto);
    
// }  // namespace bnc

// #endif /* PARETO_H */
