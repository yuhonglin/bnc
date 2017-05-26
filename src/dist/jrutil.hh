#ifndef JRUTIL_H
#define JRUTIL_H

#include <type_traits>
#include <limits>
#include <cmath>
#include <cfloat>

#include <util/logger.hh>
#include <common/math.hh>
#include <dist/jrconst.hh>

namespace bnc {
    
#define ML_POSINF	numeric_limits<double>::infinity()
#define ML_NEGINF	-numeric_limits<double>::infinity()
#define ML_NAN		numeric_limits<double>::quiet_NaN()
#define ML_ERR_return_NAN { LOG_WARNING("return NaN"); return ML_NAN; }
    

    template <class RNGType>
    inline double unif_rand(RNGType *rng)
    {
	return rng->uniform();
    }

    template <class RNGType>
    inline double exp_rand(RNGType *rng)
    {
	return rng->exponential();
    }

    template <class RNGType>
    inline double norm_rand(RNGType *rng)
    {
	return rng->normal();
    }

    inline bool R_FINITE(const double &x) {
	return isfinite(x);
    }

    inline bool ISNAN(const double &x) {
	return isnan(x);
    }

    template<class T>
    inline T R_forceint(const T &x) {
	return nearbyint(x);
    }

    static double myfmod(double x1, double x2)
    {
	double q = x1 / x2;
	return x1 - floor(q) * x2;
    }

    double JR_pow(double x, double y) /* = x ^ y */
    {
	if(x == 1. || y == 0.)
	    return(1.);
	if(x == 0.) {
	    if(y > 0.) return(0.);
	    /* y < 0 */return(ML_POSINF);
	}
	if (R_FINITE(x) && R_FINITE(y))
	    return(pow(x,y));
	if (ISNAN(x) || ISNAN(y)) {
	    return(x + y);
	}
	if(!R_FINITE(x)) {
	    if(x > 0)		/* Inf ^ y */
		return((y < 0.)? 0. : ML_POSINF);
	    else {			/* (-Inf) ^ y */
		if(R_FINITE(y) && y == floor(y)) /* (-Inf) ^ n */
		    return((y < 0.) ? 0. : (myfmod(y,2.) ? x  : -x));
	    }
	}
	if(!R_FINITE(y)) {
	    if(x >= 0) {
		if(y > 0)		/* y == +Inf */
		    return((x >= 1)? ML_POSINF : 0.);
		else		/* y == -Inf */
		    return((x < 1) ? ML_POSINF : 0.);
	    }
	}
	return(ML_NAN);		/* all other cases: (-Inf)^{+-Inf,
				   non-int}; (neg)^{+-Inf} */
    }

    double JR_pow_di(double x, int n)
    {
	double pow = 1.0;

	if (ISNAN(x)) return x;
	if (n != 0) {
	    if (!R_FINITE(x)) return JR_pow(x, (double)n);
	    if (n < 0) { n = -n; x = 1/x; }
	    for(;;) {
		if(n & 01) pow *= x;
		if(n >>= 1) x *= x; else break;
	    }
	}
	return pow;
    }

#define R_RFUNC_INTERFACE_2ARG(FUNC)					\
    template<class RNGType>						\
    Vector FUNC(const int&n, const double& a, const double& b,		\
		 RNGType *rng)						\
    {									\
        Vector ret(n);							\
									\
	for (int i=0; i<n; i++) {					\
	    ret(i) = FUNC(a, b, rng);					\
	}								\
	return ret;							\
									\
    }									\
									\
    template<class RNGType, class AType, class BType>			\
    Vector FUNC(const int&n, const AType& a, const BType& b,		\
		RNGType *rng)						\
    {									\
        Vector ret(n);							\
	int ai=0, bi=0;							\
	for (int i=0; i<n; i++) {					\
	    ret(i) = FUNC(*(a.data()+ai), *(b.data()+bi),		\
		       rng);						\
	    if (++ai >= a.size()) ai = 0;				\
	    if (++bi >= b.size()) bi = 0;				\
	}								\
	return ret;							\
    }									\
									\
    template<class RNGType, class BType>				\
    Vector FUNC(const int&n, const double &a, const BType &b,		\
		 RNGType *rng)						\
    {									\
	Vector ret(n);							\
	int bi=0;							\
	for (int i=0; i<n; i++) {					\
	    ret(i) = FUNC(a, *(b.data()+bi), rng);			\
	    if (++bi >= b.size()) bi = 0;				\
	}								\
	return ret;							\
    }									\
									\
    template<class RNGType, class AType>				\
    Vector FUNC(const int&n, const AType &a, const double &b,		\
		 RNGType *rng)						\
    {									\
	Vector ret(n);							\
	int ai=0;							\
	for (int i=0; i<n; i++) {					\
	    ret(i) = FUNC(*(a.data()+ai), b, rng);			\
	    if (++ai >= a.size()) ai = 0;				\
	}								\
	return ret;							\
    }									\

}  // namespace bnc

#endif /* JRUTIL_H */
