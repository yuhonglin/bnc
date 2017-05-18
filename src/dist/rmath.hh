#ifndef RMATH_H
#define RMATH_H

#include <common/math.hh>

namespace bnc {
    
/* Implmentation helpers */
#ifndef M_E
#define M_E		2.718281828459045235360287471353	/* e */
#endif

#ifndef M_LOG2E
#define M_LOG2E		1.442695040888963407359924681002	/* log2(e) */
#endif

#ifndef M_LOG10E
#define M_LOG10E	0.434294481903251827651128918917	/* log10(e) */
#endif

#ifndef M_LN2
#define M_LN2		0.693147180559945309417232121458	/* ln(2) */
#endif

#ifndef M_LN10
#define M_LN10		2.302585092994045684017991454684	/* ln(10) */
#endif

#ifndef M_PI
#define M_PI		3.141592653589793238462643383280	/* pi */
#endif

#ifndef M_2PI
#define M_2PI		6.283185307179586476925286766559	/* 2*pi */
#endif

#ifndef M_PI_2
#define M_PI_2		1.570796326794896619231321691640	/* pi/2 */
#endif

#ifndef M_PI_4
#define M_PI_4		0.785398163397448309615660845820	/* pi/4 */
#endif

#ifndef M_1_PI
#define M_1_PI		0.318309886183790671537767526745	/* 1/pi */
#endif

#ifndef M_2_PI
#define M_2_PI		0.636619772367581343075535053490	/* 2/pi */
#endif

#ifndef M_2_SQRTPI
#define M_2_SQRTPI	1.128379167095512573896158903122	/* 2/sqrt(pi) */
#endif

#ifndef M_SQRT2
#define M_SQRT2		1.414213562373095048801688724210	/* sqrt(2) */
#endif

#ifndef M_SQRT1_2
#define M_SQRT1_2	0.707106781186547524400844362105	/* 1/sqrt(2) */
#endif

/* R-Specific Constants */

#ifndef M_SQRT_3
#define M_SQRT_3	1.732050807568877293527446341506	/* sqrt(3) */
#endif

#ifndef M_SQRT_32
#define M_SQRT_32	5.656854249492380195206754896838	/* sqrt(32) */
#endif

#ifndef M_LOG10_2
#define M_LOG10_2	0.301029995663981195213738894724	/* log10(2) */
#endif

#ifndef M_SQRT_PI
#define M_SQRT_PI	1.772453850905516027298167483341	/* sqrt(pi) */
#endif

#ifndef M_1_SQRT_2PI
#define M_1_SQRT_2PI	0.398942280401432677939946059934	/* 1/sqrt(2pi) */
#endif

#ifndef M_SQRT_2dPI
#define M_SQRT_2dPI	0.797884560802865355879892119869	/* sqrt(2/pi) */
#endif


#ifndef M_LN_2PI
#define M_LN_2PI	1.837877066409345483560659472811	/* log(2*pi) */
#endif

#ifndef M_LN_SQRT_PI
#define M_LN_SQRT_PI	0.572364942924700087071713675677	/* log(sqrt(pi))
								   == log(pi)/2 */
#endif

#ifndef M_LN_SQRT_2PI
#define M_LN_SQRT_2PI	0.918938533204672741780329736406	/* log(sqrt(2*pi))
								   == log(2*pi)/2 */
#endif

#ifndef M_LN_SQRT_PId2
#define M_LN_SQRT_PId2	0.225791352644727432363097614947	/* log(sqrt(pi/2))
								   == log(pi/2)/2 */
#endif


#define SET_0(x) x = ((scale_p==NOR_P)?				\
		      0. :					\
		      -std::numeric_limits<double>::infinity())

#define SET_1(x) x = ((scale_p==NOR_P)? 1. : 0.)


#define RET_0 return ((scale_p==NOR_P)?				\
		      0. :					\
		      -std::numeric_limits<double>::infinity())

#define RET_1 return ((scale_p==NOR_P)? 1. : 0.)

#define RET_POSINF return std::numeric_limits<double>::infinity()

#define RET_NEGINF return -std::numeric_limits<double>::infinity()

#define RET_NAN(s) { LOG_WARNING(s);				\
	return std::numeric_limits<double>::quiet_NaN(); }

#define SET_NAN_RT(x, s) { LOG_WARNING(s);		\
	x = std::numeric_limits<double>::quiet_NaN();	\
	return; }						


#define R_D_Lval(p)	((tail_p==TAIL_LOWER) ? (p) : (0.5 - (p) + 0.5))	/*  p  */
#define R_D_Cval(p)	((tail_p==TAIL_LOWER) ? (0.5 - (p) + 0.5) : (p))	/*  1 - p */

#define R_DT_CIv(p)	((scale_p==LOG_P) ? ((tail_p==TAIL_LOWER) ?	\
					     -expm1(p) : exp(p))	\
			 : R_D_Cval(p))


#define R_DT_qIv(p)	((scale_p==LOG_P) ? ((tail_p==TAIL_LOWER)	\
					     ? exp(p) : - expm1(p))	\
			 : R_D_Lval(p))					\
			 

/* Do the boundaries exactly for q*() functions :
 * Often  _LEFT_ = ML_NEGINF , and very often _RIGHT_ = ML_POSINF;
 *
 * R_Q_P01_boundaries(p, _LEFT_, _RIGHT_)  :<==>
 *
 *     R_Q_P01_check(p);
 *     if (p == R_DT_0) return _LEFT_ ;
 *     if (p == R_DT_1) return _RIGHT_;
 *
 * the following implementation should be more efficient (less tests):
 */
#define R_Q_P01_boundaries(p, _LEFT_, _RIGHT_)			\
    if (scale_p==LOG_P) {					\
	if(p > 0)						\
	    RET_NAN("log(p)>0, return NaN");			\
	if(p == 0) /* upper bound*/				\
	    return (tail_p==TAIL_LOWER) ? _RIGHT_ : _LEFT_;	\
	if(p == -numeric_limits<double>::infinity())		\
	    return (tail_p==TAIL_LOWER) ? _LEFT_ : _RIGHT_;	\
    }								\
    else { /* !log_p */						\
	if(p < 0 || p > 1)					\
	    RET_NAN("p<0 or p>1, return NaN");			\
	if(p == 0)						\
	    return (tail_p==TAIL_LOWER) ? _LEFT_ : _RIGHT_;	\
	if(p == 1)						\
	    return (tail_p==TAIL_LOWER) ? _RIGHT_ : _LEFT_;	\
    }

/* additions for density functions (C.Loader) */
#define R_D_fexp(f,x)     ((scale_p==LOG_P) ? -0.5*log(f)+(x) : exp(x)/sqrt(f))

    
#define ML_ERR_return_NAN RET_NAN("return NaN")
#define ISNAN(x)   isnan(x)
#define ML_POSINF  (numeric_limits<double>::infinity())
#define ML_NEGINF  (-numeric_limits<double>::infinity())
#define ML_NAN  (numeric_limits<double>::quiet_NaN())
#define R_D__0     ((scale_p==NOR_P)?				\
		    0. :					\
		    -std::numeric_limits<double>::infinity())
#define R_D__1     ((scale_p==NOR_P)? 1. : 0.)
#ifndef DBL_MIN
#define DBL_MIN    (numeric_limits<double>::min())
#endif
#ifndef DBL_EPSILON
#define DBL_EPSILON numeric_limits<double>::epsilon()
#endif
#define R_Log1_Exp(x)   ((x) > -M_LN2 ? log(-expm1(x)) : log1p(-exp(x)))    

#define R_D_exp(x) ((scale_p==LOG_P)  ?  (x)	 : exp(x))	/* exp(x) */

#define ISNAN(x)   isnan(x)

#define R_DT_0	(tail_p==TAIL_LOWER ? R_D__0 : R_D__1)		/* 0 */
#define R_DT_1	(tail_p==TAIL_LOWER ? R_D__1 : R_D__0)		/* 1 */

#define R_P_bounds_01(x, x_min, x_max)		\
    if(x <= x_min) return R_DT_0;		\
    if(x >= x_max) return R_DT_1

}

#endif /* RMATH_H */
