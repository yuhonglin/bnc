#ifndef GAMMA_H
#define GAMMA_H

#include <limits>
#include <cmath>
#include <iostream>

using namespace std;

#include <rng/rng.hh>
#include <util/logger.hh>
#include <dist/dist.hh>
#include <dist/rmath.hh>
#include <dist/lgamma.hh>
#include <dist/chebyshev.hh>
#include <dist/dpois.hh>

#define ISNAN(x)   isnan(x)
#define give_log   (scale_p==LOG_P)
#define log_p      (scale_p==LOG_P)
#define lower_tail (tail_p==TAIL_LOWER)

namespace bnc {

    template<SCALE_P scale_p>
    double dpois_raw(double x, double lambda);

    template<TAIL_P tail_p, SCALE_P scale_p>
    double pgamma_smallx (const double& x, const double& alph);
    
    /* Density function */
    template<SCALE_P scale_p=NOR_P>
    double dgamma(double x, double shape, double scale)
    {
	double pr;
#ifdef IEEE_754
	if (ISNAN(x) || ISNAN(shape) || ISNAN(scale))
	    return x + shape + scale;
#endif
	if (shape < 0 || scale <= 0) ML_ERR_return_NAN;
	if (x < 0)
	    return R_D__0;
	if (shape == 0) /* point mass at 0 */
	    return (x == 0)? ML_POSINF : R_D__0;
	if (x == 0) {
	    if (shape < 1) return ML_POSINF;
	    if (shape > 1) return R_D__0;
	    /* else */
	    return give_log ? -log(scale) : 1 / scale;
	}

	if (shape < 1) {
	    pr = dpois_raw<scale_p>(shape, x/scale);
	    return give_log ?  pr + log(shape/x) : pr*shape/x;
	}
	/* else  shape >= 1 */
	pr = dpois_raw<scale_p>(shape-1, x/scale);
	return give_log ? pr - log(scale) : pr/scale;
    }

    /* probability function */
#define SQR(x) ((x)*(x))
    static const double scalefactor = SQR(SQR(SQR(4294967296.0)));
#undef SQR

/* If |x| > |k| * M_cutoff,  then  log[ exp(-x) * k^x ]	 =~=  -x */
static const double M_cutoff = M_LN2 * DBL_MAX_EXP / DBL_EPSILON;/*=3.196577e18*/

/* Continued fraction for calculation of
 *    1/i + x/(i+d) + x^2/(i+2*d) + x^3/(i+3*d) + ... = sum_{k=0}^Inf x^k/(i+k*d)
 *
 * auxilary in log1pmx() and lgamma1p()
 */
    double
    logcf (const double& x, const double& i, const double& d,
	   const double& eps /* ~ relative tolerance */)
    {
	double c1 = 2 * d;
	double c2 = i + d;
	double c4 = c2 + d;
	double a1 = c2;
	double b1 = i * (c2 - i * x);
	double b2 = d * d * x;
	double a2 = c4 * c2 - b2;

	b2 = c4 * b1 - i * b2;

	while (fabs(a2 * b1 - a1 * b2) > fabs(eps * b1 * b2)) {
	    double c3 = c2*c2*x;
	    c2 += d;
	    c4 += d;
	    a1 = c4 * a2 - c3 * a1;
	    b1 = c4 * b2 - c3 * b1;

	    c3 = c1 * c1 * x;
	    c1 += d;
	    c4 += d;
	    a2 = c4 * a1 - c3 * a2;
	    b2 = c4 * b1 - c3 * b2;

	    if (fabs (b2) > scalefactor) {
		a1 /= scalefactor;
		b1 /= scalefactor;
		a2 /= scalefactor;
		b2 /= scalefactor;
	    } else if (fabs (b2) < 1 / scalefactor) {
		a1 *= scalefactor;
		b1 *= scalefactor;
		a2 *= scalefactor;
		b2 *= scalefactor;
	    }
	}

	return a2 / b2;
    }

/* Accurate calculation of log(1+x)-x, particularly for small x.  */
    double log1pmx (double x)
    {
	static const double minLog1Value = -0.79149064;

	if (x > 1 || x < minLog1Value)
	    return log1p(x) - x;
	else { /* -.791 <=  x <= 1  -- expand in  [x/(2+x)]^2 =: y :
		* log(1+x) - x =  x/(2+x) * [ 2 * y * S(y) - x],  with
		* ---------------------------------------------
		* S(y) = 1/3 + y/5 + y^2/7 + ... = \sum_{k=0}^\infty  y^k / (2k + 3)
		*/
	    double r = x / (2 + x), y = r * r;
	    if (fabs(x) < 1e-2) {
		static const double two = 2;
		return r * ((((two / 9 * y + two / 7) * y + two / 5) * y +
			     two / 3) * y - x);
	    } else {
		static const double tol_logcf = 1e-14;
		return r * (2 * y * logcf (y, 3, 2, tol_logcf) - x);
	    }
	}
    }

/* Compute  log(gamma(a+1))  accurately also for small a (0 < a < 0.5). */
    double lgamma1p (const double& a)
    {
	const double eulers_const =	 0.5772156649015328606065120900824024;

	/* coeffs[i] holds (zeta(i+2)-1)/(i+2) , i = 0:(N-1), N = 40 : */
	#undef N
	const int N = 40;
	static const double coeffs[40] = {
	    0.3224670334241132182362075833230126e-0,/* = (zeta(2)-1)/2 */
	    0.6735230105319809513324605383715000e-1,/* = (zeta(3)-1)/3 */
	    0.2058080842778454787900092413529198e-1,
	    0.7385551028673985266273097291406834e-2,
	    0.2890510330741523285752988298486755e-2,
	    0.1192753911703260977113935692828109e-2,
	    0.5096695247430424223356548135815582e-3,
	    0.2231547584535793797614188036013401e-3,
	    0.9945751278180853371459589003190170e-4,
	    0.4492623673813314170020750240635786e-4,
	    0.2050721277567069155316650397830591e-4,
	    0.9439488275268395903987425104415055e-5,
	    0.4374866789907487804181793223952411e-5,
	    0.2039215753801366236781900709670839e-5,
	    0.9551412130407419832857179772951265e-6,
	    0.4492469198764566043294290331193655e-6,
	    0.2120718480555466586923135901077628e-6,
	    0.1004322482396809960872083050053344e-6,
	    0.4769810169363980565760193417246730e-7,
	    0.2271109460894316491031998116062124e-7,
	    0.1083865921489695409107491757968159e-7,
	    0.5183475041970046655121248647057669e-8,
	    0.2483674543802478317185008663991718e-8,
	    0.1192140140586091207442548202774640e-8,
	    0.5731367241678862013330194857961011e-9,
	    0.2759522885124233145178149692816341e-9,
	    0.1330476437424448948149715720858008e-9,
	    0.6422964563838100022082448087644648e-10,
	    0.3104424774732227276239215783404066e-10,
	    0.1502138408075414217093301048780668e-10,
	    0.7275974480239079662504549924814047e-11,
	    0.3527742476575915083615072228655483e-11,
	    0.1711991790559617908601084114443031e-11,
	    0.8315385841420284819798357793954418e-12,
	    0.4042200525289440065536008957032895e-12,
	    0.1966475631096616490411045679010286e-12,
	    0.9573630387838555763782200936508615e-13,
	    0.4664076026428374224576492565974577e-13,
	    0.2273736960065972320633279596737272e-13,
	    0.1109139947083452201658320007192334e-13/* = (zeta(40+1)-1)/(40+1) */
	};

	const double c = 0.2273736845824652515226821577978691e-12;/* zeta(N+2)-1 */
	const double tol_logcf = 1e-14;
	double lgam;
	int i;

	if (fabs (a) >= 0.5)
	    return lgammafn (a + 1);

	/* Abramowitz & Stegun 6.1.33 : for |x| < 2,
	 * <==> log(gamma(1+x)) = -(log(1+x) - x) - gamma*x + x^2 * \sum_{n=0}^\infty c_n (-x)^n
	 * where c_n := (Zeta(n+2) - 1)/(n+2)  = coeffs[n]
	 *
	 * Here, another convergence acceleration trick is used to compute
	 * lgam(x) :=  sum_{n=0..Inf} c_n (-x)^n
	 */
	lgam = c * logcf(-a / 2, N + 2, 1, tol_logcf);
	for (i = N - 1; i >= 0; i--)
	    lgam = coeffs[i] - a * lgam;

	return (a * lgam - eulers_const) * a - log1pmx (a);
    } /* lgamma1p */

/*
 * Compute the log of a sum from logs of terms, i.e.,
 *
 *     log (exp (logx) + exp (logy))
 *
 * without causing overflows and without throwing away large handfuls
 * of accuracy.
 */
    double logspace_add (const double& logx, const double& logy)
    {
	return fmax2 (logx, logy) + log1p (exp (-fabs (logx - logy)));
    }

/*
 * Compute the log of a difference from logs of terms, i.e.,
 *
 *     log (exp (logx) - exp (logy))
 *
 * without causing overflows and without throwing away large handfuls
 * of accuracy.
 */
    double logspace_sub (const double& logx, const double& logy)
    {
	return logx + R_Log1_Exp(logy - logx);
    }

/* dpois_wrap (x_P_1,  lambda, g_log) ==
 *   dpois (x_P_1 - 1, lambda, g_log) :=  exp(-L)  L^k / gamma(k+1) ,  k := x_P_1 - 1
 */
    template<SCALE_P scale_p=NOR_P>
    double dpois_wrap (const double& x_plus_1, const double& lambda)
    {
	if (isinf(lambda))
	    return R_D__0;
	if (x_plus_1 > 1)
	    return dpois_raw<scale_p> (x_plus_1 - 1, lambda);
	if (lambda > fabs(x_plus_1 - 1) * M_cutoff)
	    return R_D_exp(-lambda - lgammafn(x_plus_1));
	else {
	    double d = dpois_raw<scale_p> (x_plus_1, lambda);
	    return give_log
		? d + log (x_plus_1 / lambda)
		: d * (x_plus_1 / lambda);
	}
    }

/*
 * Abramowitz and Stegun 6.5.29 [right]
 */
    template<TAIL_P tail_p=TAIL_LOWER, SCALE_P scale_p=NOR_P>
    double pgamma_smallx (const double& x, const double& alph)
    {
	double sum = 0, c = alph, n = 0, term;

	/*
	 * Relative to 6.5.29 all terms have been multiplied by alph
	 * and the first, thus being 1, is omitted.
	 */

	do {
	    n++;
	    c *= -x / n;
	    term = c / (alph + n);
	    sum += term;
	} while (fabs (term) > DBL_EPSILON * fabs (sum));

	if (lower_tail) {
	    double f1 = log_p ? log1p (sum) : 1 + sum;
	    double f2;
	    if (alph > 1) {
		f2 = dpois_raw<scale_p> (alph, x);
		f2 = log_p ? f2 + x : f2 * exp (x);
	    } else if (log_p)
		f2 = alph * log (x) - lgamma1p (alph);
	    else
		f2 = pow (x, alph) / exp (lgamma1p (alph));
	    return log_p ? f1 + f2 : f1 * f2;
	} else {
	    double lf2 = alph * log (x) - lgamma1p (alph);
	    if (log_p)
		return R_Log1_Exp (log1p (sum) + lf2);
	    else {
		double f1m1 = sum;
		double f2m1 = expm1 (lf2);
		return -(f1m1 + f2m1 + f1m1 * f2m1);
	    }
	}
    } /* pgamma_smallx() */
    
    template<SCALE_P scale_p=NOR_P>
    double pd_upper_series (const double& x, double y)
    {
	double term = x / y;
	double sum = term;

	do {
	    y++;
	    term *= x / y;
	    sum += term;
	} while (term > sum * DBL_EPSILON);

	/* sum =  \sum_{n=1}^ oo  x^n     / (y*(y+1)*...*(y+n-1))
	 *	   =  \sum_{n=0}^ oo  x^(n+1) / (y*(y+1)*...*(y+n))
	 *	   =  x/y * (1 + \sum_{n=1}^oo	x^n / ((y+1)*...*(y+n)))
	 *	   ~  x/y +  o(x/y)   {which happens when alph -> Inf}
	 */
	return log_p ? log (sum) : sum;
    }
/* Continued fraction for calculation of
 *    scaled upper-tail F_{gamma}
 *  ~=  (y / d) * [1 +  (1-y)/d +  O( ((1-y)/d)^2 ) ]
 */
    double pd_lower_cf (double y, double d)
    {
	double f= 0.0 /* -Wall */, of, f0;
	double i, c2, c3, c4,  a1, b1,  a2, b2;

#define	NEEDED_SCALE				\
	(b2 > scalefactor) {			\
	    a1 /= scalefactor;			\
	    b1 /= scalefactor;			\
	    a2 /= scalefactor;			\
	    b2 /= scalefactor;			\
	}

#define max_it 200000

	if (y == 0) return 0;

	f0 = y/d;
	/* Needed, e.g. for  pgamma(10^c(100,295), shape= 1.1, log=TRUE): */
	if(fabs(y - 1) < fabs(d) * DBL_EPSILON) { /* includes y < d = Inf */
	    return (f0);
	}

	if(f0 > 1.) f0 = 1.;
	c2 = y;
	c4 = d; /* original (y,d), *not* potentially scaled ones!*/

	a1 = 0; b1 = 1;
	a2 = y; b2 = d;

        while NEEDED_SCALE

	    i = 0; of = -1.; /* far away */
	while (i < max_it) {

	    i++;	c2--;	c3 = i * c2;	c4 += 2;
	    /* c2 = y - i,  c3 = i(y - i),  c4 = d + 2i,  for i odd */
	    a1 = c4 * a2 + c3 * a1;
	    b1 = c4 * b2 + c3 * b1;

	    i++;	c2--;	c3 = i * c2;	c4 += 2;
	    /* c2 = y - i,  c3 = i(y - i),  c4 = d + 2i,  for i even */
	    a2 = c4 * a1 + c3 * a2;
	    b2 = c4 * b1 + c3 * b2;

	    if NEEDED_SCALE

		if (b2 != 0) {
		    f = a2 / b2;
		    /* convergence check: relative; "absolute" for very small f : */
		    if (fabs (f - of) <= DBL_EPSILON * fmax2(f0, fabs(f))) {
			return f;
		    }
		    of = f;
		}
	}
	LOG_WARNING("NON-convergence");
	return f;/* should not happen ... */
    } /* pd_lower_cf() */
#undef NEEDED_SCALE


    static double
    pd_lower_series (const double& lambda, double y)
    {
	double term = 1, sum = 0;

	while (y >= 1 && term > sum * DBL_EPSILON) {
	    term *= y / lambda;
	    sum += term;
	    y--;
	}
	/* sum =  \sum_{n=0}^ oo  y*(y-1)*...*(y - n) / lambda^(n+1)
	 *	   =  y/lambda * (1 + \sum_{n=1}^Inf  (y-1)*...*(y-n) / lambda^n)
	 *	   ~  y/lambda + o(y/lambda)
	 */
	if (y != floor (y)) {
	    /*
	     * The series does not converge as the terms start getting
	     * bigger (besides flipping sign) for y < -lambda.
	     */
	    double f;
	    /* FIXME: in quite few cases, adding  term*f  has no effect (f too small)
	     *	  and is unnecessary e.g. for pgamma(4e12, 121.1) */
	    f = pd_lower_cf (y, lambda + 1 - y);
	    sum += term * f;
	}

	return sum;
    } /* pd_lower_series() */

/*
 * Compute the following ratio with higher accuracy that would be had
 * from doing it directly.
 *
 *		 dnorm (x, 0, 1, FALSE)
 *	   ----------------------------------
 *	   pnorm (x, 0, 1, lower_tail, FALSE)
 *
 * Abramowitz & Stegun 26.2.12
 */
    template<TAIL_P tail_p=TAIL_LOWER>
    static double
    dpnorm (double x, double lp)
    {
	/*
	 * So as not to repeat a pnorm call, we expect
	 *
	 *	 lp == pnorm (x, 0, 1, lower_tail, TRUE)
	 *
	 * but use it only in the non-critical case where either x is small
	 * or p==exp(lp) is close to 1.
	 */
	bool lowertail = (tail_p==TAIL_LOWER);
	if (x < 0) {
	    x = -x;
	    lowertail = !lowertail;
	}

	if (x > 10 && !lowertail) {
	    double term = 1 / x;
	    double sum = term;
	    double x2 = x * x;
	    double i = 1;

	    do {
		term *= -i / x2;
		sum += term;
		i += 2;
	    } while (fabs (term) > DBL_EPSILON * sum);

	    return 1 / sum;
	} else {
	    double d = dnorm <LOG_P> (x, 0., 1.);
	    return d / exp (lp);
	}
    }

/*
 * Asymptotic expansion to calculate the probability that Poisson variate
 * has value <= x.
 * Various assertions about this are made (without proof) at
 * http://members.aol.com/iandjmsmith/PoissonApprox.htm
 */
    template<TAIL_P tail_p=TAIL_LOWER, SCALE_P scale_p=NOR_P>
    static double
    ppois_asymp (double x, double lambda)
    {
	static const double coefs_a[8] = {
	    -1e99, /* placeholder used for 1-indexing */
	    2/3.,
	    -4/135.,
	    8/2835.,
	    16/8505.,
	    -8992/12629925.,
	    -334144/492567075.,
	    698752/1477701225.
	};

	static const double coefs_b[8] = {
	    -1e99, /* placeholder */
	    1/12.,
	    1/288.,
	    -139/51840.,
	    -571/2488320.,
	    163879/209018880.,
	    5246819/75246796800.,
	    -534703531/902961561600.
	};

	double elfb, elfb_term;
	double res12, res1_term, res1_ig, res2_term, res2_ig;
	double dfm, pt_, s2pt, f, np;
	int i;

	dfm = lambda - x;
	/* If lambda is large, the distribution is highly concentrated
	   about lambda.  So representation error in x or lambda can lead
	   to arbitrarily large values of pt_ and hence divergence of the
	   coefficients of this approximation.
	*/
	pt_ = - log1pmx (dfm / x);
	s2pt = sqrt (2 * x * pt_);
	if (dfm < 0) s2pt = -s2pt;

	res12 = 0;
	res1_ig = res1_term = sqrt (x);
	res2_ig = res2_term = s2pt;
	for (i = 1; i < 8; i++) {
	    res12 += res1_ig * coefs_a[i];
	    res12 += res2_ig * coefs_b[i];
	    res1_term *= pt_ / i ;
	    res2_term *= 2 * pt_ / (2 * i + 1);
	    res1_ig = res1_ig / x + res1_term;
	    res2_ig = res2_ig / x + res2_term;
	}

	elfb = x;
	elfb_term = 1;
	for (i = 1; i < 8; i++) {
	    elfb += elfb_term * coefs_b[i];
	    elfb_term /= x;
	}
	if (!lower_tail) elfb = -elfb;

	f = res12 / elfb;

	np = pnorm<TAIL_UPPER, scale_p> (s2pt, 0.0, 1.0);

	if (log_p) {
	    double n_d_over_p;
	    if (tail_p==TAIL_LOWER)
		n_d_over_p = dpnorm<TAIL_UPPER> (s2pt, np);
	    else
		n_d_over_p = dpnorm<TAIL_LOWER> (s2pt, np);
	    return np + log1p (f * n_d_over_p);
	} else {
	    double nd = dnorm<scale_p> (s2pt, 0., 1.);
	    return np + f * nd;
	}
    } /* ppois_asymp() */

    template<TAIL_P tail_p=TAIL_LOWER, SCALE_P scale_p=NOR_P>
    double pgamma_raw (double x, double alph)
    {
/* Here, assume that  (x,alph) are not NA  &  alph > 0 . */

	double res;

	R_P_bounds_01(x, 0., ML_POSINF);

	if (x < 1) {
	    res = pgamma_smallx<tail_p, scale_p> (x, alph);
	} else if (x <= alph - 1 && x < 0.8 * (alph + 50)) {
	    /* incl. large alph compared to x */
	    double sum = pd_upper_series<scale_p> (x, alph);/* = x/alph + o(x/alph) */
	    double d = dpois_wrap<scale_p> (alph, x);

	    if (!lower_tail)
		res = log_p
		    ? R_Log1_Exp (d + sum)
		    : 1 - d * sum;
	    else
		res = log_p ? sum + d : sum * d;
	} else if (alph - 1 < x && alph < 0.8 * (x + 50)) {
	    /* incl. large x compared to alph */
	    double sum;
	    double d = dpois_wrap<scale_p> (alph, x);

	    if (alph < 1) {
		if (x * DBL_EPSILON > 1 - alph)
		    sum = R_D__1;
		else {
		    double f = pd_lower_cf (alph, x - (alph - 1)) * x / alph;
		    /* = [alph/(x - alph+1) + o(alph/(x-alph+1))] * x/alph = 1 + o(1) */
		    sum = log_p ? log (f) : f;
		}
	    } else {
		sum = pd_lower_series (x, alph - 1);/* = (alph-1)/x + o((alph-1)/x) */
		sum = log_p ? log1p (sum) : 1 + sum;
	    }

	    if (!lower_tail)
		res = log_p ? sum + d : sum * d;
	    else
		res = log_p
		    ? R_Log1_Exp (d + sum)
		    : 1 - d * sum;
	} else { /* x >= 1 and x fairly near alph. */

	    if (tail_p==TAIL_LOWER)
		res = ppois_asymp<TAIL_UPPER, scale_p> (alph - 1, x);
	    else
		res = ppois_asymp<TAIL_LOWER, scale_p> (alph - 1, x);
	}

	/*
	 * We lose a fair amount of accuracy to underflow in the cases
	 * where the final result is very close to DBL_MIN.	 In those
	 * cases, simply redo via log space.
	 */
	if (!log_p && res < DBL_MIN / DBL_EPSILON) {
	    /* with(.Machine, double.xmin / double.eps) #|-> 1.002084e-292 */

	    return exp (pgamma_raw<tail_p, scale_p> (x, alph));
	} else
	    return res;
    }

    template<TAIL_P tail_p=TAIL_LOWER, SCALE_P scale_p=NOR_P>
    double pgamma(double x, const double& alph, const double& scale)
    {
	if (isnan(x) || isnan(alph) || isnan(scale))
	    return x + alph + scale;

	if(alph < 0. || scale <= 0.)
	    ML_ERR_return_NAN;
	x /= scale;

	if (isnan(x)) /* eg. original x = scale = +Inf */
	    return x;

	if(alph == 0.) /* limit case; useful e.g. in pnchisq() */
	    return (x <= 0) ? R_DT_0: R_DT_1; /* <= assert  pgamma(0,0) ==> 0 */

	return pgamma_raw<tail_p, scale_p> (x, alph);
    }
/* Comment To ABOVE
 * From: terra@gnome.org (Morten Welinder)
 * To: R-bugs@biostat.ku.dk
 * Cc: maechler@stat.math.ethz.ch
 * Subject: Re: [Rd] pgamma discontinuity (PR#7307)
 * Date: Tue, 11 Jan 2005 13:57:26 -0500 (EST)

 * this version of pgamma appears to be quite good and certainly a vast
 * improvement over current R code.  (I last looked at 2.0.1)  Apart from
 * type naming, this is what I have been using for Gnumeric 1.4.1.

 * This could be included into R as-is, but you might want to benefit from
 * making logcf, log1pmx, lgamma1p, and possibly logspace_add/logspace_sub
 * available to other parts of R.

 * MM: I've not (yet?) taken  logcf(), but the other four
 */


/* Quantile function */
    template<TAIL_P tail_p=TAIL_LOWER, SCALE_P scale_p=NOR_P>
    double qchisq_appr(const double &p, const double &nu,
		       const double &g /* = log Gamma(nu/2) */,
		       const double &tol /* EPS1 */)
    {
#define C7	4.67
#define C8	6.66
#define C9	6.73
#define C10	13.32

	double alpha, a, c, ch, p1;
	double p2, q, t, x;

	/* test arguments and initialise */

	if (isnan(p) || isnan(nu))
	    return p + nu;

	R_Q_P01_check(p);
	if (nu <= 0) ML_ERR_return_NAN;

	alpha = 0.5 * nu;/* = [pq]gamma() shape */
	c = alpha-1;

	if(nu < (-1.24)*(p1 = R_DT_log(p))) {	/* for small chi-squared */
	    /* log(alpha) + g = log(alpha) + log(gamma(alpha)) =
	     *        = log(alpha*gamma(alpha)) = lgamma(alpha+1) suffers from
	     *  catastrophic cancellation when alpha << 1
	     */
	    double lgam1pa = (alpha < 0.5) ? lgamma1p(alpha) : (log(alpha) + g);
	    ch = exp((lgam1pa + p1)/alpha + M_LN2);

	} else if(nu > 0.32) {	/*  using Wilson and Hilferty estimate */

	    x = qnorm<tail_p, scale_p>(p, 0, 1);
	    p1 = 2./(9*nu);
	    ch = nu*pow(x*sqrt(p1) + 1-p1, 3);

	    /* approximation for p tending to 1: */
	    if( ch > 2.2*nu + 6 )
		ch = -2*(R_DT_Clog(p) - c*log(0.5*ch) + g);

	} else { /* "small nu" : 1.24*(-log(p)) <= nu <= 0.32 */

	    ch = 0.4;
	    a = R_DT_Clog(p) + g + c*M_LN2;

	    do {
		q = ch;
		p1 = 1. / (1+ch*(C7+ch));
		p2 = ch*(C9+ch*(C8+ch));
		t = -0.5 +(C7+2*ch)*p1 - (C9+ch*(C10+3*ch))/p2;
		ch -= (1- exp(a+0.5*ch)*p2*p1)/t;
	    } while(fabs(q - ch) > tol * fabs(ch));
	}

	return ch;
    }

    template<TAIL_P tail_p=TAIL_LOWER, SCALE_P scale_p=NOR_P>
    double qgamma(double p, const double &alpha, const double &scale)
/*			shape = alpha */
    {
#define EPS1 1e-2
#define EPS2 5e-7/* final precision of AS 91 */
#define EPS_N 1e-15/* precision of Newton step / iterations */
#define LN_EPS -36.043653389117156 /* = log(.Machine$double.eps) iff IEEE_754 */

#define MAXIT 1000/* was 20 */

#define pMIN 1e-100   /* was 0.000002 = 2e-6 */
#define pMAX (1-1e-14)/* was (1-1e-12) and 0.999998 = 1 - 2e-6 */

	const static double
	    i420  = 1./ 420.,
	    i2520 = 1./ 2520.,
	    i5040 = 1./ 5040;

	double p_, a, b, c, g, ch, ch0, p1;
	double p2, q, s1, s2, s3, s4, s5, s6, t, x;
	int i, max_it_Newton = 1;

	/* test arguments and initialise */

	if (ISNAN(p) || ISNAN(alpha) || ISNAN(scale))
	    return p + alpha + scale;

	R_Q_P01_boundaries(p, 0., ML_POSINF);

	if (alpha < 0 || scale <= 0) ML_ERR_return_NAN;

	if (alpha == 0) /* all mass at 0 : */ return 0.;

	if (alpha < 1e-10) {
	    /* Warning seems unnecessary now: */
	    LOG_WARNING("value of shape is extremely small: results may be unreliable");
	    max_it_Newton = 7;/* may still be increased below */
	}

	p_ = R_DT_qIv(p);/* lower_tail prob (in any case) */

	g = lgammafn(alpha);/* log Gamma(v/2) */

	/*----- Phase I : Starting Approximation */
	ch = qchisq_appr<tail_p, scale_p>(p, /* nu= 'df' =  */ 2*alpha, /* lgamma(nu/2)= */ g,
					  /* tol= */ EPS1);
	if(!isfinite(ch)) {
	    /* forget about all iterations! */
	    max_it_Newton = 0; goto END;
	}
	if(ch < EPS2) {/* Corrected according to AS 91; MM, May 25, 1999 */
	    max_it_Newton = 20;
	    goto END;/* and do Newton steps */
	}

	/* FIXME: This (cutoff to {0, +Inf}) is far from optimal
	 * -----  when log_p or !lower_tail, but NOT doing it can be even worse */
	if(p_ > pMAX || p_ < pMIN) {
	    /* did return ML_POSINF or 0.;	much better: */
	    max_it_Newton = 20;
	    goto END;/* and do Newton steps */
	}

/*----- Phase II: Iteration
 *	Call pgamma() [AS 239]	and calculate seven term taylor series
 */
	c = alpha-1;
	s6 = (120+c*(346+127*c)) * i5040; /* used below, is "const" */

	ch0 = ch;/* save initial approx. */
	for(i=1; i <= MAXIT; i++ ) {
	    q = ch;
	    p1 = 0.5*ch;
	    p2 = p_ - pgamma_raw<TAIL_LOWER, NOR_P>(p1, alpha);

	    if(!R_FINITE(p2) || ch <= 0)
	    { ch = ch0; max_it_Newton = 27; goto END; }/*was  return ML_NAN;*/

	    t = p2*exp(alpha*M_LN2+g+p1-c*log(ch));
	    b = t/ch;
	    a = 0.5*t - b*c;
	    s1 = (210+ a*(140+a*(105+a*(84+a*(70+60*a))))) * i420;
	    s2 = (420+ a*(735+a*(966+a*(1141+1278*a)))) * i2520;
	    s3 = (210+ a*(462+a*(707+932*a))) * i2520;
	    s4 = (252+ a*(672+1182*a) + c*(294+a*(889+1740*a))) * i5040;
	    s5 = (84+2264*a + c*(1175+606*a)) * i2520;

	    ch += t*(1+0.5*t*s1-b*c*(s1-b*(s2-b*(s3-b*(s4-b*(s5-b*s6))))));
	    if(fabs(q - ch) < EPS2*ch)
		goto END;
	    if(fabs(q - ch) > 0.1*ch) {/* diverging? -- also forces ch > 0 */
		if(ch < q) ch = 0.9 * q; else ch = 1.1 * q;
	    }
	}
/* no convergence in MAXIT iterations -- but we add Newton now... */
	LOG_WARNING("not converged");
/* was
 *    ML_ERROR(ME_PRECISION, "qgamma");
 * does nothing in R !*/

    END:
/* PR# 2214 :	 From: Morten Welinder <terra@diku.dk>, Fri, 25 Oct 2002 16:50
   --------	 To: R-bugs@biostat.ku.dk     Subject: qgamma precision

   * With a final Newton step, double accuracy, e.g. for (p= 7e-4; nu= 0.9)
   *
   * Improved (MM): - only if rel.Err > EPS_N (= 1e-15);
   *		    - also for lower_tail = FALSE	 or log_p = TRUE
   * 		    - optionally *iterate* Newton
   */
	x = 0.5*scale*ch;
    	        
	if(max_it_Newton) {
	    /* always use log scale */
	    if (scale_p!=LOG_P) {
		p = log(p);
	    }

	    if(x == 0) {
		const double _1_p = 1. + 1e-7;
		const double _1_m = 1. - 1e-7;
		x = DBL_MIN;
		p_ = pgamma<tail_p, LOG_P>(x, alpha, scale);
		if(( lower_tail && p_ > p * _1_p) ||
		   (!lower_tail && p_ < p * _1_m))
		    return(0.);
		/* else:  continue, using x = DBL_MIN instead of  0  */
	    }
	    else
		p_ = pgamma<tail_p, LOG_P>(x, alpha, scale);

	    if(p_ == ML_NEGINF) return 0; /* PR#14710 */
	    for(i = 1; i <= max_it_Newton; i++) {
		p1 = p_ - p;

		if(fabs(p1) < fabs(EPS_N * p))
		    break;
		/* else */
		if((g = dgamma<LOG_P>(x, alpha, scale)) == R_D__0) {
		    break;
		}
		/* else :
		 * delta x = f(x)/f'(x);
		 * if(log_p) f(x) := log P(x) - p; f'(x) = d/dx log P(x) = P' / P
		 * ==> f(x)/f'(x) = f*P / P' = f*exp(p_) / P' (since p_ = log P(x))
		 */
		t = log_p ? p1*exp(p_ - g) : p1/g ;/* = "delta x" */
		t = lower_tail ? x - t : x + t;
		p_ = pgamma<tail_p, LOG_P>(t, alpha, scale);
		if (fabs(p_ - p) > fabs(p1) ||
		    (i > 1 && fabs(p_ - p) == fabs(p1)) /* <- against flip-flop */) {
		    /* no improvement */
		    break;
		} /* else : */

		x = t;
	    }
	}

	return x;
    }    

} // namespace bnc

#endif /* GAMMA_H */
