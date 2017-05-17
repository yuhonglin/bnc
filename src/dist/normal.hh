#ifndef NORMAL_H
#define NORMAL_H

#include <limits>
#include <cmath>

using namespace std;

#include <dist/dist.hh>
#include <dist/rmath.hh>
#include <rng/rng.hh>
#include <util/logger.hh>


namespace bnc {

    // forward declaration
    template<TAIL_P tail_p=TAIL_LOWER, SCALE_P scale_p=NOR_P>
    void pnorm_both(double x, double *cum, double *ccum);
    
    /* Density function */
    template<SCALE_P scale_p=NOR_P>
    double dnorm(const double& xi, const double &mu, const double &sigma) {
	if (isnan(xi) || isnan(mu) || isnan(sigma))
	    RET_NAN("inputs contain NaN, return NaN.");
	    
	if (sigma < 0)
	    RET_NAN("sd<0, return NaN");
	
	if (isinf(xi) && mu==xi)
	    RET_NAN("x and mu are infs, return NaN");
	
	if (sigma == 0) { if (xi == mu) { RET_POSINF; } else { RET_0; } }

	double x = (xi - mu)/sigma;
	
	if (isinf(x)) RET_0;    // sd is very small

	x = abs(x);
	if (x >= 2*sqrt(numeric_limits<double>::max())) RET_0;
	
	if (scale_p==LOG_P)
	    return -(M_LN_SQRT_2PI + 0.5 * x * x + log(sigma));

#ifdef BNC_FAST_dnorm
	// and for R <= 3.0.x and R-devel upto 2014-01-01:
	return M_1_SQRT_2PI * exp(-0.5 * x * x) / sigma;
#else
	// more accurate, less fast :
	if (x < 5) return M_1_SQRT_2PI * exp(-0.5 * x * x) / sigma;

	/* ELSE:

	 * x*x  may lose upto about two digits accuracy for "large" x
	 * Morten Welinder's proposal for PR#15620
	 * https://bugs.r-project.org/bugzilla/show_bug.cgi?id=15620

	 * -- 1 --  No hoop jumping when we underflow to zero anyway:

	 *  -x^2/2 <         log(2)*.Machine$double.min.exp  <==>
	 *     x   > sqrt(-2*log(2)*.Machine$double.min.exp) =IEEE= 37.64031
	 * but "thanks" to denormalized numbers, underflow happens a bit later,
	 *  effective.D.MIN.EXP <- with(.Machine, double.min.exp + double.ulp.digits)
	 * for IEEE, DBL_MIN_EXP is -1022 but "effective" is -1074
	 * ==> boundary = sqrt(-2*log(2)*(.Machine$double.min.exp + .Machine$double.ulp.digits))
	 *              =IEEE=  38.58601
	 * [on one x86_64 platform, effective boundary a bit lower: 38.56804]
	 */
	if (x > sqrt(-2*M_LN2*(DBL_MIN_EXP + 1-DBL_MANT_DIG))) return 0.;

	/* Now, to get full accurary, split x into two parts,
	 *  x = x1+x2, such that |x2| <= 2^-16.
	 * Assuming that we are using IEEE doubles, that means that
	 * x1*x1 is error free for x<1024 (but we have x < 38.6 anyway).

	 * If we do not have IEEE this is still an improvement over the naive formula.
	 */
	double x1 = //  R_forceint(x * 65536) / 65536 =
	    ldexp( nearbyint(ldexp(x, 16)), -16);
	double x2 = x - x1;
	return M_1_SQRT_2PI / sigma *
	    (exp(-0.5 * x1 * x1) * exp( (-0.5 * x2 - x1) * x2 ) );
#endif
    }

    /* Probability function */
    template<TAIL_P tail_p=TAIL_LOWER, SCALE_P scale_p=NOR_P>
    double pnorm(const double& xi, const double &mu, const double &sigma) {
	if (isnan(xi) || isnan(mu) || isnan(sigma))
	    RET_NAN("inputs contain NaN, return NaN");
	
	double x, p, cp;
	if (isinf(xi) && mu==xi) RET_NAN("x and mu are infs, return NaN");
	if (sigma <= 0) {
	    if (sigma < 0) RET_NAN("sigma<0, return NaN");
	    if (xi < mu) RET_0; else  RET_1;
	}
	p = (xi - mu) / sigma;
	if (isinf(p))
	    if (xi < mu) RET_0; RET_1; // sigma is very small
	
	x = p;
	pnorm_both<tail_p, scale_p>(x, &p, &cp);
    } // function dnorm
    
#define SIXTEN	16 /* Cutoff allowing exact "*" and "/" */

    /* helper function: pnorm_both
       Compute both lower and upper tail cum 
       if(TAIL_LOWER) return  *cum := P[X <= x]
       if(TAIL_UPPER) return *ccum := P[X >  x] = 1 - P[X <= x]
     */
    template<TAIL_P tail_p, SCALE_P scale_p>  // defaults are defined above
    void pnorm_both(double x, double *cum, double *ccum)
    {
	const static double a[5] = {
	    2.2352520354606839287,
	    161.02823106855587881,
	    1067.6894854603709582,
	    18154.981253343561249,
	    0.065682337918207449113
	};
	const static double b[4] = {
	    47.20258190468824187,
	    976.09855173777669322,
	    10260.932208618978205,
	    45507.789335026729956
	};
	const static double c[9] = {
	    0.39894151208813466764,
	    8.8831497943883759412,
	    93.506656132177855979,
	    597.27027639480026226,
	    2494.5375852903726711,
	    6848.1904505362823326,
	    11602.651437647350124,
	    9842.7148383839780218,
	    1.0765576773720192317e-8
	};
	const static double d[8] = {
	    22.266688044328115691,
	    235.38790178262499861,
	    1519.377599407554805,
	    6485.558298266760755,
	    18615.571640885098091,
	    34900.952721145977266,
	    38912.003286093271411,
	    19685.429676859990727
	};
	const static double p[6] = {
	    0.21589853405795699,
	    0.1274011611602473639,
	    0.022235277870649807,
	    0.001421619193227893466,
	    2.9112874951168792e-5,
	    0.02307344176494017303
	};
	const static double q[5] = {
	    1.28426009614491121,
	    0.468238212480865118,
	    0.0659881378689285515,
	    0.00378239633202758244,
	    7.29751555083966205e-5
	};

	double xden, xnum, temp, del, xsq, y;
	int i;
	if (isnan(x)) {*cum = *ccum = x; return;}

	/* Consider changing these : */
	const double eps = numeric_limits<double>::epsilon() * 0.5;

	/* i_tail in {0,1,2} =^= {lower, upper, both} */
	const bool lower = (tail_p != TAIL_UPPER);
	const bool upper = (tail_p != TAIL_LOWER);

	y = fabs(x);
	if (y <= 0.67448975) { /* qnorm(3/4) = .6744.... -- earlier had 0.66291 */
	    if (y > eps) {
		xsq = x * x;
		xnum = a[4] * xsq;
		xden = xsq;
		for (i = 0; i < 3; ++i) {
		    xnum = (xnum + a[i]) * xsq;
		    xden = (xden + b[i]) * xsq;
		}
	    } else xnum = xden = 0.0;

	    temp = x * (xnum + a[3]) / (xden + b[3]);
	    if(lower)  *cum = 0.5 + temp;
	    if(upper) *ccum = 0.5 - temp;
	    if(scale_p==LOG_P) {
		if(lower)  *cum = log(*cum);
		if(upper) *ccum = log(*ccum);
	    }
	}
	else if (y <= M_SQRT_32) {

	    /* Evaluate pnorm for 0.674.. = qnorm(3/4) < |x| <= sqrt(32) ~= 5.657 */

	    xnum = c[8] * y;
	    xden = y;
	    for (i = 0; i < 7; ++i) {
		xnum = (xnum + c[i]) * y;
		xden = (xden + d[i]) * y;
	    }
	    temp = (xnum + c[7]) / (xden + d[7]);

#define do_del(X)							\
	    xsq = trunc(X * SIXTEN) / SIXTEN;				\
	    del = (X - xsq) * (X + xsq);				\
	    if(scale_p==LOG_P) {							\
		*cum = (-xsq * xsq * 0.5) + (-del * 0.5) + log(temp);	\
		if((lower && x > 0.) || (upper && x <= 0.))		\
		    *ccum = log1p(-exp(-xsq * xsq * 0.5) *		\
				  exp(-del * 0.5) * temp);		\
	    }								\
	    else {							\
		*cum = exp(-xsq * xsq * 0.5) * exp(-del * 0.5) * temp;	\
		*ccum = 1.0 - *cum;					\
	    }

#define swap_tail							\
	    if (x > 0.) {/* swap  ccum <--> cum */			\
		temp = *cum; if(lower) *cum = *ccum; *ccum = temp;	\
	    }

	    do_del(y);
	    swap_tail;
	}
	
/* else	  |x| > sqrt(32) = 5.657 :
 * the next two case differentiations were really for lower=T, log=F
 * Particularly	 *not*	for  log_p !

 * Cody had (-37.5193 < x  &&  x < 8.2924) ; R originally had y < 50
 *
 * Note that we do want symmetry(0), lower/upper -> hence use y
 */
	else if(((scale_p==LOG_P) && y < 1e170) /* avoid underflow below */
		/*  ^^^^^ MM FIXME: can speedup for log_p and much larger |x| !
		 * Then, make use of  Abramowitz & Stegun, 26.2.13, something like

		 xsq = x*x;

		 if(xsq * DBL_EPSILON < 1.)
		 del = (1. - (1. - 5./(xsq+6.)) / (xsq+4.)) / (xsq+2.);
		 else
		 del = 0.;
		 *cum = -.5*xsq - M_LN_SQRT_2PI - log(x) + log1p(-del);
		 *ccum = log1p(-exp(*cum)); /.* ~ log(1) = 0 *./

		 swap_tail;

		 [Yes, but xsq might be infinite.]

		*/
		|| (lower && -37.5193 < x  &&  x < 8.2924)
		|| (upper && -8.2924  < x  &&  x < 37.5193)
	    ) {

	    /* Evaluate pnorm for x in (-37.5, -5.657) union (5.657, 37.5) */
	    xsq = 1.0 / (x * x); /* (1./x)*(1./x) might be better */
	    xnum = p[5] * xsq;
	    xden = xsq;
	    for (i = 0; i < 4; ++i) {
		xnum = (xnum + p[i]) * xsq;
		xden = (xden + q[i]) * xsq;
	    }
	    temp = xsq * (xnum + p[4]) / (xden + q[4]);
	    temp = (M_1_SQRT_2PI - temp) / y;

	    do_del(x);
	    swap_tail;
	} else { /* large x such that probs are 0 or 1 */
	    if(x > 0) { SET_1(*cum); SET_0(*ccum);	}
	    else {	SET_0(*cum); SET_1(*ccum);	}
	}
	
	/* do not return "denormalized" -- we do in R */
	if(scale_p==LOG_P) {
	    if(*cum > -numeric_limits<double>::min())	 *cum = -0.;
	    if(*ccum > -numeric_limits<double>::min())   *ccum = -0.;
	}
	else {
	    if(*cum < numeric_limits<double>::min())	 *cum = 0.;
	    if(*ccum < numeric_limits<double>::min())	 *ccum = 0.;
	}	
    } // function pnorm_both


    /* Quantile function */
    template<TAIL_P tail_p=TAIL_LOWER, SCALE_P scale_p=NOR_P>
    double qnorm(double p, double mu, double sigma)
    {
	double p_, q, r, val;

	if (isnan(p) || isnan(mu) || isnan(sigma))
	    RET_NAN("inputs contain NaN, return NaN");

	R_Q_P01_boundaries(p, -numeric_limits<double>::infinity(),
			   numeric_limits<double>::infinity());

	if(sigma  < 0)	RET_NAN("sigma<0, return NaN");
	if(sigma == 0)	return mu;

	p_ = R_DT_qIv(p);/* real (tail_p==TAIL_LOWER) prob. p */
	q = p_ - 0.5;

/*-- use AS 241 --- */
/* double ppnd16_(double *p, long *ifault)*/
/*      ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3

        Produces the normal deviate Z corresponding to a given lower
        tail area of P; Z is accurate to about 1 part in 10**16.

        (original fortran code used PARAMETER(..) for the coefficients
	and provided hash codes for checking them...)
*/
	if (fabs(q) <= .425) {/* 0.075 <= p <= 0.925 */
	    r = .180625 - q * q;
	    val =
		q * (((((((r * 2509.0809287301226727 +
			   33430.575583588128105) * r + 67265.770927008700853) * r +
			 45921.953931549871457) * r + 13731.693765509461125) * r +
		       1971.5909503065514427) * r + 133.14166789178437745) * r +
		     3.387132872796366608)
		/ (((((((r * 5226.495278852854561 +
			 28729.085735721942674) * r + 39307.89580009271061) * r +
		       21213.794301586595867) * r + 5394.1960214247511077) * r +
		     687.1870074920579083) * r + 42.313330701600911252) * r + 1.);
	}
	else { /* closer than 0.075 from {0,1} boundary */

	    /* r = min(p, 1-p) < 0.075 */
	    if (q > 0)
		r = R_DT_CIv(p);/* 1-p */
	    else
		r = p_;/* = R_DT_Iv(p) ^=  p */

	    r = sqrt(- (((scale_p==LOG_P) &&
			 (((tail_p==TAIL_LOWER) && q <= 0) || (!(tail_p==TAIL_LOWER) && q > 0))) ?
			p : /* else */ log(r)));
	    /* r = sqrt(-log(r))  <==>  min(p, 1-p) = exp( - r^2 ) */

	    if (r <= 5.) { /* <==> min(p,1-p) >= exp(-25) ~= 1.3888e-11 */
		r += -1.6;
		val = (((((((r * 7.7454501427834140764e-4 +
			     .0227238449892691845833) * r + .24178072517745061177) *
			   r + 1.27045825245236838258) * r +
			  3.64784832476320460504) * r + 5.7694972214606914055) *
			r + 4.6303378461565452959) * r +
		       1.42343711074968357734)
		    / (((((((r *
			     1.05075007164441684324e-9 + 5.475938084995344946e-4) *
			    r + .0151986665636164571966) * r +
			   .14810397642748007459) * r + .68976733498510000455) *
			 r + 1.6763848301838038494) * r +
			2.05319162663775882187) * r + 1.);
	    }
	    else { /* very close to  0 or 1 */
		r += -5.;
		val = (((((((r * 2.01033439929228813265e-7 +
			     2.71155556874348757815e-5) * r +
			    .0012426609473880784386) * r + .026532189526576123093) *
			  r + .29656057182850489123) * r +
			 1.7848265399172913358) * r + 5.4637849111641143699) *
		       r + 6.6579046435011037772)
		    / (((((((r *
			     2.04426310338993978564e-15 + 1.4215117583164458887e-7)*
			    r + 1.8463183175100546818e-5) * r +
			   7.868691311456132591e-4) * r + .0148753612908506148525)
			 * r + .13692988092273580531) * r +
			.59983220655588793769) * r + 1.);
	    }

	    if(q < 0.0)
		val = -val;
	    /* return (q >= 0.)? r : -r ;*/
	}
	return mu + sigma * val;	
    } // function qnorm

    template<class RNGType>
    void rnorm(double *x, double mu, double sigma, RNGType *rng)
    {
	if (isnan(mu) || isinf(sigma) || sigma < 0.)
	    SET_NAN_RT(*x, "invalid inputs, return NaN");
	if (sigma == 0. || isinf(mu))
	    *x = mu; /* includes mu = +/- Inf with finite sigma */
	else
	    *x = mu + sigma * rng->normal();
    }

    
} // namespace bnc

#endif /* NORMAL_H */
