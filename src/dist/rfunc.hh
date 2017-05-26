#ifndef RFUNC_H
#define RFUNC_H

#include <type_traits>
#include <limits>
#include <cmath>
#include <cfloat>

#include <matrix/matrix.hh>
#include <util/logger.hh>
#include <common/math.hh>
#include "dconst.hh"

namespace bnc {

#define ML_POSINF	numeric_limits<double>::infinity()
#define ML_NEGINF	-numeric_limits<double>::infinity()
#define ML_NAN		numeric_limits<double>::quiet_NaN()

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

    
    /////////////
    // Uniform //
    /////////////
    template<class RNGType>
    double runif(const double& a, const double& b, RNGType *rng)
    {
	if (!isfinite(a) || !isfinite(b) || b < a)
	{
	    LOG_ERROR("invalid input, return NaN");
	    return numeric_limits<double>::quiet_NaN();
	}

	if (a == b)
	    return a;
	else {
	    double u;
	    /* This is true of all builtin generators, but protect against
	       user-supplied ones */
	    do { u = rng->uniform(); } while (u <= 0 || u >= 1);
	    return a + (b - a) * u;
	}
    }

    template<class RNGType, class AType, class BType>
    Vector runif(const int& n, const AType& a, const BType& b, RNGType *rng)
    {
	Vector ret(n);
	int ai=0, bi=0;
	for (int i=0; i<n; i++) {
	    ret(i) = runif(*(a.data()+ai), *(b.data()+bi), rng);
	    if (++ai >= a.size()) ai = 0;
	    if (++bi >= b.size()) bi = 0;	    
	}
	return ret;
    }    

    template<class RNGType, class BType>
    Vector runif(const int& n, const double& a, const BType& b, RNGType *rng)
    {
	Vector ret(n);
	int bi=0;
	for (int i=0; i<n; i++) {
	    ret(i) = runif(a, *(b.data()+bi), rng);
	    if (++bi >= b.size()) bi = 0;	    
	}
	return ret;
    }    

    template<class RNGType, class AType>
    Vector runif(const int& n, const AType& a, const double& b, RNGType *rng)
    {
	Vector ret(n);
	int ai=0;
	for (int i=0; i<n; i++) {
	    ret(i) = runif(*(a.data()+ai), b, rng);
	    if (++ai >= a.size()) ai = 0;
	}
	return ret;
    }    


    ///////////
    // Gamma //
    ///////////

    template<class RNGType>
    double rgamma(const double &a, const double &scale, RNGType *rng)
    {
/* Constants : */
	const static double sqrt32 = 5.656854;
	const static double exp_m1 = 0.36787944117144232159;/* exp(-1) = 1/e */

	/* Coefficients q[k] - for q0 = sum(q[k]*a^(-k))
	 * Coefficients a[k] - for q = q0+(t*t/2)*sum(a[k]*v^k)
	 * Coefficients e[k] - for exp(q)-1 = sum(e[k]*q^k)
	 */
	const static double q1 = 0.04166669;
	const static double q2 = 0.02083148;
	const static double q3 = 0.00801191;
	const static double q4 = 0.00144121;
	const static double q5 = -7.388e-5;
	const static double q6 = 2.4511e-4;
	const static double q7 = 2.424e-4;

	const static double a1 = 0.3333333;
	const static double a2 = -0.250003;
	const static double a3 = 0.2000062;
	const static double a4 = -0.1662921;
	const static double a5 = 0.1423657;
	const static double a6 = -0.1367177;
	const static double a7 = 0.1233795;

	/* State variables [FIXME for threading!] :*/
	static double aa = 0.;
	static double aaa = 0.;
	static double s, s2, d;    /* no. 1 (step 1) */
	static double q0, b, si, c;/* no. 2 (step 4) */

	double e, p, q, r, t, u, v, w, x, ret_val;

	if (!isfinite(a) || !isfinite(scale) || a < 0.0 || scale <= 0.0) {
	    if(scale == 0.) return 0.;
	    LOG_ERROR("invalid input, return NaN");
	    return numeric_limits<double>::quiet_NaN();
	}

	if (a < 1.) { /* GS algorithm for parameters a < 1 */
	    if(a == 0)
		return 0.;
	    e = 1.0 + exp_m1 * a;
	    for(;;) {
		p = e * rng->uniform();
		if (p >= 1.0) {
		    x = -log((e - p) / a);
		    if (rng->exponential() >= (1.0 - a) * log(x))
			break;
		} else {
		    x = exp(log(p) / a);
		    if (rng->exponential() >= x)
			break;
		}
	    }
	    return scale * x;
	}

	/* --- a >= 1 : GD algorithm --- */

	/* Step 1: Recalculations of s2, s, d if a has changed */
	if (a != aa) {
	    aa = a;
	    s2 = a - 0.5;
	    s = sqrt(s2);
	    d = sqrt32 - s * 12.0;
	}
	/* Step 2: t = standard normal deviate,
	   x = (s,1/2) -normal deviate. */

	/* immediate acceptance (i) */
	t = rng->normal();
	x = s + 0.5 * t;
	ret_val = x * x;
	if (t >= 0.0)
	    return scale * ret_val;

	/* Step 3: u = 0,1 - uniform sample. squeeze acceptance (s) */
	u = rng->uniform();
	if (d * u <= t * t * t)
	    return scale * ret_val;

	/* Step 4: recalculations of q0, b, si, c if necessary */

	if (a != aaa) {
	    aaa = a;
	    r = 1.0 / a;
	    q0 = ((((((q7 * r + q6) * r + q5) * r + q4) * r + q3) * r
		   + q2) * r + q1) * r;

	    /* Approximation depending on size of parameter a */
	    /* The constants in the expressions for b, si and c */
	    /* were established by numerical experiments */

	    if (a <= 3.686) {
		b = 0.463 + s + 0.178 * s2;
		si = 1.235;
		c = 0.195 / s - 0.079 + 0.16 * s;
	    } else if (a <= 13.022) {
		b = 1.654 + 0.0076 * s2;
		si = 1.68 / s + 0.275;
		c = 0.062 / s + 0.024;
	    } else {
		b = 1.77;
		si = 0.75;
		c = 0.1515 / s;
	    }
	}
	/* Step 5: no quotient test if x not positive */

	if (x > 0.0) {
	    /* Step 6: calculation of v and quotient q */
	    v = t / (s + s);
	    if (fabs(v) <= 0.25)
		q = q0 + 0.5 * t * t * ((((((a7 * v + a6) * v + a5) * v + a4) * v
					  + a3) * v + a2) * v + a1) * v;
	    else
		q = q0 - s * t + 0.25 * t * t + (s2 + s2) * log(1.0 + v);


	    /* Step 7: quotient acceptance (q) */
	    if (log(1.0 - u) <= q)
		return scale * ret_val;
	}

	for(;;) {
	    /* Step 8: e = standard exponential deviate
	     *	u =  0,1 -uniform deviate
	     *	t = (b,si)-double exponential (laplace) sample */
	    e = rng->exponential();
	    u = rng->uniform();
	    u = u + u - 1.0;
	    if (u < 0.0)
		t = b - si * e;
	    else
		t = b + si * e;
	    /* Step	 9:  rejection if t < tau(1) = -0.71874483771719 */
	    if (t >= -0.71874483771719) {
		/* Step 10:	 calculation of v and quotient q */
		v = t / (s + s);
		if (fabs(v) <= 0.25)
		    q = q0 + 0.5 * t * t *
			((((((a7 * v + a6) * v + a5) * v + a4) * v + a3) * v
			  + a2) * v + a1) * v;
		else
		    q = q0 - s * t + 0.25 * t * t + (s2 + s2) * log(1.0 + v);
		/* Step 11:	 hat acceptance (h) */
		/* (if q not positive go to step 8) */
		if (q > 0.0) {
		    w = expm1(q);
		    /*  ^^^^^ original code had approximation with rel.err < 2e-7 */
		    /* if t is rejected sample again at step 8 */
		    if (c * fabs(u) <= w * exp(e - 0.5 * t * t))
			break;
		}
	    }
	} /* repeat .. until  `t' is accepted */
	x = s + 0.5 * t;
	return scale * x * x;
    }

    template<class RNGType, class AType, class ScaleType>
    Vector rgamma(const int &n, const AType &a, const ScaleType &scale, RNGType *rng) {
	Vector ret(n);
	int ai=0, bi=0;
	for (int i=0; i<n; i++) {
	    ret(i) = rgamma(*(a.data()+ai), *(scale.data()+bi), rng);
	    if (++ai >= a.size()) ai = 0;
	    if (++bi >= scale.size()) bi = 0;
	}
	return ret;
    }

    template<class RNGType, class ScaleType>
    Vector rgamma(const int &n, const double &a, const ScaleType &scale, RNGType *rng) {
	Vector ret(n);
	int bi=0;
	for (int i=0; i<n; i++) {
	    ret(i) = rgamma(a, *(scale.data()+bi), rng);
	    if (++bi >= scale.size()) bi = 0;	    
	}
	return ret;
    }

    template<class RNGType, class AType>
    Vector rgamma(const int &n, const AType &a, const double &scale, RNGType *rng) {
	Vector ret(n);
	int ai=0;
	for (int i=0; i<n; i++) {
	    ret(i) = rgamma(*(a.data()+ai), scale, rng);
	    if (++ai >= a.size()) ai = 0;
	}
	return ret;
    }


    //////////
    // Beta //
    //////////
    
#define expmax	(DBL_MAX_EXP * M_LN2)/* = log(DBL_MAX) */    
    template <class RNGType>
    double rbeta(const double &aa, const double &bb, RNGType *rng)
    {
	if (aa < 0. || bb < 0.)
	{
	    LOG_ERROR("invalid input, return NaN");
	    return numeric_limits<double>::quiet_NaN();
	}
	if (!R_FINITE(aa) && !R_FINITE(bb)) // a = b = Inf : all mass at 1/2
	    return 0.5;
	if (aa == 0. && bb == 0.) // point mass 1/2 at each of {0,1} :
	    return (unif_rand(rng) < 0.5) ? 0. : 1.;
	// now, at least one of a, b is finite and positive
	if (!R_FINITE(aa) || bb == 0.)
	    return 1.0;
	if (!R_FINITE(bb) || aa == 0.)
	    return 0.0;

	double a, b, alpha;
	double r, s, t, u1, u2, v, w, y, z;
	int qsame;
	/* FIXME:  Keep Globals (properly) for threading */
	/* Uses these GLOBALS to save time when many rv's are generated : */
	static double beta, gamma, delta, k1, k2;
	static double olda = -1.0;
	static double oldb = -1.0;

	/* Test if we need new "initializing" */
	qsame = (olda == aa) && (oldb == bb);
	if (!qsame) { olda = aa; oldb = bb; }

	a = fmin2(aa, bb);
	b = fmax2(aa, bb); /* a <= b */
	alpha = a + b;

#define v_w_from__u1_bet(AA) 			\
	v = beta * log(u1 / (1.0 - u1));	\
	if (v <= expmax) {			\
	    w = AA * exp(v);			\
	    if(!R_FINITE(w)) w = DBL_MAX;	\
	} else					\
	    w = DBL_MAX


	if (a <= 1.0) {	/* --- Algorithm BC --- */

	    /* changed notation, now also a <= b (was reversed) */

	    if (!qsame) { /* initialize */
		beta = 1.0 / a;
		delta = 1.0 + b - a;
		k1 = delta * (0.0138889 + 0.0416667 * a) / (b * beta - 0.777778);
		k2 = 0.25 + (0.5 + 0.25 / delta) * a;
	    }
	    /* FIXME: "do { } while()", but not trivially because of "continue"s:*/
	    for(;;) {
		u1 = unif_rand(rng);
		u2 = unif_rand(rng);
		if (u1 < 0.5) {
		    y = u1 * u2;
		    z = u1 * y;
		    if (0.25 * u2 + z - y >= k1)
			continue;
		} else {
		    z = u1 * u1 * u2;
		    if (z <= 0.25) {
			v_w_from__u1_bet(b);
			break;
		    }
		    if (z >= k2)
			continue;
		}

		v_w_from__u1_bet(b);

		if (alpha * (log(alpha / (a + w)) + v) - 1.3862944 >= log(z))
		    break;
	    }
	    return (aa == a) ? a / (a + w) : w / (a + w);

	}
	else {		/* Algorithm BB */

	    if (!qsame) { /* initialize */
		beta = sqrt((alpha - 2.0) / (2.0 * a * b - alpha));
		gamma = a + 1.0 / beta;
	    }
	    do {
		u1 = unif_rand(rng);
		u2 = unif_rand(rng);

		v_w_from__u1_bet(a);

		z = u1 * u1 * u2;
		r = gamma * v - 1.3862944;
		s = a + r - w;
		if (s + 2.609438 >= 5.0 * z)
		    break;
		t = log(z);
		if (s > t)
		    break;
	    }
	    while (r + alpha * log(alpha / (b + w)) < t);

	    return (aa != a) ? b / (b + w) : w / (b + w);
	}
    }
#undef    expmax
#undef    v_w_from__u1_bet

    template <class RNGType, class AType, class BType>
    Vector rbeta(const int &n, const AType &aa, const BType &bb, RNGType *rng)
    {
	Vector ret(n);
	int ai=0, bi=0;
	for (int i=0; i<n; i++) {
	    ret(i) = rbeta(*(aa.data()+ai), *(bb.data()+bi), rng);
	    if (++ai >= aa.size()) ai = 0;
	    if (++bi >= bb.size()) bi = 0;
	}
	return ret;
    }

    template <class RNGType, class BType>
    Vector rbeta(const int &n, const double &aa, const BType &bb, RNGType *rng)
    {
	Vector ret(n);
	int bi=0;
	for (int i=0; i<n; i++) {
	    ret(i) = rbeta(aa, *(bb.data()+bi), rng);
	    if (++bi >= bb.size()) bi = 0;
	}
	return ret;
    }

    template <class RNGType, class AType>
    Vector rbeta(const int &n, const AType &aa, const double &bb, RNGType *rng)
    {
	Vector ret(n);
	int ai=0;
	for (int i=0; i<n; i++) {
	    ret(i) = rbeta(*(aa.data()+ai), bb, rng);
	    if (++ai >= aa.size()) ai = 0;
	}
	return ret;
    }

    
    //////////////
    // Binomial //
    //////////////

    template<class RNGType>
    double rbinom(const double &nin, const double &pp, RNGType *rng)
    {
	/* FIXME: These should become THREAD_specific globals : */

	static double c, fm, npq, p1, p2, p3, p4, qn;
	static double xl, xll, xlr, xm, xr;

	static double psave = -1.0;
	static int nsave = -1;
	static int m;

	double f, f1, f2, u, v, w, w2, x, x1, x2, z, z2;
	double p, q, np, g, r, al, alv, amaxp, ffm, ynorm;
	int i, ix, k, n;

	if (!R_FINITE(nin))
	{
	    LOG_ERROR("invalid input, return NaN");
	    return numeric_limits<double>::quiet_NaN();
	}

	r = R_forceint(nin);
	if (r != nin)
	{
	    LOG_ERROR("invalid input, return NaN");
	    return numeric_limits<double>::quiet_NaN();
	}

	if (!R_FINITE(pp) ||
	    /* n=0, p=0, p=1 are not errors <TSL>*/
	    r < 0 || pp < 0. || pp > 1.)
	{
	    LOG_ERROR("invalid input, return NaN");
	    return numeric_limits<double>::quiet_NaN();
	}

	if (r == 0 || pp == 0.) return 0;
	if (pp == 1.) return r;

	if (r >= INT_MAX)/* evade integer overflow,
			    and r == INT_MAX gave only even values */
	    return qbinom(unif_rand(rng), r, pp, /*lower_tail*/ 0, /*log_p*/ 0);
	/* else */
	n = (int) r;

	p = fmin2(pp, 1. - pp);
	q = 1. - p;
	np = n * p;
	r = p / q;
	g = r * (n + 1);

	/* Setup, perform only when parameters change [using static (globals): */

	/* FIXING: Want this thread safe
	   -- use as little (thread globals) as possible
	*/
	if (pp != psave || n != nsave) {
	    psave = pp;
	    nsave = n;
	    if (np < 30.0) {
		/* inverse cdf logic for mean less than 30 */
		qn = JR_pow_di(q, n);
		goto L_np_small;
	    } else {
		ffm = np + p;
		m = (int) ffm;
		fm = m;
		npq = np * q;
		p1 = (int)(2.195 * sqrt(npq) - 4.6 * q) + 0.5;
		xm = fm + 0.5;
		xl = xm - p1;
		xr = xm + p1;
		c = 0.134 + 20.5 / (15.3 + fm);
		al = (ffm - xl) / (ffm - xl * p);
		xll = al * (1.0 + 0.5 * al);
		al = (xr - ffm) / (xr * q);
		xlr = al * (1.0 + 0.5 * al);
		p2 = p1 * (1.0 + c + c);
		p3 = p2 + c / xll;
		p4 = p3 + c / xlr;
	    }
	} else if (n == nsave) {
	    if (np < 30.0)
		goto L_np_small;
	}

	/*-------------------------- np = n*p >= 30 : ------------------- */
	for(;;) {
	    u = unif_rand(rng) * p4;
	    v = unif_rand(rng);
	    /* triangular region */
	    if (u <= p1) {
		ix = (int)(xm - p1 * v + u);
		goto finis;
	    }
	    /* parallelogram region */
	    if (u <= p2) {
		x = xl + (u - p1) / c;
		v = v * c + 1.0 - fabs(xm - x) / p1;
		if (v > 1.0 || v <= 0.)
		    continue;
		ix = (int) x;
	    } else {
		if (u > p3) {	/* right tail */
		    ix = (int)(xr - log(v) / xlr);
		    if (ix > n)
			continue;
		    v = v * (u - p3) * xlr;
		} else {/* left tail */
		    ix = (int)(xl + log(v) / xll);
		    if (ix < 0)
			continue;
		    v = v * (u - p2) * xll;
		}
	    }
	    /* determine appropriate way to perform accept/reject test */
	    k = abs(ix - m);
	    if (k <= 20 || k >= npq / 2 - 1) {
		/* explicit evaluation */
		f = 1.0;
		if (m < ix) {
		    for (i = m + 1; i <= ix; i++)
			f *= (g / i - r);
		} else if (m != ix) {
		    for (i = ix + 1; i <= m; i++)
			f /= (g / i - r);
		}
		if (v <= f)
		    goto finis;
	    } else {
		/* squeezing using upper and lower bounds on log(f(x)) */
		amaxp = (k / npq) * ((k * (k / 3. + 0.625) + 0.1666666666666) / npq + 0.5);
		ynorm = -k * k / (2.0 * npq);
		alv = log(v);
		if (alv < ynorm - amaxp)
		    goto finis;
		if (alv <= ynorm + amaxp) {
		    /* stirling's formula to machine accuracy */
		    /* for the final acceptance/rejection test */
		    x1 = ix + 1;
		    f1 = fm + 1.0;
		    z = n + 1 - fm;
		    w = n - ix + 1.0;
		    z2 = z * z;
		    x2 = x1 * x1;
		    f2 = f1 * f1;
		    w2 = w * w;
		    if (alv <= xm * log(f1 / x1) + (n - m + 0.5) * log(z / w) + (ix - m) * log(w * p / (x1 * q)) + (13860.0 - (462.0 - (132.0 - (99.0 - 140.0 / f2) / f2) / f2) / f2) / f1 / 166320.0 + (13860.0 - (462.0 - (132.0 - (99.0 - 140.0 / z2) / z2) / z2) / z2) / z / 166320.0 + (13860.0 - (462.0 - (132.0 - (99.0 - 140.0 / x2) / x2) / x2) / x2) / x1 / 166320.0 + (13860.0 - (462.0 - (132.0 - (99.0 - 140.0 / w2) / w2) / w2) / w2) / w / 166320.)
			goto finis;
		}
	    }
	}

    L_np_small:
	/*---------------------- np = n*p < 30 : ------------------------- */

	for(;;) {
	    ix = 0;
	    f = qn;
	    u = unif_rand(rng);
	    for(;;) {
		if (u < f)
		    goto finis;
		if (ix > 110)
		    break;
		u -= f;
		ix++;
		f *= (g / ix - r);
	    }
	}
    finis:
	if (psave > 0.5)
	    ix = n - ix;
	return (double)ix;
    }

    template<class RNGType>
    Vector rbinom(const int &n, const double &nin, const double &pp, RNGType *rng)
    {
	Vector ret(n);
	int ai=0, bi=0;
	for (int i=0; i<n; i++) {
	    ret(i) = rbinom(*(nin.data()+ai), *(pp.data()+bi), rng);
	    if (++ai >= nin.size()) ai = 0;	    
	    if (++bi >= pp.size()) bi = 0;
	}
	return ret;
    }

    template<class RNGType, class AType, class BType>
    Vector rbinom(const int &n, const AType &nin, const BType &pp, RNGType *rng)
    {
	Vector ret(n);
	int ai=0, bi=0;
	for (int i=0; i<n; i++) {
	    ret(i) = rbinom(*(nin.data()+ai), *(pp.data()+bi), rng);
	    if (++ai >= nin.size()) ai = 0;	    
	    if (++bi >= pp.size()) bi = 0;
	}
	return ret;
    }

    template<class RNGType, class BType>
    Vector rbinom(const int &n, const double &nin, const BType &pp, RNGType *rng)
    {
	Vector ret(n);
	int bi=0;
	for (int i=0; i<n; i++) {
	    ret(i) = rbinom(nin, *(pp.data()+bi), rng);
	    if (++bi >= pp.size()) bi = 0;
	}
	return ret;
    }

    
#undef ML_POSINF    
#undef ML_NEGINF    
#undef ML_NAN

}  // namespace bnc

#endif /* RFUNC_H */
