#ifndef POIS_H
#define POIS_H

#include <dist/brutil.hh>

extern "C" {
    double br_ppois(double x, double lambda, int lower_tail, int log_p);
    double br_dpois(double x, double lambda, int give_log);
    double br_qpois(double p, double lambda, int lower_tail, int log_p);
}

namespace bnc {
    
#define a0	-0.5
#define a1	 0.3333333
#define a2	-0.2500068
#define a3	 0.2000118
#define a4	-0.1661269
#define a5	 0.1421878
#define a6	-0.1384794
#define a7	 0.1250060

#define one_7	0.1428571428571428571
#define one_12	0.0833333333333333333
#define one_24	0.0416666666666666667

#define repeat for(;;)

    template<class RNGType>
    double rpois(double mu, RNGType *rng)
    {
	/* Factorial Table (0:9)! */
	const static double fact[10] =
	    {
		1., 1., 2., 6., 24., 120., 720., 5040., 40320., 362880.
	    };

	/* These are static --- persistent between calls for same mu : */
	static int l, m;

	static double b1, b2, c, c0, c1, c2, c3;
	static double pp[36], p0, p, q, s, d, omega;
	static double big_l;/* integer "w/o overflow" */
	static double muprev = 0., muprev2 = 0.;/*, muold	 = 0.*/

	/* Local Vars  [initialize some for -Wall]: */
	double del, difmuk= 0., E= 0., fk= 0., fx, fy, g, px, py, t, u= 0., v, x;
	double pois = -1.;
	int k, kflag, big_mu, new_big_mu = false;

	if (!R_FINITE(mu) || mu < 0)
	    ML_ERR_return_NAN;

	if (mu <= 0.)
	    return 0.;

	big_mu = mu >= 10.;
	if(big_mu)
	    new_big_mu = false;

	if (!(big_mu && mu == muprev)) {/* maybe compute new persistent par.s */

	    if (big_mu) {
		new_big_mu = true;
		/* Case A. (recalculation of s,d,l	because mu has changed):
		 * The poisson probabilities pk exceed the discrete normal
		 * probabilities fk whenever k >= m(mu).
		 */
		muprev = mu;
		s = sqrt(mu);
		d = 6. * mu * mu;
		big_l = floor(mu - 1.1484);
		/* = an upper bound to m(mu) for all mu >= 10.*/
	    }
	    else { /* Small mu ( < 10) -- not using normal approx. */

		/* Case B. (start new table and calculate p0 if necessary) */

		/*muprev = 0.;-* such that next time, mu != muprev ..*/
		if (mu != muprev) {
		    muprev = mu;
		    m = imax2(1, (int) mu);
		    l = 0; /* pp[] is already ok up to pp[l] */
		    q = p0 = p = exp(-mu);
		}

		repeat {
		    /* Step U. uniform sample for inversion method */
		    u = unif_rand(rng);
		    if (u <= p0)
			return 0.;

		    /* Step T. table comparison until the end pp[l] of the
		       pp-table of cumulative poisson probabilities
		       (0.458 > ~= pp[9](= 0.45792971447) for mu=10 ) */
		    if (l != 0) {
			for (k = (u <= 0.458) ? 1 : imin2(l, m);  k <= l; k++)
			    if (u <= pp[k])
				return (double)k;
			if (l == 35) /* u > pp[35] */
			    continue;
		    }
		    /* Step C. creation of new poisson
		       probabilities p[l..] and their cumulatives q =: pp[k] */
		    l++;
		    for (k = l; k <= 35; k++) {
			p *= mu / k;
			q += p;
			pp[k] = q;
			if (u <= q) {
			    l = k;
			    return (double)k;
			}
		    }
		    l = 35;
		} /* end(repeat) */
	    }/* mu < 10 */

	} /* end {initialize persistent vars} */

/* Only if mu >= 10 : ----------------------- */

	/* Step N. normal sample */
	g = mu + s * norm_rand(rng);/* norm_rand(rng) ~ N(0,1), standard normal */

	if (g >= 0.) {
	    pois = floor(g);
	    /* Step I. immediate acceptance if pois is large enough */
	    if (pois >= big_l)
		return pois;
	    /* Step S. squeeze acceptance */
	    fk = pois;
	    difmuk = mu - fk;
	    u = unif_rand(rng); /* ~ U(0,1) - sample */
	    if (d * u >= difmuk * difmuk * difmuk)
		return pois;
	}

	/* Step P. preparations for steps Q and H.
	   (recalculations of parameters if necessary) */

	if (new_big_mu || mu != muprev2) {
	    /* Careful! muprev2 is not always == muprev
	       because one might have exited in step I or S
	    */
	    muprev2 = mu;
	    omega = M_1_SQRT_2PI / s;
	    /* The quantities b1, b2, c3, c2, c1, c0 are for the Hermite
	     * approximations to the discrete normal probabilities fk. */

	    b1 = one_24 / mu;
	    b2 = 0.3 * b1 * b1;
	    c3 = one_7 * b1 * b2;
	    c2 = b2 - 15. * c3;
	    c1 = b1 - 6. * b2 + 45. * c3;
	    c0 = 1. - b1 + 3. * b2 - 15. * c3;
	    c = 0.1069 / mu; /* guarantees majorization by the 'hat'-function. */
	}

	if (g >= 0.) {
	    /* 'Subroutine' F is called (kflag=0 for correct return) */
	    kflag = 0;
	    goto Step_F;
	}


	repeat {
	    /* Step E. Exponential Sample */

	    E = exp_rand(rng);	/* ~ Exp(1) (standard exponential) */

	    /*  sample t from the laplace 'hat'
		(if t <= -0.6744 then pk < fk for all mu >= 10.) */
	    u = 2 * unif_rand(rng) - 1.;
	    t = 1.8 + fsign(E, u);
	    if (t > -0.6744) {
		pois = floor(mu + s * t);
		fk = pois;
		difmuk = mu - fk;

		/* 'subroutine' F is called (kflag=1 for correct return) */
		kflag = 1;

	    Step_F: /* 'subroutine' F : calculation of px,py,fx,fy. */

		if (pois < 10) { /* use factorials from table fact[] */
		    px = -mu;
		    py = pow(mu, pois) / fact[(int)pois];
		}
		else {
		    /* Case pois >= 10 uses polynomial approximation
		       a0-a7 for accuracy when advisable */
		    del = one_12 / fk;
		    del = del * (1. - 4.8 * del * del);
		    v = difmuk / fk;
		    if (fabs(v) <= 0.25)
			px = fk * v * v * (((((((a7 * v + a6) * v + a5) * v + a4) *
					      v + a3) * v + a2) * v + a1) * v + a0)
			    - del;
		    else /* |v| > 1/4 */
			px = fk * log(1. + v) - difmuk - del;
		    py = M_1_SQRT_2PI / sqrt(fk);
		}
		x = (0.5 - difmuk) / s;
		x *= x;/* x^2 */
		fx = -0.5 * x;
		fy = omega * (((c3 * x + c2) * x + c1) * x + c0);
		if (kflag > 0) {
		    /* Step H. Hat acceptance (E is repeated on rejection) */
		    if (c * fabs(u) <= py * exp(px + E) - fy * exp(fx + E))
			break;
		} else
		    /* Step Q. Quotient acceptance (rare case) */
		    if (fy - u * fy <= py * exp(px - fx))
			break;
	    }/* t > -.67.. */
	}
	return pois;
    }
    
#undef a0
#undef a1
#undef a2
#undef a3
#undef a4
#undef a5
#undef a6
#undef a7

#undef one_7
#undef one_12
#undef one_24

#undef repeat

    R_RFUNC_INTERFACE_1ARG(rpois);

    /* D functions */
    R_DFUNC_INTERFACE_3ARG(dpois, br_dpois);

    /* P functions */
    R_PFUNC_INTERFACE_4ARG(ppois, br_ppois);

    /* Q functions */
    R_QFUNC_INTERFACE_4ARG(qpois, br_qpois);
    
    
}  // namespace bnc

#endif /* POIS_H */
