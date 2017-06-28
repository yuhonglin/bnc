/*
 * Modified from the mvt.f in mvtnorm function.
 * Provide interface:
 *   - ```mvt```
 *
 * The code is based on the result of f2c
 * To enable multithread, I have tried to remove the ```static```
 * declaration for non-constant local variables
 * 
 */


#ifndef MVT_H
#define MVT_H

#include <algorithm>
#include <cmath>
#include <memory>
#include <iostream>

#include <util/logger.hh>
#include <util/constant.hh>
#include <matrix/matrix.hh>

extern "C" {
    double br_qchisq(double, double, int, int);
}

namespace bnc {
    /* 
     *  I decide to wrap the mvt functions into a 
     *  class because there lots of static arrays used
     *  in these functions. This will make multi thread
     *  not available.
     */
    class Mvt {
	// The member variables are prefixed by the functions that uses it.

	// function mvsubr_0_
	double *mvsubr_0_a;
	double *mvsubr_0_b;
	double mvsubr_0_r__;
	double *mvsubr_0_y;
	double mvsubr_0_di;
	double mvsubr_0_ei;
	double *mvsubr_0_dl;
	double *mvsubr_0_cov;
	int *mvsubr_0_infi;
	int mvsubr_0_nu;
	int mvsubr_0_ny;
	double mvsubr_0_snu;

	// function mvspcl_
	double mvspcl_r__;

	// function mvvlsb_

	// function mvsort_

	// function mvkbrv_
	double *mvkbrv_r__;
	double *mvkbrv_x;
	double *mvkbrv_fs;
	int mvkbrv_np;
	double *mvkbrv_vk;
	int *mvkbrv_pr;
	double *mvkbrv_finval;
	int mvkbrv_sampls;
	double *mvkbrv_values;
	double *mvkbrv_varest;
	double *mvkbrv_varsqr;
	
    public:
	friend int mvsubr_(int *n, double *w, int *nf, 
			   double *f);
	Mvt() {
	    mvsubr_0_a = new double[1000];
	    mvsubr_0_b = new double[1000];
	    mvsubr_0_y = new double[1000];
	    mvsubr_0_dl = new double[1000];
	    mvsubr_0_cov = new double[500500];
	    mvsubr_0_infi = new int[1000];

	    mvkbrv_r__ = new double[1000];
	    mvkbrv_x   = new double[1000];
	    mvkbrv_fs  = new double[5000];
	    mvkbrv_vk  = new double[1000];
	    mvkbrv_pr  = new int[1000];
	    mvkbrv_finval = new double[5000];
	    mvkbrv_values = new double[5000];
	    mvkbrv_varest = new double[5000];
	    mvkbrv_varsqr = new double[5000];
	}

	template <class T>
	inline void safedelete(T* p) {
	    if (p!=NULL)
		delete[] p;
	}
	
	~Mvt() {
	    safedelete(mvsubr_0_a);
	    safedelete(mvsubr_0_b);
	    safedelete(mvsubr_0_y);
	    safedelete(mvsubr_0_dl);
	    safedelete(mvsubr_0_cov);
	    safedelete(mvsubr_0_infi);
	}
	
	
	// src: www.netlib.org/clapack/CLAPACK-3.1.1/F2CLIBS/libf2c/pow_dd.c
	double pow_dd(double *ap, double *bp)
	{
	    return(std::pow(*ap, *bp) );
	}
    
	// src: www.netlib.org/clapack/CLAPACK-3.1.1/F2CLIBS/libf2c/d_mod.c
	double d_mod(double *x, double *y)
	{
	    double quotient;
	    if( (quotient = *x / *y) >= 0)
		quotient = std::floor(quotient);
	    else
		quotient = -std::floor(-quotient);
	    return(*x - (*y) * quotient );
	}
    
	// src: www.netlib.org/clapack/CLAPACK-3.1.1/F2CLIBS/libf2c/d_sign.c
	double d_sign(double *a, double *b)
	{
	    double x;
	    x = (*a >= 0 ? *a : - *a);
	    return( *b >= 0 ? x : -x);
	}
    
	// src: www.netlib.org/clapack/CLAPACK-3.1.1/F2CLIBS/libf2c/pow_di.c
	double pow_di(double *ap, int *bp)
	{
	    double pow, x;
	    int n;
	    unsigned long u;

	    pow = 1;
	    x = *ap;
	    n = *bp;

	    if(n != 0)
	    {
		if(n < 0)
		{
		    n = -n;
		    x = 1/x;
		}
		for(u = n; ; )
		{
		    if(u & 01)
			pow *= x;
		    if(u >>= 1)
			x *= x;
		    else
			break;
		}
	    }
	    return(pow);
	}
    
	double sqrtqchisqint_(int *n, double *p) {
	    return(std::sqrt(::br_qchisq(p[0], (double) n[0], 0, 0)));
	}
    
	/* Common Block Declarations */
    
	int ivls;

	/* Table of constant values */

	int c__1 = 1;
	int c_true = true;
	int c__0 = 0;
	double c_b24 = 1.;
	double c_b29 = 2.;
	int c__100 = 100;

	/* Subroutine */ int mvsubr_(int *n, double *w, int *nf, 
					  double *f)
	{
	    return mvsubr_0_(0, n, w, nf, f, (int *)0, (double *)0, (
				 double *)0, (double *)0, (double *)0, (int *)0, (
				     int *)0, (double *)0, (double *)0, (int *)0);
	}

	/*    $Id: mvt.f 297 2014-12-16 17:24:38Z thothorn $ */

	template <class RNGType>
	/* Subroutine */ int mvtdst_(int *n, int *nu, double *lower, 
					  double *upper, int *infin, double *correl, double *
					  delta, int *maxpts, double *abseps, double *releps, 
					  double *error, double *value, int *inform__, RNGType *rng)
	{
	    double e[1], v[1];
	    int nd;

	    /*     A subroutine for computing non-central multivariate t probabilities. */
	    /*     This subroutine uses an algorithm (QRSVN) described in the paper */
	    /*     "Comparison of Methods for the Computation of Multivariate */
	    /*         t-Probabilities", by Alan Genz and Frank Bretz */
	    /*         J. Comp. Graph. Stat. 11 (2002), pp. 950-971. */

	    /*          Alan Genz */
	    /*          Department of Mathematics */
	    /*          Washington State University */
	    /*          Pullman, WA 99164-3113 */
	    /*          Email : AlanGenz@wsu.edu */

	    /* 	Original source available from */
	    /* 	http://www.math.wsu.edu/faculty/genz/software/fort77/mvtdstpack.f */

	    /* 	This is version 28/05/2013 */

	    /*  Parameters */

	    /*     N      INT, the number of variables. */
	    /*     NU     INT, the number of degrees of freedom. */
	    /*            If NU < 1, then an MVN probability is computed. */
	    /*     LOWER  DOUBLE PRECISION, array of lower integration limits. */
	    /*     UPPER  DOUBLE PRECISION, array of upper integration limits. */
	    /*     INFIN  INT, array of integration limits flags: */
	    /*             if INFIN(I) < 0, Ith limits are (-infinity, infinity); */
	    /*             if INFIN(I) = 0, Ith limits are (-infinity, UPPER(I)]; */
	    /*             if INFIN(I) = 1, Ith limits are [LOWER(I), infinity); */
	    /*             if INFIN(I) = 2, Ith limits are [LOWER(I), UPPER(I)]. */
	    /*     CORREL DOUBLE PRECISION, array of correlation coefficients; */
	    /*            the correlation coefficient in row I column J of the */
	    /*            correlation matrixshould be stored in */
	    /*               CORREL( J + ((I-2)*(I-1))/2 ), for J < I. */
	    /*            The correlation matrix must be positive semi-definite. */
	    /*     DELTA  DOUBLE PRECISION, array of non-centrality parameters. */
	    /*     MAXPTS INT, maximum number of function values allowed. This */
	    /*            parameter can be used to limit the time. A sensible */
	    /*            strategy is to start with MAXPTS = 1000*N, and then */
	    /*            increase MAXPTS if ERROR is too large. */
	    /*     ABSEPS DOUBLE PRECISION absolute error tolerance. */
	    /*     RELEPS DOUBLE PRECISION relative error tolerance. */
	    /*     ERROR  DOUBLE PRECISION estimated absolute error, */
	    /*            with 99% confidence level. */
	    /*     VALUE  DOUBLE PRECISION estimated value for the integral */
	    /*     INFORM INT, termination status parameter: */
	    /*            if INFORM = 0, normal completion with ERROR < EPS; */
	    /*            if INFORM = 1, completion with ERROR > EPS and MAXPTS */
	    /*                           function vaules used; increase MAXPTS to */
	    /*                           decrease ERROR; */
	    /*            if INFORM = 2, N > 1000 or N < 1. */
	    /*            if INFORM = 3, correlation matrix not positive semi-definite. */

	    /* Parameter adjustments */
	    --delta;
	    --correl;
	    --infin;
	    --upper;
	    --lower;
      
	    /* Function Body */
	    ivls = 0;
	    if (*n > 1000 || *n < 1) {
		*value = 0.;
		*error = 1.;
		*inform__ = 2;
	    } else {
		mvints_(n, nu, &correl[1], &lower[1], &upper[1], &delta[1], &infin[1],
			&nd, value, error, inform__);
		if (*inform__ == 0 && nd > 0) {

		    /*           Call the lattice rule integration subroutine */
	  
		    mvkbrv_(&nd, &ivls, maxpts, &c__1, abseps,
			    releps, e, v, inform__, rng);
		    *error = e[0];
		    *value = v[0];
		}
	    }
	    return 0;
	} /* mvtdst_ */


	int mvsubr_0_(int n__, int *n, double *w, int *
			   nf, double *f, int *nuin, double *correl, double *
			   lower, double *upper, double *delta, int *infin, int *
			   nd, double *vl, double *er, int *inform__)
	{
	    /* System generated locals */
	    int i__1;

	    /* Local variables */

	    /*     Integrand subroutine */

	    /* Parameter adjustments */
	    if (w) {
		--w;
	    }
	    if (f) {
		--f;
	    }
	    if (correl) {
		--correl;
	    }
	    if (lower) {
		--lower;
	    }
	    if (upper) {
		--upper;
	    }
	    if (delta) {
		--delta;
	    }
	    if (infin) {
		--infin;
	    }

	    /* Function Body */
	    switch(n__) {
	    case 1: goto L_mvints;
	    }

	    if (mvsubr_0_nu <= 0) {
		mvsubr_0_r__ = 1.;
		i__1 = *n + 1;
		mvvlsb_(&i__1, &w[1], &mvsubr_0_r__, mvsubr_0_dl, mvsubr_0_infi, mvsubr_0_a, mvsubr_0_b, mvsubr_0_cov, mvsubr_0_y, &mvsubr_0_di, &mvsubr_0_ei, &mvsubr_0_ny, &f[
			    1]);
	    } else {
		mvsubr_0_r__ = mvchnv_(&mvsubr_0_nu, &w[*n]) / mvsubr_0_snu;
		mvvlsb_(n, &w[1], &mvsubr_0_r__, mvsubr_0_dl, mvsubr_0_infi, mvsubr_0_a, mvsubr_0_b, mvsubr_0_cov, mvsubr_0_y, &mvsubr_0_di, &mvsubr_0_ei, &mvsubr_0_ny, &f[1]);
	    }
	    return 0;

	    /*     Entry point for intialization. */


	L_mvints:

	    /*     Initialization and computation of covariance Cholesky factor. */

	    mvsort_(n, &lower[1], &upper[1], &delta[1], &correl[1], &infin[1], mvsubr_0_y, &
		    c_true, nd, mvsubr_0_a, mvsubr_0_b, mvsubr_0_dl, mvsubr_0_cov, mvsubr_0_infi, inform__);
	    mvsubr_0_nu = *nuin;
	    mvspcl_(nd, &mvsubr_0_nu, mvsubr_0_a, mvsubr_0_b, mvsubr_0_dl, mvsubr_0_cov, mvsubr_0_infi, &mvsubr_0_snu, vl, er, inform__);
	    return 0;
	} /* mvsubr_ */

	/* Subroutine */ int mvints_(int *n, int *nuin, double *correl, 
					  double *lower, double *upper, double *delta, int *
					  infin, int *nd, double *vl, double *er, int *inform__)
	{
	    return mvsubr_0_(1, n, (double *)0, (int *)0, (double *)0, 
			     nuin, correl, lower, upper, delta, infin, nd, vl, er, inform__);
	}


	int mvspcl_(int *nd, int *nu, double *a, 
			 double *b, double *dl, double *cov, int *infi, 
			 double *snu, double *vl, double *er, int *inform__)
	{
	    /* System generated locals */
	    double d__1;

	    /* Local variables */

	    /*     Special cases subroutine */

	    /* Parameter adjustments */
	    --infi;
	    --cov;
	    --dl;
	    --b;
	    --a;

	    /* Function Body */
	    if (*inform__ > 0) {
		*vl = 0.;
		*er = 1.;
	    } else {

		/*        Special cases */

		if (*nd == 0) {
		    *er = 0.;
		    /*  Code added to fix ND = 0 bug, 24/03/2009 -> */
		    *vl = 1.;
		    /*  <- Code added to fix ND = 0 bug, 24/03/2009 */
		} else if (*nd == 1 && (*nu < 1 || std::abs(dl[1]) == 0.)) {

		    /*           1-d case for normal or central t */

		    *vl = 1.;
		    if (infi[1] != 1) {
			d__1 = b[1] - dl[1];
			*vl = mvstdt_(nu, &d__1);
		    }
		    if (infi[1] != 0) {
			d__1 = a[1] - dl[1];
			*vl -= mvstdt_(nu, &d__1);
		    }
		    if (*vl < 0.) {
			*vl = 0.;
		    }
		    *er = 2e-16;
		    *nd = 0;
		} else if (*nd == 2 && (*nu < 1 || std::abs(dl[1]) + std::abs(dl[2]) == 0.)) {

		    /*           2-d case for normal or central t */

		    if (infi[1] != 0) {
			a[1] -= dl[1];
		    }
		    if (infi[1] != 1) {
			b[1] -= dl[1];
		    }
		    if (infi[2] != 0) {
			a[2] -= dl[2];
		    }
		    if (infi[2] != 1) {
			b[2] -= dl[2];
		    }
		    if (std::abs(cov[3]) > 0.) {

			/*              2-d nonsingular case */

			/* Computing 2nd power */
			d__1 = cov[2];
			mvspcl_r__ = std::sqrt(d__1 * d__1 + 1);
			if (infi[2] != 0) {
			    a[2] /= mvspcl_r__;
			}
			if (infi[2] != 1) {
			    b[2] /= mvspcl_r__;
			}
			cov[2] /= mvspcl_r__;
			*vl = mvbvt_(nu, &a[1], &b[1], &infi[1], &cov[2]);
			*er = 1e-15;
		    } else {

			/*              2-d singular case */

			if (infi[1] != 0) {
			    if (infi[2] != 0) {
				a[1] = std::max(a[1],a[2]);
			    }
			} else {
			    if (infi[2] != 0) {
				a[1] = a[2];
			    }
			}
			if (infi[1] != 1) {
			    if (infi[2] != 1) {
				b[1] = std::min(b[1],b[2]);
			    }
			} else {
			    if (infi[2] != 1) {
				b[1] = b[2];
			    }
			}
			if (infi[1] != infi[2]) {
			    infi[1] = 2;
			}
			*vl = 1.;
			/*  A(1), B(1) Bug Fixed, 28/05/2013 */
			if (infi[1] != 1) {
			    *vl = mvstdt_(nu, &b[1]);
			}
			if (infi[1] != 0) {
			    *vl -= mvstdt_(nu, &a[1]);
			}
			if (*vl < 0.) {
			    *vl = 0.;
			}
			*er = 2e-16;
		    }
		    *nd = 0;
		} else {
		    if (*nu > 0) {
			*snu = std::sqrt((double) (*nu));
		    } else {
			--(*nd);
		    }
		}
	    }
	    return 0;
	} /* mvspcl_ */


	int mvvlsb_(int *n, double *w, double *r__, 
			 double *dl, int *infi, double *a, double *b, 
			 double *cov, double *y, double *di, double *ei, 
			 int *nd, double *value)
	{
	    /* System generated locals */
	    int i__1, i__2;
	    double d__1, d__2;

	    /* Local variables */
	    int mvvlsb_i__;
	    int mvvlsb_j;
	    double mvvlsb_ai;
	    double mvvlsb_bi;
	    int mvvlsb_ij;
	    double mvvlsb_sum;
	    int mvvlsb_infa;
	    int mvvlsb_infb;

	    /*     Integrand subroutine */

	    /* Parameter adjustments */
	    --y;
	    --cov;
	    --b;
	    --a;
	    --infi;
	    --dl;
	    --w;

	    /* Function Body */
	    *value = 1.;
	    mvvlsb_infa = 0;
	    mvvlsb_infb = 0;
	    *nd = 0;
	    mvvlsb_ij = 0;
	    i__1 = *n;
	    for (mvvlsb_i__ = 1; mvvlsb_i__ <= i__1; ++mvvlsb_i__) {
		mvvlsb_sum = dl[mvvlsb_i__];
		i__2 = mvvlsb_i__ - 1;
		for (mvvlsb_j = 1; mvvlsb_j <= i__2; ++mvvlsb_j) {
		    ++mvvlsb_ij;
		    if (mvvlsb_j <= *nd) {
			mvvlsb_sum += cov[mvvlsb_ij] * y[mvvlsb_j];
		    }
		}
		if (infi[mvvlsb_i__] != 0) {
		    if (mvvlsb_infa == 1) {
			/* Computing MAX */
			d__1 = mvvlsb_ai, d__2 = *r__ * a[mvvlsb_i__] - mvvlsb_sum;
			mvvlsb_ai = std::max(d__1,d__2);
		    } else {
			mvvlsb_ai = *r__ * a[mvvlsb_i__] - mvvlsb_sum;
			mvvlsb_infa = 1;
		    }
		}
		if (infi[mvvlsb_i__] != 1) {
		    if (mvvlsb_infb == 1) {
			/* Computing MIN */
			d__1 = mvvlsb_bi, d__2 = *r__ * b[mvvlsb_i__] - mvvlsb_sum;
			mvvlsb_bi = std::min(d__1,d__2);
		    } else {
			mvvlsb_bi = *r__ * b[mvvlsb_i__] - mvvlsb_sum;
			mvvlsb_infb = 1;
		    }
		}
		++mvvlsb_ij;
		if (mvvlsb_i__ == *n || cov[mvvlsb_ij + *nd + 2] > 0.) {
		    i__2 = mvvlsb_infa + mvvlsb_infa + mvvlsb_infb - 1;
		    mvlims_(&mvvlsb_ai, &mvvlsb_bi, &i__2, di, ei);
		    if (*di >= *ei) {
			*value = 0.;
			return 0;
		    } else {
			*value *= *ei - *di;
			++(*nd);
			if (mvvlsb_i__ < *n) {
			    d__1 = *di + w[*nd] * (*ei - *di);
			    y[*nd] = mvphnv_(&d__1);
			}
			mvvlsb_infa = 0;
			mvvlsb_infb = 0;
		    }
		}
	    }
	    return 0;
	} /* mvvlsb_ */


	int mvsort_(int *n, double *lower, double *upper,
			 double *delta, double *correl, int *infin, double *y,
			 int *pivot, int *nd, double *a, double *b, 
			 double *dl, double *cov, int *infi, int *inform__)
	{
	    /* System generated locals */
	    int i__1, i__2, i__3, i__4;
	    double d__1;

	    /* Local variables */
	    double d__, e;
	    int i__, j, k, l, m;
	    double aj, bj;
	    int ii, ij, il, jl;
	    double sum, amin, bmin;
	    int jmin;
	    double epsi, demin, sumsq, cvdiag;

	    /*     Subroutine to sort integration limits and determine Cholesky factor. */

	    /* Parameter adjustments */
	    --infi;
	    --cov;
	    --dl;
	    --b;
	    --a;
	    --y;
	    --infin;
	    --correl;
	    --delta;
	    --upper;
	    --lower;

	    /* Function Body */
	    *inform__ = 0;
	    ij = 0;
	    ii = 0;
	    *nd = *n;
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		a[i__] = 0.;
		b[i__] = 0.;
		dl[i__] = 0.;
		infi[i__] = infin[i__];
		if (infi[i__] < 0) {
		    --(*nd);
		} else {
		    if (infi[i__] != 0) {
			a[i__] = lower[i__];
		    }
		    if (infi[i__] != 1) {
			b[i__] = upper[i__];
		    }
		    dl[i__] = delta[i__];
		}
		i__2 = i__ - 1;
		for (j = 1; j <= i__2; ++j) {
		    ++ij;
		    ++ii;
		    cov[ij] = correl[ii];
		}
		++ij;
		cov[ij] = 1.;
	    }

	    /*     First move any doubly infinite limits to innermost positions. */

	    if (*nd > 0) {
		i__1 = *nd + 1;
		for (i__ = *n; i__ >= i__1; --i__) {
		    if (infi[i__] >= 0) {
			i__2 = i__ - 1;
			for (j = 1; j <= i__2; ++j) {
			    if (infi[j] < 0) {
				mvswap_(&j, &i__, &a[1], &b[1], &dl[1], &infi[1], n, &
					cov[1]);
				goto L10;
			    }
			}
		    }
		L10:
		    ;
		}

		/*     Sort remaining limits and determine Cholesky factor. */

		ii = 0;
		jl = *nd;
		i__1 = *nd;
		for (i__ = 1; i__ <= i__1; ++i__) {

		    /*        Determine the integration limits for variable with minimum */
		    /*        expected probability and interchange that variable with Ith. */

		    demin = 1.;
		    jmin = i__;
		    cvdiag = 0.;
		    ij = ii;
		    epsi = i__ * 1e-10;
		    if (! (*pivot)) {
			jl = i__;
		    }
		    i__2 = jl;
		    for (j = i__; j <= i__2; ++j) {
			if (cov[ij + j] > epsi) {
			    sumsq = std::sqrt(cov[ij + j]);
			    sum = dl[j];
			    i__3 = i__ - 1;
			    for (k = 1; k <= i__3; ++k) {
				sum += cov[ij + k] * y[k];
			    }
			    aj = (a[j] - sum) / sumsq;
			    bj = (b[j] - sum) / sumsq;
			    mvlims_(&aj, &bj, &infi[j], &d__, &e);
			    if (demin >= e - d__) {
				jmin = j;
				amin = aj;
				bmin = bj;
				demin = e - d__;
				cvdiag = sumsq;
			    }
			}
			ij += j;
		    }
		    if (jmin > i__) {
			mvswap_(&i__, &jmin, &a[1], &b[1], &dl[1], &infi[1], n, &cov[
				    1]);
		    }
		    if (cov[ii + i__] < -epsi) {
			*inform__ = 3;
		    }
		    cov[ii + i__] = cvdiag;

		    /*        Compute Ith column of Cholesky factor. */
		    /*        Compute expected value for Ith integration variable and */
		    /*         scale Ith covariance matrix row and limits. */

		    if (cvdiag > 0.) {
			il = ii + i__;
			i__2 = *nd;
			for (l = i__ + 1; l <= i__2; ++l) {
			    cov[il + i__] /= cvdiag;
			    ij = ii + i__;
			    i__3 = l;
			    for (j = i__ + 1; j <= i__3; ++j) {
				cov[il + j] -= cov[il + i__] * cov[ij + i__];
				ij += j;
			    }
			    il += l;
			}

			/*              Expected Y = -( density(b) - density(a) )/( b - a ) */

			if (demin > epsi) {
			    y[i__] = 0.;
			    if (infi[i__] != 0) {
				y[i__] = mvtdns_(&c__0, &amin);
			    }
			    if (infi[i__] != 1) {
				y[i__] -= mvtdns_(&c__0, &bmin);
			    }
			    y[i__] /= demin;
			} else {
			    if (infi[i__] == 0) {
				y[i__] = bmin;
			    }
			    if (infi[i__] == 1) {
				y[i__] = amin;
			    }
			    if (infi[i__] == 2) {
				y[i__] = (amin + bmin) / 2;
			    }
			}
			i__2 = i__;
			for (j = 1; j <= i__2; ++j) {
			    ++ii;
			    cov[ii] /= cvdiag;
			}
			a[i__] /= cvdiag;
			b[i__] /= cvdiag;
			dl[i__] /= cvdiag;
		    } else {
			il = ii + i__;
			i__2 = *nd;
			for (l = i__ + 1; l <= i__2; ++l) {
			    cov[il + i__] = 0.;
			    il += l;
			}

			/*        If the covariance matrix diagonal entry is zero, */
			/*         permute limits and rows, if necessary. */


			for (j = i__ - 1; j >= 1; --j) {
			    if ((d__1 = cov[ii + j], std::abs(d__1)) > epsi) {
				a[i__] /= cov[ii + j];
				b[i__] /= cov[ii + j];
				dl[i__] /= cov[ii + j];
				if (cov[ii + j] < 0.) {
				    mvsswp_(&a[i__], &b[i__]);
				    if (infi[i__] != 2) {
					infi[i__] = 1 - infi[i__];
				    }
				}
				i__2 = j;
				for (l = 1; l <= i__2; ++l) {
				    cov[ii + l] /= cov[ii + j];
				}
				i__2 = i__ - 1;
				for (l = j + 1; l <= i__2; ++l) {
				    if (cov[(l - 1) * l / 2 + j + 1] > 0.) {
					ij = ii;
					i__3 = l;
					for (k = i__ - 1; k >= i__3; --k) {
					    i__4 = k;
					    for (m = 1; m <= i__4; ++m) {
						mvsswp_(&cov[ij - k + m], &cov[ij + m]
						    );
					    }
					    mvsswp_(&a[k], &a[k + 1]);
					    mvsswp_(&b[k], &b[k + 1]);
					    mvsswp_(&dl[k], &dl[k + 1]);
					    m = infi[k];
					    infi[k] = infi[k + 1];
					    infi[k + 1] = m;
					    ij -= k;
					}
					goto L20;
				    }
				}
				goto L20;
			    }
			    cov[ii + j] = 0.;
			}
		    L20:
			ii += i__;
			y[i__] = 0.;
		    }
		}
	    }
	    return 0;
	} /* mvsort_ */


	double mvtdns_(int *nu, double *x)
	{
	    /* System generated locals */
	    int i__1;
	    double ret_val, d__1;

	    /* Local variables */
	    int i__;
	    double prod;

	    ret_val = 0.;
	    if (*nu > 0) {
		prod = 1 / std::sqrt((double) (*nu));
		for (i__ = *nu - 2; i__ >= 1; i__ += -2) {
		    prod = prod * (i__ + 1) / i__;
		}
		if (*nu % 2 == 0) {
		    prod /= 2;
		} else {
		    prod /= 3.141592653589793;
		}
		d__1 = std::sqrt(*x * *x / *nu + 1);
		i__1 = *nu + 1;
		ret_val = prod / pow_di(&d__1, &i__1);
	    } else {
		if (std::abs(*x) < 10.) {
		    ret_val = std::exp(-(*x) * *x / 2) / 2.506628274631001;
		}
	    }
	    return ret_val;
	} /* mvtdns_ */


	/* Subroutine */ int mvlims_(double *a, double *b, int *infin, 
					  double *lower, double *upper)
	{

	    *lower = 0.;
	    *upper = 1.;
	    if (*infin >= 0) {
		if (*infin != 0) {
		    *lower = mvphi_(a);
		}
		if (*infin != 1) {
		    *upper = mvphi_(b);
		}
	    }
	    *upper = std::max(*upper,*lower);
	    return 0;
	} /* mvlims_ */


	/* Subroutine */ int mvsswp_(double *x, double *y)
	{
	    double t;

	    t = *x;
	    *x = *y;
	    *y = t;
	    return 0;
	} /* mvsswp_ */


	/* Subroutine */ int mvswap_(int *p, int *q, double *a, 
					  double *b, double *d__, int *infin, int *n, 
					  double *c__)
	{
	    /* System generated locals */
	    int i__1;

	    /* Local variables */
	    int i__, j, ii, jj;

	    /*     Swaps rows and columns P and Q in situ, with P <= Q. */

	    /* Parameter adjustments */
	    --c__;
	    --infin;
	    --d__;
	    --b;
	    --a;

	    /* Function Body */
	    mvsswp_(&a[*p], &a[*q]);
	    mvsswp_(&b[*p], &b[*q]);
	    mvsswp_(&d__[*p], &d__[*q]);
	    j = infin[*p];
	    infin[*p] = infin[*q];
	    infin[*q] = j;
	    jj = *p * (*p - 1) / 2;
	    ii = *q * (*q - 1) / 2;
	    mvsswp_(&c__[jj + *p], &c__[ii + *q]);
	    i__1 = *p - 1;
	    for (j = 1; j <= i__1; ++j) {
		mvsswp_(&c__[jj + j], &c__[ii + j]);
	    }
	    jj += *p;
	    i__1 = *q - 1;
	    for (i__ = *p + 1; i__ <= i__1; ++i__) {
		mvsswp_(&c__[jj + *p], &c__[ii + i__]);
		jj += i__;
	    }
	    ii += *q;
	    i__1 = *n;
	    for (i__ = *q + 1; i__ <= i__1; ++i__) {
		mvsswp_(&c__[ii + *p], &c__[ii + *q]);
		ii += i__;
	    }
	    return 0;
	} /* mvswap_ */


	double mvphi_(double *z__)
	{
	    /* Initialized data */

	    static const double a[44] = { .610143081923200417926465815756,
				    -.434841272712577471828182820888,.176351193643605501125840298123,
				    -.060710795609249414860051215825,.017712068995694114486147141191,
				    -.004321119385567293818599864968,8.54216676887098678819832055e-4,
				    -1.2715509060916274262889394e-4,1.1248167243671189468847072e-5,
				    3.13063885421820972630152e-7,-2.70988068537762022009086e-7,
				    3.0737622701407688440959e-8,2.515620384817622937314e-9,
				    -1.02892992132031912759e-9,2.9944052119949939363e-11,
				    2.605178968726693629e-11,-2.634839924171969386e-12,
				    -6.43404509890636443e-13,1.12457401801663447e-13,
				    1.7281533389986098e-14,-4.264101694942375e-15,
				    -5.45371977880191e-16,1.58697607761671e-16,2.0899837844334e-17,
				    -5.900526869409e-18,-9.41893387554e-19,2.1497735647e-19,
				    4.6660985008e-20,-7.243011862e-21,-2.387966824e-21,1.91177535e-22,
				    1.20482568e-22,-6.72377e-25,-5.747997e-24,-4.28493e-25,
				    2.44856e-25,4.3793e-26,-8.151e-27,-3.089e-27,9.3e-29,1.74e-28,
				    1.6e-29,-8e-30,-2e-30 };

	    /* System generated locals */
	    double ret_val;

	    /* Local variables */
	    double b;
	    int i__;
	    double p, t, bm, bp, xa;


	    /*     Normal distribution probabilities accurate to 1d-15. */
	    /*     Reference: J.L. Schonfelder, Math Comp 32(1978), pp 1232-1240. */


	    xa = std::abs(*z__) / 1.414213562373095048801688724209;
	    if (xa > 100.) {
		p = 0.;
	    } else {
		t = (xa * 8 - 30) / (xa * 4 + 15);
		bm = 0.;
		b = 0.;
		for (i__ = 24; i__ >= 0; --i__) {
		    bp = b;
		    b = bm;
		    bm = t * b - bp + a[i__];
		}
		p = std::exp(-xa * xa) * (bm - bp) / 4;
	    }
	    if (*z__ > 0.) {
		p = 1 - p;
	    }
	    ret_val = p;
	    return ret_val;
	} /* mvphi_ */


	double mvphnv_(double *p)
	{
	    /* System generated locals */
	    double ret_val, d__1, d__2;

	    /* Local variables */
	    double q, r__;


	    /* 	ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3 */

	    /* 	Produces the normal deviate Z corresponding to a given lower */
	    /* 	tail area of P. */

	    /* 	The hash sums below are the sums of the mantissas of the */
	    /* 	coefficients.   They are included for use in checking */
	    /* 	transcription. */


	    /*     Coefficients for P close to 0.5 */

	    /*     HASH SUM AB    55.88319 28806 14901 4439 */

	    /*     Coefficients for P not close to 0, 0.5 or 1. */

	    /*     HASH SUM CD    49.33206 50330 16102 89036 */

	    /* 	Coefficients for P near 0 or 1. */

	    /*     HASH SUM EF    47.52583 31754 92896 71629 */

	    q = (*p * 2 - 1) / 2;
	    if (std::abs(q) <= .425) {
		r__ = .180625 - q * q;
		ret_val = q * (((((((r__ * 2509.0809287301226727 + 
				     33430.575583588128105) * r__ + 67265.770927008700853) * r__ + 
				   45921.953931549871457) * r__ + 13731.693765509461125) * r__ + 
				 1971.5909503065514427) * r__ + 133.14166789178437745) * r__ + 
			       3.387132872796366608) / (((((((r__ * 5226.495278852854561 + 
							      28729.085735721942674) * r__ + 39307.89580009271061) * r__ + 
							    21213.794301586595867) * r__ + 5394.1960214247511077) * r__ + 
							  687.1870074920579083) * r__ + 42.313330701600911252) * r__ + 
							1);
	    } else {
		/* Computing MIN */
		d__1 = *p, d__2 = 1 - *p;
		r__ = std::min(d__1,d__2);
		if (r__ > 0.) {
		    r__ = std::sqrt(-std::log(r__));
		    if (r__ <= 5.) {
			r__ += -1.6;
			ret_val = (((((((r__ * 7.7454501427834140764e-4 + 
					 .0227238449892691845833) * r__ + 
					.24178072517745061177) * r__ + 1.27045825245236838258)
				      * r__ + 3.64784832476320460504) * r__ + 
				     5.7694972214606914055) * r__ + 4.6303378461565452959) 
				   * r__ + 1.42343711074968357734) / (((((((r__ * 
									    1.05075007164441684324e-9 + 5.475938084995344946e-4) *
									   r__ + .0151986665636164571966) * r__ + 
									  .14810397642748007459) * r__ + .68976733498510000455) 
									* r__ + 1.6763848301838038494) * r__ + 
								       2.05319162663775882187) * r__ + 1);
		    } else {
			r__ += -5.;
			ret_val = (((((((r__ * 2.01033439929228813265e-7 + 
					 2.71155556874348757815e-5) * r__ + 
					.0012426609473880784386) * r__ + 
				       .026532189526576123093) * r__ + .29656057182850489123)
				     * r__ + 1.7848265399172913358) * r__ + 
				    5.4637849111641143699) * r__ + 6.6579046435011037772) 
			    / (((((((r__ * 2.04426310338993978564e-15 + 
				     1.4215117583164458887e-7) * r__ + 
				    1.8463183175100546818e-5) * r__ + 
				   7.868691311456132591e-4) * r__ + 
				  .0148753612908506148525) * r__ + 
				 .13692988092273580531) * r__ + .59983220655588793769) 
			       * r__ + 1);
		    }
		} else {
		    ret_val = 9.;
		}
		if (q < 0.) {
		    ret_val = -ret_val;
		}
	    }
	    return ret_val;
	} /* mvphnv_ */

	double mvbvn_(double *lower, double *upper, int *infin, 
		      double *correl)
	{
	    /* System generated locals */
	    double ret_val, d__1, d__2, d__3, d__4;

	    /* Local variables */

	    /*     A function for computing bivariate normal probabilities. */

	    /*  Parameters */

	    /*     LOWER  REAL, array of lower integration limits. */
	    /*     UPPER  REAL, array of upper integration limits. */
	    /*     INFIN  INT, array of integration limits flags: */
	    /*            if INFIN(I) = 0, Ith limits are (-infinity, UPPER(I)]; */
	    /*            if INFIN(I) = 1, Ith limits are [LOWER(I), infinity); */
	    /*            if INFIN(I) = 2, Ith limits are [LOWER(I), UPPER(I)]. */
	    /*     CORREL REAL, correlation coefficient. */

	    /* Parameter adjustments */
	    --infin;
	    --upper;
	    --lower;

	    /* Function Body */
	    if (infin[1] == 2 && infin[2] == 2) {
		ret_val = mvbvu_(&lower[1], &lower[2], correl) - mvbvu_(&upper[1], &
									lower[2], correl) - mvbvu_(&lower[1], &upper[2], correl) + 
		    mvbvu_(&upper[1], &upper[2], correl);
	    } else if (infin[1] == 2 && infin[2] == 1) {
		ret_val = mvbvu_(&lower[1], &lower[2], correl) - mvbvu_(&upper[1], &
									lower[2], correl);
	    } else if (infin[1] == 1 && infin[2] == 2) {
		ret_val = mvbvu_(&lower[1], &lower[2], correl) - mvbvu_(&lower[1], &
									upper[2], correl);
	    } else if (infin[1] == 2 && infin[2] == 0) {
		d__1 = -upper[1];
		d__2 = -upper[2];
		d__3 = -lower[1];
		d__4 = -upper[2];
		ret_val = mvbvu_(&d__1, &d__2, correl) - mvbvu_(&d__3, &d__4, correl);
	    } else if (infin[1] == 0 && infin[2] == 2) {
		d__1 = -upper[1];
		d__2 = -upper[2];
		d__3 = -upper[1];
		d__4 = -lower[2];
		ret_val = mvbvu_(&d__1, &d__2, correl) - mvbvu_(&d__3, &d__4, correl);
	    } else if (infin[1] == 1 && infin[2] == 0) {
		d__1 = -upper[2];
		d__2 = -(*correl);
		ret_val = mvbvu_(&lower[1], &d__1, &d__2);
	    } else if (infin[1] == 0 && infin[2] == 1) {
		d__1 = -upper[1];
		d__2 = -(*correl);
		ret_val = mvbvu_(&d__1, &lower[2], &d__2);
	    } else if (infin[1] == 1 && infin[2] == 1) {
		ret_val = mvbvu_(&lower[1], &lower[2], correl);
	    } else if (infin[1] == 0 && infin[2] == 0) {
		d__1 = -upper[1];
		d__2 = -upper[2];
		ret_val = mvbvu_(&d__1, &d__2, correl);
	    } else {
		ret_val = 1.;
	    }
	    return ret_val;
	} /* mvbvn_ */

	double mvbvu_(double *sh, double *sk, double *r__)
	{
	    /* Initialized data */

	    static const struct {
		double e_1[3];
		double fill_2[7];
		double e_3[6];
		double fill_4[4];
		double e_5[10];
	    } equiv_83 = { .1713244923791705, .3607615730481384, 
			   .4679139345726904, {0}, .04717533638651177, .1069393259953183,
			   .1600783285433464, .2031674267230659, .2334925365383547, 
			   .2491470458134029, {0}, .01761400713915212, 
			   .04060142980038694, .06267204833410906, .08327674157670475, 
			   .1019301198172404, .1181945319615184, .1316886384491766, 
			   .1420961093183821, .1491729864726037, .1527533871307259 };

#define w ((double *)&equiv_83)

	    static const struct {
		double e_1[3];
		double fill_2[7];
		double e_3[6];
		double fill_4[4];
		double e_5[10];
	    } equiv_84 = { -.9324695142031522, -.6612093864662647, 
			   -.238619186083197, {0}, -.9815606342467191, -.904117256370475,
			   -.769902674194305, -.5873179542866171, -.3678314989981802, 
			   -.1252334085114692, {0}, -.9931285991850949, 
			   -.9639719272779138, -.9122344282513259, -.8391169718222188, 
			   -.7463319064601508, -.636053680726515, -.5108670019508271, 
			   -.3737060887154196, -.2277858511416451, -.07652652113349733 };

#define x ((double *)&equiv_84)


	    /* System generated locals */
	    int i__1;
	    double ret_val, d__1, d__2;

	    /* Local variables */
	    double a, b, c__, d__, h__;
	    int i__;
	    double k;
	    int lg;
	    double as;
	    int ng;
	    double bs, hk, hs, sn, rs, xs, bvn, asr;

	    /*     A function for computing bivariate normal probabilities; */
	    /*       developed using */
	    /*         Drezner, Z. and Wesolowsky, G. O. (1989), */
	    /*         On the Computation of the Bivariate Normal Integral, */
	    /*         J. Stat. Comput. Simul.. 35 pp. 101-107. */
	    /*       with extensive modications for double precisions by */
	    /*         Alan Genz and Yihong Ge */
	    /*         Department of Mathematics */
	    /*         Washington State University */
	    /*         Pullman, WA 99164-3113 */
	    /*         Email : alangenz@wsu.edu */

	    /* BVN - calculate the probability that X is larger than SH and Y is */
	    /*       larger than SK. */

	    /* Parameters */

	    /*   SH  REAL, integration limit */
	    /*   SK  REAL, integration limit */
	    /*   R   REAL, correlation coefficient */
	    /*   LG  INT, number of Gauss Rule Points and Weights */

	    /*     Gauss Legendre Points and Weights, N =  6 */
	    /*     Gauss Legendre Points and Weights, N = 12 */
	    /*     Gauss Legendre Points and Weights, N = 20 */
	    if (std::abs(*r__) < .3f) {
		ng = 1;
		lg = 3;
	    } else if (std::abs(*r__) < .75f) {
		ng = 2;
		lg = 6;
	    } else {
		ng = 3;
		lg = 10;
	    }
	    h__ = *sh;
	    k = *sk;
	    hk = h__ * k;
	    bvn = 0.;
	    if (std::abs(*r__) < .925f) {
		hs = (h__ * h__ + k * k) / 2;
		asr = std::asin(*r__);
		i__1 = lg;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    sn = std::sin(asr * (x[i__ + ng * 10 - 11] + 1) / 2);
		    bvn += w[i__ + ng * 10 - 11] * std::exp((sn * hk - hs) / (1 - sn * sn))
			;
		    sn = std::sin(asr * (-x[i__ + ng * 10 - 11] + 1) / 2);
		    bvn += w[i__ + ng * 10 - 11] * std::exp((sn * hk - hs) / (1 - sn * sn))
			;
		}
		d__1 = -h__;
		d__2 = -k;
		bvn = bvn * asr / 12.566370614359172 + mvphi_(&d__1) * mvphi_(&d__2);
	    } else {
		if (*r__ < 0.) {
		    k = -k;
		    hk = -hk;
		}
		if (std::abs(*r__) < 1.) {
		    as = (1 - *r__) * (*r__ + 1);
		    a = std::sqrt(as);
		    /* Computing 2nd power */
		    d__1 = h__ - k;
		    bs = d__1 * d__1;
		    c__ = (4 - hk) / 8;
		    d__ = (12 - hk) / 16;
		    bvn = a * std::exp(-(bs / as + hk) / 2) * (1 - c__ * (bs - as) * (1 - 
										      d__ * bs / 5) / 3 + c__ * d__ * as * as / 5);
		    if (hk > -160.) {
			b = std::sqrt(bs);
			d__1 = -b / a;
			bvn -= std::exp(-hk / 2) * std::sqrt(6.283185307179586) * mvphi_(&d__1) 
			    * b * (1 - c__ * bs * (1 - d__ * bs / 5) / 3);
		    }
		    a /= 2;
		    i__1 = lg;
		    for (i__ = 1; i__ <= i__1; ++i__) {
			/* Computing 2nd power */
			d__1 = a * (x[i__ + ng * 10 - 11] + 1);
			xs = d__1 * d__1;
			rs = std::sqrt(1 - xs);
			bvn += a * w[i__ + ng * 10 - 11] * (std::exp(-bs / (xs * 2) - hk / 
								     (rs + 1)) / rs - std::exp(-(bs / xs + hk) / 2) * (c__ * xs 
														       * (d__ * xs + 1) + 1));
			/* Computing 2nd power */
			d__1 = -x[i__ + ng * 10 - 11] + 1;
			xs = as * (d__1 * d__1) / 4;
			rs = std::sqrt(1 - xs);
			/* Computing 2nd power */
			d__1 = rs + 1;
			bvn += a * w[i__ + ng * 10 - 11] * std::exp(-(bs / xs + hk) / 2) * 
			    (std::exp(-hk * xs / (d__1 * d__1 * 2)) / rs - (c__ * xs * 
									    (d__ * xs + 1) + 1));
		    }
		    bvn = -bvn / 6.283185307179586;
		}
		if (*r__ > 0.) {
		    d__1 = -std::max(h__,k);
		    bvn += mvphi_(&d__1);
		} else {
		    bvn = -bvn;
		    if (k > h__) {
			if (h__ < 0.) {
			    bvn = bvn + mvphi_(&k) - mvphi_(&h__);
			} else {
			    d__1 = -h__;
			    d__2 = -k;
			    bvn = bvn + mvphi_(&d__1) - mvphi_(&d__2);
			}
		    }
		}
	    }
	    ret_val = bvn;
	    return ret_val;
	} /* mvbvu_ */

#undef x
#undef w



	double mvstdt_(int *nu, double *t)
	{
	    /* System generated locals */
	    double ret_val;

	    /* Local variables */
	    int j;
	    double rn, ts, tt, csthe, snthe;
	    double polyn;


	    /*     Student t Distribution Function */

	    /*                       T */
	    /*         TSTDNT = C   I  ( 1 + y*y/NU )**( -(NU+1)/2 ) dy */
	    /*                   NU -INF */

	    if (*nu < 1) {
		ret_val = mvphi_(t);
	    } else if (*nu == 1) {
		ret_val = (std::atan(*t) * 2 / 3.141592653589793 + 1) / 2;
	    } else if (*nu == 2) {
		ret_val = (*t / std::sqrt(*t * *t + 2) + 1) / 2;
	    } else {
		tt = *t * *t;
		csthe = *nu / (*nu + tt);
		polyn = 1.;
		for (j = *nu - 2; j >= 2; j += -2) {
		    polyn = (j - 1) * csthe * polyn / j + 1;
		}
		if (*nu % 2 == 1) {
		    rn = (double) (*nu);
		    ts = *t / std::sqrt(rn);
		    ret_val = ((atan(ts) + ts * csthe * polyn) * 2 / 
			       3.141592653589793 + 1) / 2;
		} else {
		    snthe = *t / std::sqrt(*nu + tt);
		    ret_val = (snthe * polyn + 1) / 2;
		}
		if (ret_val < 0.) {
		    ret_val = 0.;
		}
	    }
	    return ret_val;
	} /* mvstdt_ */


	double mvbvt_(int *nu, double *lower, double *upper, int *
		      infin, double *correl)
	{
	    /* System generated locals */
	    double ret_val, d__1, d__2, d__3, d__4;

	    /* Local variables */

	    /*     A function for computing bivariate normal and t probabilities. */

	    /*  Parameters */

	    /*     NU     INT degrees of freedom parameter; NU < 1 gives normal case. */
	    /*     LOWER  REAL, array of lower integration limits. */
	    /*     UPPER  REAL, array of upper integration limits. */
	    /*     INFIN  INT, array of integration limits flags: */
	    /*            if INFIN(I) = 0, Ith limits are (-infinity, UPPER(I)]; */
	    /*            if INFIN(I) = 1, Ith limits are [LOWER(I), infinity); */
	    /*            if INFIN(I) = 2, Ith limits are [LOWER(I), UPPER(I)]. */
	    /*     CORREL REAL, correlation coefficient. */

	    /* Parameter adjustments */
	    --infin;
	    --upper;
	    --lower;

	    /* Function Body */
	    if (*nu < 1) {
		ret_val = mvbvn_(&lower[1], &upper[1], &infin[1], correl);
	    } else {
		if (infin[1] == 2 && infin[2] == 2) {
		    ret_val = mvbvtl_(nu, &upper[1], &upper[2], correl) - mvbvtl_(nu, 
										  &upper[1], &lower[2], correl) - mvbvtl_(nu, &lower[1], &
															  upper[2], correl) + mvbvtl_(nu, &lower[1], &lower[2], 
																		      correl);
		} else if (infin[1] == 2 && infin[2] == 1) {
		    d__1 = -lower[1];
		    d__2 = -lower[2];
		    d__3 = -upper[1];
		    d__4 = -lower[2];
		    ret_val = mvbvtl_(nu, &d__1, &d__2, correl) - mvbvtl_(nu, &d__3, &
									  d__4, correl);
		} else if (infin[1] == 1 && infin[2] == 2) {
		    d__1 = -lower[1];
		    d__2 = -lower[2];
		    d__3 = -lower[1];
		    d__4 = -upper[2];
		    ret_val = mvbvtl_(nu, &d__1, &d__2, correl) - mvbvtl_(nu, &d__3, &
									  d__4, correl);
		} else if (infin[1] == 2 && infin[2] == 0) {
		    ret_val = mvbvtl_(nu, &upper[1], &upper[2], correl) - mvbvtl_(nu, 
										  &lower[1], &upper[2], correl);
		} else if (infin[1] == 0 && infin[2] == 2) {
		    ret_val = mvbvtl_(nu, &upper[1], &upper[2], correl) - mvbvtl_(nu, 
										  &upper[1], &lower[2], correl);
		} else if (infin[1] == 1 && infin[2] == 0) {
		    d__1 = -lower[1];
		    d__2 = -(*correl);
		    ret_val = mvbvtl_(nu, &d__1, &upper[2], &d__2);
		} else if (infin[1] == 0 && infin[2] == 1) {
		    d__1 = -lower[2];
		    d__2 = -(*correl);
		    ret_val = mvbvtl_(nu, &upper[1], &d__1, &d__2);
		} else if (infin[1] == 1 && infin[2] == 1) {
		    d__1 = -lower[1];
		    d__2 = -lower[2];
		    ret_val = mvbvtl_(nu, &d__1, &d__2, correl);
		} else if (infin[1] == 0 && infin[2] == 0) {
		    ret_val = mvbvtl_(nu, &upper[1], &upper[2], correl);
		} else {
		    ret_val = 1.;
		}
	    }
	    return ret_val;
	} /* mvbvt_ */


	double mvbvtc_(int *nu, double *l, double *u, int *infin, 
		       double *rho)
	{
	    /* System generated locals */
	    double ret_val;

	    /* Local variables */
	    double b;
	    int i__;
	    double lw[2], up[2];
	    int inf[2];

	    /*     A function for computing complementary bivariate normal and t */
	    /*       probabilities. */

	    /*  Parameters */

	    /*     NU     INT degrees of freedom parameter. */
	    /*     L      REAL, array of lower integration limits. */
	    /*     U      REAL, array of upper integration limits. */
	    /*     INFIN  INT, array of integration limits flags: */
	    /*            if INFIN(1) INFIN(2),        then MVBVTC computes */
	    /*                 0         0              P( X>U(1), Y>U(2) ) */
	    /*                 1         0              P( X<L(1), Y>U(2) ) */
	    /*                 0         1              P( X>U(1), Y<L(2) ) */
	    /*                 1         1              P( X<L(1), Y<L(2) ) */
	    /*                 2         0      P( X>U(1), Y>U(2) ) + P( X<L(1), Y>U(2) ) */
	    /*                 2         1      P( X>U(1), Y<L(2) ) + P( X<L(1), Y<L(2) ) */
	    /*                 0         2      P( X>U(1), Y>U(2) ) + P( X>U(1), Y<L(2) ) */
	    /*                 1         2      P( X<L(1), Y>U(2) ) + P( X<L(1), Y<L(2) ) */
	    /*                 2         2      P( X>U(1), Y<L(2) ) + P( X<L(1), Y<L(2) ) */
	    /*                               +  P( X>U(1), Y>U(2) ) + P( X<L(1), Y>U(2) ) */

	    /*     RHO    REAL, correlation coefficient. */


	    /* Parameter adjustments */
	    --infin;
	    --u;
	    --l;

	    /* Function Body */
	    for (i__ = 1; i__ <= 2; ++i__) {
		if (infin[i__] % 2 == 0) {
		    inf[i__ - 1] = 1;
		    lw[i__ - 1] = u[i__];
		} else {
		    inf[i__ - 1] = 0;
		    up[i__ - 1] = l[i__];
		}
	    }
	    b = mvbvt_(nu, lw, up, inf, rho);
	    for (i__ = 1; i__ <= 2; ++i__) {
		if (infin[i__] == 2) {
		    inf[i__ - 1] = 0;
		    up[i__ - 1] = l[i__];
		    b += mvbvt_(nu, lw, up, inf, rho);
		}
	    }
	    if (infin[1] == 2 && infin[2] == 2) {
		inf[0] = 1;
		lw[0] = u[1];
		b += mvbvt_(nu, lw, up, inf, rho);
	    }
	    ret_val = b;
	    return ret_val;
	} /* mvbvtc_ */


	double mvbvtl_(int *nu, double *dh, double *dk, double *r__)
	{
	    /* System generated locals */
	    int i__1;
	    double ret_val, d__1, d__2, d__3;

	    /* Local variables */
	    int j, hs, ks;
	    double hkn, hpk, hrk, krh, bvt, ors, snu, gmph, gmpk, hkrn, 
		qhrk, xnkh, xnhk, btnckh, btnchk, btpdkh, btpdhk;


	    /*     a function for computing bivariate t probabilities. */

	    /*       Alan Genz */
	    /*       Department of Mathematics */
	    /*       Washington State University */
	    /*       Pullman, Wa 99164-3113 */
	    /*       Email : alangenz@wsu.edu */

	    /*    this function is based on the method described by */
	    /*        Dunnett, C.W. and M. Sobel, (1954), */
	    /*        A bivariate generalization of Student's t-distribution */
	    /*        with tables for certain special cases, */
	    /*        Biometrika 41, pp. 153-169. */

	    /* mvbvtl - calculate the probability that x < dh and y < dk. */

	    /* parameters */

	    /*   nu number of degrees of freedom */
	    /*   dh 1st lower integration limit */
	    /*   dk 2nd lower integration limit */
	    /*   r   correlation coefficient */

	    snu = std::sqrt((double) (*nu));
	    ors = 1 - *r__ * *r__;
	    hrk = *dh - *r__ * *dk;
	    krh = *dk - *r__ * *dh;
	    if (std::abs(hrk) + ors > 0.) {
		/* Computing 2nd power */
		d__1 = hrk;
		/* Computing 2nd power */
		d__2 = hrk;
		/* Computing 2nd power */
		d__3 = *dk;
		xnhk = d__1 * d__1 / (d__2 * d__2 + ors * (*nu + d__3 * d__3));
		/* Computing 2nd power */
		d__1 = krh;
		/* Computing 2nd power */
		d__2 = krh;
		/* Computing 2nd power */
		d__3 = *dh;
		xnkh = d__1 * d__1 / (d__2 * d__2 + ors * (*nu + d__3 * d__3));
	    } else {
		xnhk = 0.;
		xnkh = 0.;
	    }
	    d__1 = *dh - *r__ * *dk;
	    hs = (int) d_sign(&c_b24, &d__1);
	    d__1 = *dk - *r__ * *dh;
	    ks = (int) d_sign(&c_b24, &d__1);
	    if (*nu % 2 == 0) {
		bvt = std::atan2(std::sqrt(ors), -(*r__)) / 6.2831853071795862;
		/* Computing 2nd power */
		d__1 = *dh;
		gmph = *dh / std::sqrt((*nu + d__1 * d__1) * 16);
		/* Computing 2nd power */
		d__1 = *dk;
		gmpk = *dk / std::sqrt((*nu + d__1 * d__1) * 16);
		btnckh = std::atan2(std::sqrt(xnkh), std::sqrt(1 - xnkh)) * 2 / 
		    3.14159265358979323844;
		btpdkh = std::sqrt(xnkh * (1 - xnkh)) * 2 / 3.14159265358979323844;
		btnchk = std::atan2(std::sqrt(xnhk), std::sqrt(1 - xnhk)) * 2 / 
		    3.14159265358979323844;
		btpdhk = std::sqrt(xnhk * (1 - xnhk)) * 2 / 3.14159265358979323844;
		i__1 = *nu / 2;
		for (j = 1; j <= i__1; ++j) {
		    bvt += gmph * (ks * btnckh + 1);
		    bvt += gmpk * (hs * btnchk + 1);
		    btnckh += btpdkh;
		    btpdkh = (j << 1) * btpdkh * (1 - xnkh) / ((j << 1) + 1);
		    btnchk += btpdhk;
		    btpdhk = (j << 1) * btpdhk * (1 - xnhk) / ((j << 1) + 1);
		    /* Computing 2nd power */
		    d__1 = *dh;
		    gmph = gmph * ((j << 1) - 1) / ((j << 1) * (d__1 * d__1 / *nu + 1)
			);
		    /* Computing 2nd power */
		    d__1 = *dk;
		    gmpk = gmpk * ((j << 1) - 1) / ((j << 1) * (d__1 * d__1 / *nu + 1)
			);
		}
	    } else {
		/* Computing 2nd power */
		d__1 = *dh;
		/* Computing 2nd power */
		d__2 = *dk;
		qhrk = std::sqrt(d__1 * d__1 + d__2 * d__2 - *r__ * 2 * *dh * *dk + *nu * 
				 ors);
		hkrn = *dh * *dk + *r__ * *nu;
		hkn = *dh * *dk - *nu;
		hpk = *dh + *dk;
		bvt = std::atan2(-snu * (hkn * qhrk + hpk * hkrn), hkn * hkrn - *nu * hpk *
				 qhrk) / 6.2831853071795862;
		if (bvt < -1e-15) {
		    bvt += 1;
		}
		/* Computing 2nd power */
		d__1 = *dh;
		gmph = *dh / (snu * 6.2831853071795862 * (d__1 * d__1 / *nu + 1));
		/* Computing 2nd power */
		d__1 = *dk;
		gmpk = *dk / (snu * 6.2831853071795862 * (d__1 * d__1 / *nu + 1));
		btnckh = std::sqrt(xnkh);
		btpdkh = btnckh;
		btnchk = std::sqrt(xnhk);
		btpdhk = btnchk;
		i__1 = (*nu - 1) / 2;
		for (j = 1; j <= i__1; ++j) {
		    bvt += gmph * (ks * btnckh + 1);
		    bvt += gmpk * (hs * btnchk + 1);
		    btpdkh = ((j << 1) - 1) * btpdkh * (1 - xnkh) / (j << 1);
		    btnckh += btpdkh;
		    btpdhk = ((j << 1) - 1) * btpdhk * (1 - xnhk) / (j << 1);
		    btnchk += btpdhk;
		    /* Computing 2nd power */
		    d__1 = *dh;
		    gmph = (j << 1) * gmph / (((j << 1) + 1) * (d__1 * d__1 / *nu + 1)
			);
		    /* Computing 2nd power */
		    d__1 = *dk;
		    gmpk = (j << 1) * gmpk / (((j << 1) + 1) * (d__1 * d__1 / *nu + 1)
			);
		}
	    }
	    ret_val = bvt;

	    /*     end mvbvtl */

	    return ret_val;
	} /* mvbvtl_ */


	double mvchnv_(int *n, double *p)
	{
	    /* System generated locals */
	    double ret_val;

	    /* Local variables */
	    double x;

	    /*                  MVCHNV */
	    /*     P =  1 - K  I     exp(-t*t/2) t**(N-1) dt, for N >= 1. */
	    /*               N  0 */

	    x = sqrtqchisqint_(n, p);
	    ret_val = x;
	    return ret_val;
	} /* mvchnv_ */
	template <class RNGType>
	/* Subroutine */ int mvkbrv_(int *ndim, int *minvls, int *maxvls, 
					  int *nf, double *abseps, double *releps, 
					  double *abserr, double *finest, int *inform__, RNGType *rng)
	{
	    /* Initialized data */

	    static const int p[28] = { 31,47,73,113,173,263,397,593,907,1361,2053,3079,
				      4621,6947,10427,15641,23473,35221,52837,79259,118891,178349,
				      267523,401287,601943,902933,1354471,2031713 };
	    static const int c__[2772] = { 12,13,27,35,64,111,163,
					  246,347,505,794,1189,1763,2872,4309,6610,9861,10327,19540,34566,
					  31929,40701,103650,165843,130365,333459,500884,858339,9,11,28,27,
					  66,42,154,189,402,220,325,888,1018,3233,3758,6977,3647,7582,19926,
					  9579,49367,69087,125480,90647,236711,375354,566009,918142,9,17,10,
					  27,28,54,83,242,322,601,960,259,1500,1534,4034,1686,4073,7124,
					  11582,12654,10982,77576,59978,59925,110235,102417,399251,501970,
					  13,10,11,36,28,118,43,102,418,644,528,1082,432,2941,1963,3819,
					  2535,8214,11113,26856,3527,64590,46875,189541,125699,383544,
					  652979,234813,12,15,11,22,44,20,82,250,215,612,247,725,1332,2910,
					  730,2314,3430,9600,24585,37873,27066,39397,77172,67647,56483,
					  292630,355008,460565,12,15,20,29,44,31,92,250,220,160,247,811,
					  2203,393,642,5647,9865,10271,8726,38806,13226,33179,83021,74795,
					  93735,41147,430235,31996,12,15,11,29,55,31,150,102,339,206,338,
					  636,126,1796,1502,3953,2830,10193,17218,29501,56010,10858,126904,
					  68365,234469,374614,328722,753018,12,15,11,20,67,72,59,250,339,
					  206,366,965,2240,919,2246,3614,9328,10800,419,17271,18911,38935,
					  14541,167485,60549,48032,670680,256150,12,15,28,45,10,17,76,280,
					  339,206,847,497,1719,446,3834,5115,4320,9086,4918,3663,40574,
					  43129,56299,143918,1291,435453,405585,199809,12,15,13,5,10,94,76,
					  118,337,422,753,497,1284,919,1511,423,5913,2365,4918,10763,20767,
					  35468,43636,74912,93937,281493,405585,993599,12,22,13,5,10,14,47,
					  196,218,134,753,1490,878,919,1102,423,10365,4409,4918,18955,20767,
					  35468,11655,167289,245291,358168,424646,245149,12,15,28,5,10,14,
					  11,118,315,518,236,1490,1983,1117,1102,5408,8272,13812,15701,1298,
					  9686,5279,52680,75517,196061,114121,670180,794183,3,15,13,21,10,
					  11,11,191,315,134,334,392,266,103,1522,7426,3706,5661,17710,26560,
					  47603,61518,88549,8148,258647,346892,670180,121349,3,6,13,21,10,
					  14,100,215,315,134,334,1291,266,103,1522,423,6186,9344,4037,17132,
					  47603,61518,29804,172106,162489,238990,641587,150619,3,6,13,21,38,
					  14,131,121,315,518,461,508,266,103,3427,423,7806,9344,4037,17132,
					  11736,27945,101894,126159,176631,317313,215580,376952,12,6,14,21,
					  38,14,116,121,167,652,711,508,266,103,3427,487,7806,10362,15808,
					  4753,11736,70975,113675,35867,204895,164158,59048,809123,7,15,14,
					  21,10,94,116,49,167,382,652,1291,747,103,3928,6227,7806,9344,
					  11401,4753,41601,70975,48040,35867,73353,35497,633320,809123,7,15,
					  14,21,10,10,116,49,167,206,381,1291,747,103,915,2660,8610,9344,
					  19398,8713,12888,86478,113675,35867,172319,70530,81010,804319,12,
					  9,14,21,10,10,116,49,167,158,381,508,127,103,915,6227,2563,8585,
					  25950,18624,32948,86478,34987,121694,28881,70530,20789,67352,12,
					  13,14,21,10,10,116,49,361,441,381,1291,127,2311,3818,1221,11558,
					  11114,25950,13082,30801,20514,48308,52171,136787,434839,389250,
					  969594,12,2,14,21,10,10,116,49,201,179,652,508,2074,3117,3818,
					  3811,11558,13080,4454,6791,44243,20514,97926,95354,122081,24754,
					  389250,434796,12,2,14,21,49,14,138,49,124,441,381,508,127,1101,
					  3818,197,9421,13080,24987,1122,53351,73178,5475,113969,122081,
					  24754,638764,969594,12,2,14,21,49,14,138,49,124,56,381,867,2074,
					  3117,3818,4367,1181,13080,11719,19363,53351,73178,49449,113969,
					  275993,24754,638764,804319,12,13,14,21,49,14,138,49,124,559,381,
					  867,1400,3117,4782,351,9421,6949,8697,34695,16016,43098,6850,
					  76304,64673,393656,389250,391368,12,11,14,21,49,14,138,49,124,559,
					  381,867,1383,1101,4782,1281,1181,3436,1452,18770,35086,43098,
					  62545,123709,211587,118711,389250,761041,12,11,14,21,49,14,138,49,
					  124,56,381,867,1383,1101,4782,1221,1181,3436,1452,18770,35086,
					  4701,62545,123709,211587,118711,398094,754049,12,10,14,21,49,14,
					  138,49,124,56,381,934,1383,1101,3818,351,1181,3436,1452,18770,
					  32581,59979,9440,144615,211587,148227,80846,466264,3,15,14,21,49,
					  14,138,49,124,56,381,867,1383,1101,4782,351,9421,13213,1452,18770,
					  2464,59979,33242,123709,282859,271087,147776,754049,3,15,14,29,49,
					  11,138,171,124,56,226,867,1383,1101,3818,351,1181,6130,1452,15628,
					  2464,58556,9440,64958,282859,355831,147776,754049,3,15,14,17,49,
					  11,138,171,124,56,326,867,1383,2503,3818,7245,1181,6130,8697,
					  18770,49554,69916,33242,64958,211587,91034,296177,466264,12,15,14,
					  17,49,11,101,171,124,56,326,867,1383,2503,1327,1984,10574,8159,
					  8697,18770,2464,15170,9440,32377,242821,417029,398094,754049,7,15,
					  31,17,49,8,101,171,124,56,326,867,1383,2503,1327,2999,10574,8159,
					  6436,18770,2464,15170,33242,193002,256865,417029,398094,754049,7,
					  15,31,17,49,8,101,171,231,56,326,867,1383,2503,1327,2999,3534,
					  11595,21475,18770,49554,4832,9440,193002,256865,91034,147776,
					  282852,12,15,5,17,38,8,101,171,231,56,326,867,1383,2503,1327,2999,
					  3534,8159,6436,33766,49554,4832,62850,25023,256865,91034,147776,
					  429907,12,15,5,17,38,8,101,171,90,56,326,1284,1400,2503,1327,2999,
					  3534,3436,22913,20837,2464,43064,9440,40017,122203,417029,396313,
					  390017,12,15,5,17,31,8,101,171,90,56,326,1284,1383,2503,1327,2999,
					  3534,7096,6434,20837,81,71685,9440,141605,291915,91034,578233,
					  276645,12,6,31,17,4,8,101,171,90,56,126,1284,1383,2503,1327,2999,
					  3534,7096,18497,20837,27260,4832,9440,189165,122203,299843,578233,
					  994856,12,6,13,17,4,8,101,171,90,56,326,1284,1383,429,1387,3995,
					  2898,7096,11089,20837,10681,15170,90308,189165,291915,299843,
					  578233,250142,12,6,11,17,31,18,101,171,90,56,326,1284,1383,429,
					  1387,2063,2898,7096,11089,20837,2185,15170,90308,141605,291915,
					  413548,19482,144595,12,15,11,23,64,18,101,171,90,101,326,1284,
					  1383,429,1387,2063,2898,7096,11089,20837,2185,15170,90308,189165,
					  122203,413548,620706,907454,12,15,11,23,4,18,101,171,90,101,326,
					  1284,1383,429,1387,2063,3450,7096,11089,6545,2185,27679,47904,
					  189165,25639,308300,187095,689648,12,9,11,23,4,18,101,171,90,56,
					  326,1284,1383,429,1387,2063,2141,7096,3036,6545,2185,27679,47904,
					  141605,25639,413548,620706,687580,3,13,11,23,4,18,101,171,90,101,
					  326,1284,507,429,1387,1644,2141,7096,3036,6545,2185,27679,47904,
					  141605,291803,413548,187095,687580,3,2,11,23,64,113,101,171,90,
					  101,326,563,1073,429,1387,2063,2141,7096,14208,6545,2185,60826,
					  47904,141605,245397,413548,126467,687580,3,2,13,23,45,62,101,171,
					  90,101,326,563,1073,1702,1387,2077,2141,7096,14208,6545,2185,
					  60826,47904,189165,284047,308300,241663,687580,12,2,13,23,45,62,
					  101,171,90,101,326,563,1073,1702,1387,2512,2141,7096,14208,12138,
					  18086,6187,47904,127047,245397,308300,241663,978368,7,13,13,23,45,
					  45,101,171,90,101,326,563,1073,1702,2339,2512,2141,7096,14208,
					  12138,18086,6187,47904,127047,245397,308300,241663,687580,7,11,13,
					  23,45,45,101,171,90,101,195,1010,1990,184,2339,2512,2141,7096,
					  12906,12138,18086,4264,47904,127047,245397,413548,241663,552742,
					  12,11,13,23,45,113,101,171,48,101,195,1010,1990,184,2339,2077,
					  7055,7096,12906,12138,18086,4264,47904,127047,245397,308300,
					  241663,105195,12,10,13,23,45,113,101,171,48,101,55,1010,1990,184,
					  2339,2077,7055,7096,12906,12138,18086,4264,41143,127047,245397,
					  308300,241663,942843,12,15,13,23,66,113,101,171,48,193,55,208,
					  1990,184,2339,2077,7055,7096,12906,12138,17631,4264,41143,127047,
					  245397,308300,241663,768249,12,15,14,21,66,113,116,171,48,193,55,
					  838,1990,184,2339,2077,7055,7096,12906,12138,17631,4264,41143,
					  127047,245397,308300,241663,307142,12,15,14,27,66,113,116,171,90,
					  193,55,563,507,105,2339,754,7055,7096,12906,12138,18086,45567,
					  41143,127047,94241,308300,241663,307142,12,15,14,3,66,113,116,171,
					  90,193,55,563,507,105,2339,754,7055,4377,12906,12138,18086,32269,
					  41143,127047,66575,15311,241663,307142,12,15,14,3,66,113,116,171,
					  90,193,55,563,507,105,2339,754,7055,7096,12906,12138,18086,32269,
					  41143,127047,66575,15311,241663,307142,12,15,14,3,66,113,116,171,
					  90,193,55,759,507,105,2339,754,7055,4377,7614,12138,37335,32269,
					  41143,127047,217673,15311,241663,880619,12,15,14,24,66,113,116,
					  171,90,193,55,759,507,105,2339,754,7055,4377,7614,12138,37774,
					  32269,36114,127047,217673,15311,321632,880619,3,15,14,27,66,113,
					  100,171,90,101,55,564,507,105,2339,754,7055,4377,7614,12138,37774,
					  62060,36114,127047,217673,176255,23210,880619,3,15,14,27,66,113,
					  100,171,90,101,55,759,507,105,2339,754,7055,4377,7614,12138,37774,
					  62060,36114,127047,217673,176255,23210,880619,3,6,14,17,66,113,
					  100,171,90,101,55,759,507,105,3148,754,7055,4377,5021,30483,26401,
					  62060,36114,127047,217673,23613,394484,880619,12,6,14,29,66,113,
					  100,171,90,101,55,801,507,105,3148,754,7055,5410,5021,30483,26401,
					  62060,36114,127047,217673,23613,394484,880619,7,6,14,29,66,113,
					  100,171,90,101,55,801,1073,105,3148,754,7055,5410,5021,30483,
					  26401,62060,24997,127047,217673,23613,394484,880619,7,15,14,29,66,
					  113,138,161,90,101,55,801,1073,105,3148,754,7055,4377,5021,30483,
					  26401,62060,65162,127047,217673,23613,78101,117185,12,15,14,17,66,
					  113,138,161,90,101,55,801,1073,105,3148,754,2831,4377,5021,30483,
					  26401,62060,65162,127047,217673,23613,78101,117185,12,9,14,5,66,
					  113,138,161,90,101,55,759,1073,105,3148,754,8204,4377,5021,12138,
					  26401,62060,65162,127047,217673,23613,78101,117185,12,13,14,5,66,
					  63,138,161,90,101,55,759,1073,105,3148,754,8204,4377,10145,12138,
					  26401,62060,65162,127785,217673,172210,542095,117185,12,2,14,5,66,
					  63,138,161,90,101,55,759,1073,105,3148,754,8204,4377,10145,12138,
					  26401,1803,65162,127785,217673,204328,542095,117185,12,2,31,5,66,
					  53,101,161,90,101,55,759,1073,105,3148,754,8204,4377,10145,12138,
					  26401,1803,65162,127785,217673,204328,542095,117185,12,2,31,21,66,
					  63,101,161,90,101,195,759,1073,105,3148,754,8204,4377,10145,12138,
					  26401,1803,65162,127785,217673,204328,542095,117185,12,13,5,21,11,
					  67,101,161,90,101,195,563,1073,105,3148,754,8204,4377,10145,12138,
					  26401,1803,65162,127785,217673,204328,542095,117185,12,11,5,21,66,
					  67,101,14,90,101,195,563,1073,105,3148,754,8204,4377,10145,12138,
					  26401,1803,65162,127785,217673,121626,542095,117185,7,11,5,21,66,
					  67,101,14,90,101,195,563,1073,105,3148,1097,8204,4377,10145,12138,
					  26401,1803,65162,127785,217673,121626,542095,117185,3,10,11,21,66,
					  67,101,14,90,101,195,563,1073,105,3148,1097,8204,4377,10145,12138,
					  12982,1803,65162,127785,217673,121626,542095,117185,3,10,13,21,66,
					  67,101,14,90,101,195,563,1073,105,3148,754,8204,4377,10145,12138,
					  40398,1803,65162,127785,217673,121626,542095,60731,3,15,11,21,66,
					  67,101,14,90,101,195,563,1073,105,3148,754,8204,4377,10145,12138,
					  40398,1803,65162,127785,210249,121626,542095,60731,7,15,11,21,66,
					  67,101,14,243,101,132,563,1073,105,3148,754,8204,4377,10145,12138,
					  40398,1803,65162,80822,210249,200187,542095,60731,7,15,11,21,66,
					  67,101,14,243,101,132,563,1073,105,3148,754,8204,4377,10145,12138,
					  40398,1803,47650,80822,210249,200187,542095,60731,7,15,11,21,66,
					  67,101,14,243,101,132,226,1073,105,1776,248,8204,4377,10145,12138,
					  40398,1803,47650,80822,210249,200187,542095,60731,3,15,11,21,66,
					  67,101,14,243,122,132,226,22,105,1776,754,8204,4377,10145,12138,
					  40398,1803,47650,80822,210249,200187,542095,60731,3,15,11,21,45,
					  67,101,14,243,122,132,226,22,105,1776,1097,8204,4377,10145,12138,
					  3518,51108,47650,80822,210249,200187,542095,60731,3,15,11,21,11,
					  67,101,14,243,122,132,226,22,105,3354,1097,8204,4377,10145,12138,
					  3518,51108,47650,80822,210249,121551,542095,60731,3,15,13,21,7,67,
					  101,14,243,122,132,226,22,105,3354,1097,8204,4377,10145,12138,
					  3518,51108,47650,131661,210249,121551,542095,60731,3,6,13,21,3,67,
					  101,14,243,122,132,226,22,105,3354,1097,8204,4377,10145,12138,
					  37799,51108,47650,131661,210249,248492,542095,60731,3,2,11,21,2,
					  67,101,14,243,122,132,226,22,105,925,222,8204,4377,10145,9305,
					  37799,51108,40586,131661,210249,248492,542095,60731,3,3,13,17,2,
					  51,101,14,243,122,132,226,1073,105,3354,222,8204,4377,10145,11107,
					  37799,51108,40586,131661,94453,248492,277743,178309,3,2,5,17,2,51,
					  101,14,283,122,132,226,452,105,3354,222,8204,4377,10145,11107,
					  37799,51108,40586,131661,94453,248492,277743,178309,3,3,5,17,27,
					  51,38,14,283,122,387,226,452,784,925,222,8204,4377,10145,11107,
					  37799,51108,40586,131661,94453,248492,277743,178309,3,2,5,6,5,51,
					  38,10,283,122,387,226,452,784,925,754,8204,4377,10145,11107,37799,
					  51108,40586,131661,94453,248492,457259,178309,3,2,5,17,3,51,38,10,
					  283,122,387,226,452,784,925,1982,4688,4377,10145,11107,37799,
					  51108,40586,131661,94453,248492,457259,74373,3,2,14,17,3,12,38,10,
					  283,122,387,226,452,784,925,1982,4688,4377,4544,11107,37799,51108,
					  40586,131661,94453,248492,457259,74373,3,2,13,6,5,51,38,10,283,
					  122,387,226,452,784,925,1982,4688,4377,4544,11107,37799,51108,
					  38725,131661,94453,248492,457259,74373,3,2,5,3,5,12,38,10,283,122,
					  387,226,318,784,2133,1982,2831,4377,4544,11107,4721,55315,38725,
					  131661,94453,248492,457259,74373,3,2,5,6,2,51,38,10,283,122,387,
					  226,301,784,2133,1982,2831,4377,4544,11107,4721,55315,38725,
					  131661,94453,248492,457259,74373,3,2,5,6,2,5,38,103,283,122,387,
					  226,301,784,2133,1982,2831,4377,4544,11107,4721,54140,38725,
					  131661,94453,248492,457259,74373,3,2,5,3,2,3,3,10,16,122,387,226,
					  301,784,2133,1982,2831,440,4544,11107,4721,54140,88329,131661,
					  94453,13942,457259,74373,3,2,5,3,2,3,3,10,283,101,387,226,301,784,
					  2133,1982,2831,440,8394,11107,7067,54140,88329,131661,94453,13942,
					  457259,74373,3,2,5,3,2,2,3,10,16,101,387,226,86,784,2133,1982,
					  2831,1199,8394,11107,7067,54140,88329,131661,94453,13942,457259,
					  214965,3,2,5,3,2,2,3,10,283,101,387,226,86,784,2133,1982,2831,
					  1199,8394,9305,7067,54140,88329,7114,94453,13942,457259,214965,3,
					  2,5,3,2,5,3,5,283,101,387,226,15,784,2133,1982,2831,1199,8394,
					  9305,7067,13134,88329,131661,94453,13942,457259,214965 };

	    /* System generated locals */
	    int i__1, i__2;
	    double d__1, d__2, d__3;

	    /* Local variables */
	    int i__, k;
	    int kmx;
	    double difint, varprd;
	    int intvls;


	    /*  Automatic Multidimensional Integration Subroutine */

	    /*         AUTHOR: Alan Genz */
	    /*                 Department of Mathematics */
	    /*                 Washington State University */
	    /*                 Pulman, WA 99164-3113 */
	    /*                 Email: AlanGenz@wsu.edu */

	    /*         Last Change: 12/15/00 */

	    /*  MVKBRV computes an approximation to the integral */

	    /*      1  1     1 */
	    /*     I  I ... I       F(X)  dx(NDIM)...dx(2)dx(1) */
	    /*      0  0     0 */

	    /*    F(X) is a real NF-vector of integrands. */

	    /*  It uses randomized Korobov rules. The primary references are */
	    /*   "Randomization of Number Theoretic Methods for Multiple Integration" */
	    /*    R. Cranley and T.N.L. Patterson, SIAM J Numer Anal, 13, pp. 904-14, */
	    /*  and */
	    /*   "Optimal Parameters for Multidimensional Integration", */
	    /*    P. Keast, SIAM J Numer Anal, 10, pp.831-838. */
	    /*  If there are more than 100 variables, the remaining variables are */
	    /*  integrated using the rules described in the reference */
	    /*   "On a Number-Theoretical Integration Method" */
	    /*   H. Niederreiter, Aequationes Mathematicae, 8(1972), pp. 304-11. */

	    /* **************  Parameters ******************************************** */
	    /* ***** Input parameters */
	    /*  NDIM    Number of variables, must exceed 1, but not exceed 100 */
	    /*  MINVLS  Int minimum number of function evaluations allowed. */
	    /*          MINVLS must not exceed MAXVLS.  If MINVLS < 0 then the */
	    /*          routine assumes a previous call has been made with */
	    /*          the same integrands and continues that calculation. */
	    /*  MAXVLS  Int maximum number of function evaluations allowed. */
	    /*  NF      Number of integrands, must exceed 1, but not exceed 5000 */
	    /*  FUNSUB  EXTERNALly declared user defined integrand subroutine. */
	    /*          It must have parameters ( NDIM, Z, NF, FUNVLS ), where */
	    /*          Z is a real NDIM-vector and FUNVLS is a real NF-vector. */

	    /*  ABSEPS  Required absolute accuracy. */
	    /*  RELEPS  Required relative accuracy. */
	    /* ***** Output parameters */
	    /*  MINVLS  Actual number of function evaluations used. */
	    /*  ABSERR  Maximum norm of estimated absolute accuracy of FINEST. */
	    /*  FINEST  Estimated NF-vector of values of the integrals. */
	    /*  INFORM  INFORM = 0 for normal exit, when */
	    /*                     ABSERR <= MAX(ABSEPS, RELEPS*||FINEST||) */
	    /*                  and */
	    /*                     INTVLS <= MAXCLS. */
	    /*          INFORM = 1 If MAXVLS was too small to obtain the required */
	    /*          accuracy. In this case a value FINEST is returned with */
	    /*          estimated absolute accuracy ABSERR. */
	    /* *********************************************************************** */
	    /* Parameter adjustments */
	    --finest;

	    /* Function Body */
	    *inform__ = 1;
	    intvls = 0;
	    varprd = 0.;
	    if (*minvls >= 0) {
		i__1 = *nf;
		for (k = 1; k <= i__1; ++k) {
		    finest[k] = 0.;
		    mvkbrv_varest[k - 1] = 0.;
		}
		mvkbrv_sampls = 8;
		for (i__ = std::min(*ndim,10); i__ <= 28; ++i__) {
		    mvkbrv_np = i__;
		    if (*minvls < (mvkbrv_sampls << 1) * p[i__ - 1]) {
			goto L10;
		    }
		}
		/* Computing MAX */
		i__1 = 8, i__2 = *minvls / (p[mvkbrv_np - 1] << 1);
		mvkbrv_sampls = std::max(i__1,i__2);
	    }
	L10:
	    mvkbrv_vk[0] = 1. / p[mvkbrv_np - 1];
	    k = 1;
	    i__1 = *ndim;
	    for (i__ = 2; i__ <= i__1; ++i__) {
		if (i__ <= 100) {
		    /* Computing MIN */
		    i__2 = *ndim - 1;
		    d__1 = c__[mvkbrv_np + std::min(i__2,99) * 28 - 29] * (double) k;
		    d__2 = (double) p[mvkbrv_np - 1];
		    k = (int) d_mod(&d__1, &d__2);
		    mvkbrv_vk[i__ - 1] = k * mvkbrv_vk[0];
		} else {
		    d__1 = (double) (i__ - 100) / (*ndim - 99);
		    mvkbrv_vk[i__ - 1] = (double) ((int) (p[mvkbrv_np - 1] * pow_dd(&c_b29, 
									   &d__1)));
		    d__1 = mvkbrv_vk[i__ - 1] / p[mvkbrv_np - 1];
		    mvkbrv_vk[i__ - 1] = d_mod(&d__1, &c_b24);
		}
	    }
	    i__1 = *nf;
	    for (k = 1; k <= i__1; ++k) {
		mvkbrv_finval[k - 1] = 0.;
		mvkbrv_varsqr[k - 1] = 0.;
	    }

	    i__1 = mvkbrv_sampls;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		mvkrsv_(ndim, &c__100, mvkbrv_values, &p[mvkbrv_np - 1], mvkbrv_vk, nf, mvkbrv_x, 
			mvkbrv_r__, mvkbrv_pr, mvkbrv_fs, rng);
		i__2 = *nf;
		for (k = 1; k <= i__2; ++k) {
		    difint = (mvkbrv_values[k - 1] - mvkbrv_finval[k - 1]) / i__;
		    mvkbrv_finval[k - 1] += difint;
		    /* Computing 2nd power */
		    d__1 = difint;
		    mvkbrv_varsqr[k - 1] = (i__ - 2) * mvkbrv_varsqr[k - 1] / i__ + d__1 * d__1;
		}
	    }

	    intvls += (mvkbrv_sampls << 1) * p[mvkbrv_np - 1];
	    kmx = 1;
	    i__1 = *nf;

	    for (k = 1; k <= i__1; ++k) {
		varprd = mvkbrv_varest[k - 1] * mvkbrv_varsqr[k - 1];
		finest[k] += (mvkbrv_finval[k - 1] - finest[k]) / (varprd + 1);
		if (mvkbrv_varsqr[k - 1] > 0.) {
		    mvkbrv_varest[k - 1] = (varprd + 1) / mvkbrv_varsqr[k - 1];
		}
		if ((d__1 = finest[k], std::abs(d__1)) > (d__2 = finest[kmx], std::abs(d__2))) {
		    kmx = k;
		}
	    }
	    *abserr = std::sqrt(mvkbrv_varsqr[kmx - 1] / (varprd + 1)) * 7 / 2;
	    /* Computing MAX */
	    d__2 = *abseps, d__3 = (d__1 = finest[kmx], std::abs(d__1)) * *releps;
	    if (*abserr > std::max(d__2,d__3)) {
		if (mvkbrv_np < 28) {
		    ++mvkbrv_np;
		} else {
		    /* Computing MIN */
		    i__1 = mvkbrv_sampls * 3 / 2, i__2 = (*maxvls - intvls) / (p[mvkbrv_np - 1] << 
									1);
		    mvkbrv_sampls = std::min(i__1,i__2);
		    mvkbrv_sampls = std::max(8,mvkbrv_sampls);
		}
		if (intvls + (mvkbrv_sampls << 1) * p[mvkbrv_np - 1] <= *maxvls) {
		    goto L10;
		}
	    } else {
		*inform__ = 0;
	    }
	    *minvls = intvls;

	    /*    Optimal Parameters for Lattice Rules */


	    return 0;
	} /* mvkbrv_ */

	template <class RNGType>
	/* Subroutine */ int mvkrsv_(int *ndim, int *kl, double *values, 
					  const int *prime, double *vk, int *nf, double *
					  x, double *r__, int *pr, double *fs, RNGType *rng)
	{
	    /* System generated locals */
	    int i__1, i__2;
	    double d__1;

	    /* Local variables */
	    int j, k, jp;

	    /*     For lattice rule sums */

	    /* Parameter adjustments */
	    --fs;
	    --pr;
	    --r__;
	    --x;
	    --vk;
	    --values;

	    /* Function Body */
	    i__1 = *nf;
	    for (j = 1; j <= i__1; ++j) {
		values[j] = 0.;
	    }

	    /*     Determine random shifts for each variable; scramble lattice rule */

	    i__1 = *ndim;
	    for (j = 1; j <= i__1; ++j) {
		r__[j] = mvuni_(rng);
		if (j < *kl) {
		    jp = (int) (j * r__[j] + 1);
		    if (jp < j) {
			pr[j] = pr[jp];
		    }
		    pr[jp] = j;
		} else {
		    pr[j] = j;
		}
	    }

	    /*     Compute latice rule sums */

	    i__1 = *prime;
	    for (k = 1; k <= i__1; ++k) {
		i__2 = *ndim;
		for (j = 1; j <= i__2; ++j) {
		    r__[j] += vk[pr[j]];
		    if (r__[j] > 1.) {
			--r__[j];
		    }
		    x[j] = (d__1 = r__[j] * 2 - 1, std::abs(d__1));
		}
		mvsubr_(ndim, &x[1], nf, &fs[1]);
		i__2 = *nf;
		for (j = 1; j <= i__2; ++j) {
		    values[j] += (fs[j] - values[j]) / ((k << 1) - 1);
		}
		i__2 = *ndim;
		for (j = 1; j <= i__2; ++j) {
		    x[j] = 1 - x[j];
		}
		mvsubr_(ndim, &x[1], nf, &fs[1]);
		i__2 = *nf;
		for (j = 1; j <= i__2; ++j) {
		    values[j] += (fs[j] - values[j]) / (k << 1);
		}
	    }

	    return 0;
	} /* mvkrsv_ */

	template<class RNGType>
	inline double mvuni_(RNGType *rng)
	{
	    return rng->uniform();
	} /* mvuni_ */
    
	template <class RNGType>
	double mvt(Vector& lower, Vector& upper, int& df,
		   Matrix &corr, Vector &delta, double &error, int &inform, RNGType *rng) {
	    // check inputs
	    int n = corr.cols();
	    ASSERT_MSG(2<=n && n <= 1000, "Invalid dimension");

	    VectorTemplate<int> infin = VectorTemplate<int>::Constant(n, 2);
	    for (int i=0; i<n; i++) {
		if (upper(i) == INF)
		    if (lower(i) == NEGINF)
			infin(i) = -1;
		    else
			infin(i) = 1;
		else if(lower(i) == NEGINF)
		    infin(i) = 0;
	    }

	    // handling a bug in `mvtdst' not yet fixed
	    if ((infin.array() < 0).all())
		return 1.;

	    // n must >= 2
	    double *corrF = std::unique_ptr<double>(new double[(n-1)*n/2]).get();
	    int idx  = 1;
	    int idxF = 0;
	    for (int j=0; j<n-1; j++) {
		for (int i=0; i<n-1-j; i++) {
		    corrF[idxF] = corr(idx);
		    idxF ++;
		    idx  ++;
		}
		idx += j+2;
	    }

	    // call mvtdst_
	    double value=0.;
	    int maxpts = 25000;
	    double abseps = 0.001;
	    double releps = 0;
	    
	    mvtdst_(&n, &df, lower.data(), upper.data(), infin.data(), corrF,
		    delta.data(), &maxpts, &abseps, &releps, &error, &value, &inform, rng);

	    switch(inform) {
	    case 1:
		LOG_WARNING("Completion with error > abseps");
		break;
	    case 3:
		LOG_WARNING("Covariance matrix not positive semidefinite");
		break;
	    }

	    return value;
	} // function mvt
	
    };  // class Mvt
    
}  // namespace bnc

#endif /* MVT_H */
