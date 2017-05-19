/* This is MODIFIED version of R's Mathlib */

/* -*- C -*-
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 1998-2013  The R Core Team
 *  Copyright (C) 2004       The R Foundation
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation; either version 2.1 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program; if not, a copy is available at
 *  http://www.r-project.org/Licenses/
 *

 * Rmath.h  should contain ALL headers from R's C code in `src/nmath'
   -------  such that ``the Math library'' can be used by simply

   ``#include <Rmath.h> ''

   and nothing else.

   It is part of the API and supports 'standalone Rmath'.

*/
#ifndef JRMATH_H
#define JRMATH_H

/* Note that on some systems we need to include math.h before the
   defines below. */
#ifndef NO_C_HEADERS
# define __STDC_WANT_IEC_60559_TYPES_EXT__ 1
# include <math.h>
#endif

#define R_VERSION_STRING "3.1.1"

	/* Undo SGI Madness */

#ifdef __sgi
# ifdef ftrunc
#  undef ftrunc
# endif
# ifdef qexp
#  undef qexp
# endif
# ifdef qgamma
#  undef qgamma
# endif
#endif


/* ----- The following constants and entry points are part of the R API ---- */

/* 30 Decimal-place constants */
/* Computed with bc -l (scale=32; proper round) */

/* SVID & X/Open Constants */
/* Names from Solaris math.h */

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


# ifndef R_EXT_BOOLEAN_H_
/* "copy-paste" R_ext/Boolean.h if not already included: */
 #define R_EXT_BOOLEAN_H_
// #undef FALSE
// #undef TRUE
// typedef enum { FALSE = 0, TRUE } Rboolean;
#define FALSE false
#define TRUE true
#define Rboolean bool
# endif

#define bessel_i	bnc_bessel_i
#define bessel_j	bnc_bessel_j
#define bessel_k	bnc_bessel_k
#define bessel_y	bnc_bessel_y
#define bessel_i_ex	bnc_bessel_i_ex
#define bessel_j_ex	bnc_bessel_j_ex
#define bessel_k_ex	bnc_bessel_k_ex
#define bessel_y_ex	bnc_bessel_y_ex
#define beta		bnc_beta
#define choose		bnc_choose
#define dbeta		bnc_dbeta
#define dbinom		bnc_dbinom
#define dcauchy		bnc_dcauchy
#define dchisq		bnc_dchisq
#define dexp		bnc_dexp
#define dF		bnc_dF
#define dgamma		bnc_dgamma
#define dgeom		bnc_dgeom
#define dhyper		bnc_dhyper
#define digamma		bnc_digamma
#define dlnorm		bnc_dlnorm
#define dlogis		bnc_dlogis
#define dnbeta		bnc_dnbeta
#define dnbinom		bnc_dnbinom
#define dnbinom_mu	bnc_dnbinom_mu
#define dnchisq		bnc_dnchisq
#define dnf		bnc_dnf
#define dnorm4		bnc_dnorm4
#define dnt		bnc_dnt
#define dpois		bnc_dpois
#define dpsifn		bnc_dpsifn
#define dsignrank	bnc_dsignrank
#define dt		bnc_dt
#define dtukey		bnc_dtukey
#define dunif		bnc_dunif
#define dweibull	bnc_dweibull
#define dweibull2	bnc_dweibull2
#define dwilcox		bnc_dwilcox
#define fmax2		bnc_fmax2
#define fmin2		bnc_fmin2
#define fprec		bnc_fprec
#define fround		bnc_fround
#define ftrunc		bnc_ftrunc
#define fsign		bnc_fsign
#define gammafn		bnc_gammafn
#define imax2		bnc_imax2
#define imin2		bnc_imin2
#define lbeta		bnc_lbeta
#define lchoose		bnc_lchoose
#define lgammafn	bnc_lgammafn
#define lgammafn_sign	bnc_lgammafn_sign
#define lgamma1p	bnc_lgamma1p
#define log1pmx		bnc_log1pmx
#define logspace_add	bnc_logspace_add
#define logspace_sub	bnc_logspace_sub
#define pbeta		bnc_pbeta
#define pbeta_raw	bnc_pbeta_raw
#define pbinom		bnc_pbinom
#define pcauchy		bnc_pcauchy
#define pchisq		bnc_pchisq
#define pentagamma	bnc_pentagamma
#define pexp		bnc_pexp
#define pF		bnc_pF
#define pgamma		bnc_pgamma
#define pgeom		bnc_pgeom
#define phyper		bnc_phyper
#define plnorm		bnc_plnorm
#define plogis		bnc_plogis
#define pnbeta		bnc_pnbeta
#define pnbinom		bnc_pnbinom
#define pnbinom_mu     	bnc_pnbinom_mu
#define pnchisq		bnc_pnchisq
#define pnf		bnc_pnf
#define pnorm5		bnc_pnorm5
#define pnorm_both	bnc_pnorm_both
#define pnt		bnc_pnt
#define ppois		bnc_ppois
#define psignrank	bnc_psignrank
#define psigamma	bnc_psigamma
#define pt		bnc_pt
#define ptukey		bnc_ptukey
#define punif		bnc_punif
#define pythag		bnc_pythag
#define pweibull	bnc_pweibull
#define pweibull2	bnc_pweibull2
#define pwilcox		bnc_pwilcox
#define qbeta		bnc_qbeta
#define qbinom		bnc_qbinom
#define qcauchy		bnc_qcauchy
#define qchisq		bnc_qchisq
#define qchisq_appr	bnc_qchisq_appr
#define qexp		bnc_qexp
#define qF		bnc_qF
#define qgamma		bnc_qgamma
#define qgeom		bnc_qgeom
#define qhyper		bnc_qhyper
#define qlnorm		bnc_qlnorm
#define qlogis		bnc_qlogis
#define qnbeta		bnc_qnbeta
#define qnbinom		bnc_qnbinom
#define qnbinom_mu     	bnc_qnbinom_mu
#define qnchisq		bnc_qnchisq
#define qnf		bnc_qnf
#define qnorm5		bnc_qnorm5
#define qnt		bnc_qnt
#define qpois		bnc_qpois
#define qsignrank	bnc_qsignrank
#define qt		bnc_qt
#define qtukey		bnc_qtukey
#define qunif		bnc_qunif
#define qweibull	bnc_qweibull
#define qweibull2	bnc_qweibull2
#define qwilcox		bnc_qwilcox
#define rbeta		bnc_rbeta
#define rbinom		bnc_rbinom
#define rcauchy		bnc_rcauchy
#define rchisq		bnc_rchisq
#define rexp		bnc_rexp
#define rF		bnc_rF
#define rgamma		bnc_rgamma
#define rgeom		bnc_rgeom
#define rhyper		bnc_rhyper
#define rlnorm		bnc_rlnorm
#define rlogis		bnc_rlogis
#define rnbeta		bnc_rnbeta
#define rnbinom		bnc_rnbinom
#define rnchisq		bnc_rnchisq
#define rnf		bnc_rnf
#define rnorm		bnc_rnorm
#define rnt		bnc_rnt
#define rpois		bnc_rpois
#define rsignrank	bnc_rsignrank
#define rt		bnc_rt
#define rtukey		bnc_rtukey
#define runif		bnc_runif
#define rweibull	bnc_rweibull
#define rweibull2	bnc_rweibull2
#define rwilcox		bnc_rwilcox
#define sign		bnc_sign
#define tetragamma	bnc_tetragamma
#define trigamma	bnc_trigamma

#define dnorm dnorm4
#define pnorm pnorm5
#define qnorm qnorm5

#ifdef  __cplusplus
namespace bnc {
    struct RNG;
}
// typedef bnc::RNG JRNG;
extern "C" {
#else
typedef struct JRNG JRNG;
#endif
	/* R's versions with !R_FINITE checks */

double JR_pow(double x, double y);
double JR_pow_di(double, int);

	/* Normal Distribution */

double	dnorm(double, double, double, int);
double	pnorm(double, double, double, int, int);
double	qnorm(double, double, double, int, int);
void	pnorm_both(double, double *, double *, int, int);/* both tails */

	/* Uniform Distribution */

double	dunif(double, double, double, int);
double	punif(double, double, double, int, int);
double	qunif(double, double, double, int, int);

	/* Gamma Distribution */

double	dgamma(double, double, double, int);
double	pgamma(double, double, double, int, int);
double	qgamma(double, double, double, int, int);

double  log1pmx(double);
double  log1pexp(double); // <-- ../nmath/plogis.c
double  lgamma1p(double);
double  logspace_add(double, double);
double  logspace_sub(double, double);

	/* Beta Distribution */

double	dbeta(double, double, double, int);
double	pbeta(double, double, double, int, int);
double	qbeta(double, double, double, int, int);

	/* Lognormal Distribution */

double	dlnorm(double, double, double, int);
double	plnorm(double, double, double, int, int);
double	qlnorm(double, double, double, int, int);

	/* Chi-squared Distribution */

double	dchisq(double, double, int);
double	pchisq(double, double, int, int);
double	qchisq(double, double, int, int);

	/* Non-central Chi-squared Distribution */

double	dnchisq(double, double, double, int);
double	pnchisq(double, double, double, int, int);
double	qnchisq(double, double, double, int, int);

	/* F Distibution */

double	dF(double, double, double, int);
double	pF(double, double, double, int, int);
double	qF(double, double, double, int, int);

	/* Student t Distibution */

double	dt(double, double, int);
double	pt(double, double, int, int);
double	qt(double, double, int, int);

	/* Binomial Distribution */

double	dbinom(double, double, double, int);
double	pbinom(double, double, double, int, int);
double	qbinom(double, double, double, int, int);

	/* Cauchy Distribution */

double	dcauchy(double, double, double, int);
double	pcauchy(double, double, double, int, int);
double	qcauchy(double, double, double, int, int);

	/* Exponential Distribution */

double	dexp(double, double, int);
double	pexp(double, double, int, int);
double	qexp(double, double, int, int);

	/* Geometric Distribution */

double	dgeom(double, double, int);
double	pgeom(double, double, int, int);
double	qgeom(double, double, int, int);

	/* Hypergeometric Distibution */

double	dhyper(double, double, double, double, int);
double	phyper(double, double, double, double, int, int);
double	qhyper(double, double, double, double, int, int);

	/* Negative Binomial Distribution */

double	dnbinom(double, double, double, int);
double	pnbinom(double, double, double, int, int);
double	qnbinom(double, double, double, int, int);

double	dnbinom_mu(double, double, double, int);
double	pnbinom_mu(double, double, double, int, int);
double	qnbinom_mu(double, double, double, int, int);

	/* Poisson Distribution */

double	dpois(double, double, int);
double	ppois(double, double, int, int);
double	qpois(double, double, int, int);

	/* Weibull Distribution with shape-scale parameterization */

double	dweibull(double, double, double, int);
double	pweibull(double, double, double, int, int);
double	qweibull(double, double, double, int, int);

    	/* Weibull Distribution with shape-rate parameterization */

double	dweibull2(double, double, double, int);
double	pweibull2(double, double, double, int, int);
double	qweibull2(double, double, double, int, int);

	/* Logistic Distribution */

double	dlogis(double, double, double, int);
double	plogis(double, double, double, int, int);
double	qlogis(double, double, double, int, int);

	/* Non-central Beta Distribution */

double	dnbeta(double, double, double, double, int);
double	pnbeta(double, double, double, double, int, int);
double	qnbeta(double, double, double, double, int, int);

	/* Non-central F Distribution */

double  dnf(double, double, double, double, int);
double	pnf(double, double, double, double, int, int);
double	qnf(double, double, double, double, int, int);

	/* Non-central Student t Distribution */

double	dnt(double, double, double, int);
double	pnt(double, double, double, int, int);
double	qnt(double, double, double, int, int);

	/* Studentized Range Distribution */

double	ptukey(double, double, double, double, int, int);
double	qtukey(double, double, double, double, int, int);

	/* Wilcoxon Rank Sum Distribution */

double dwilcox(double, double, double, int);
double pwilcox(double, double, double, int, int);
double qwilcox(double, double, double, int, int);

	/* Wilcoxon Signed Rank Distribution */

double dsignrank(double, double, int);
double psignrank(double, double, int, int);
double qsignrank(double, double, int, int);

	/* Gamma and Related Functions */
double	gammafn(double);
double	lgammafn(double);
double	lgammafn_sign(double, int*);
void    dpsifn(double, int, int, int, double*, int*, int*);
double	psigamma(double, double);
double	digamma(double);
double	trigamma(double);
double	tetragamma(double);
double	pentagamma(double);

double	beta(double, double);
double	lbeta(double, double);

double	choose(double, double);
double	lchoose(double, double);

	/* Bessel Functions */

double	bessel_i(double, double, double);
double	bessel_j(double, double);
double	bessel_k(double, double, double);
double	bessel_y(double, double);
double	bessel_i_ex(double, double, double, double *);
double	bessel_j_ex(double, double, double *);
double	bessel_k_ex(double, double, double, double *);
double	bessel_y_ex(double, double, double *);


	/* General Support Functions */

//double 	pythag(double, double);
int	imax2(int, int);
int	imin2(int, int);
double	fmax2(double, double);
double	fmin2(double, double);
double	sign(double);
double	fprec(double, double);
double	fround(double, double);
double	fsign(double, double);
double	ftrunc(double);

double  log1pmx(double); /* Accurate log(1+x) - x, {care for small x} */
double  lgamma1p(double);/* accurate log(gamma(x+1)), small x (0 < x < 0.5) */

/* More accurate cos(pi*x), sin(pi*x), tan(pi*x)

   In future these declarations could clash with system headers if
   someone had already included math.h with
   __STDC_WANT_IEC_60559_TYPES_EXT__ defined.
   We can add a check for that via the value of
   __STDC_IEC_60559_FUNCS__ (>= 201ymmL, exact value not yet known).
*/
double cospi(double);
double sinpi(double);
double tanpi(double);

/* Compute the log of a sum or difference from logs of terms, i.e.,
 *
 *     log (exp (logx) + exp (logy))
 * or  log (exp (logx) - exp (logy))
 *
 * without causing overflows or throwing away too much accuracy:
 */
double  logspace_add(double logx, double logy);
double  logspace_sub(double logx, double logy);


/* ----------------- Private part of the header file ------------------- */

	/* old-R Compatibility */

#ifdef OLD_RMATH_COMPAT
# define snorm	norm_rand
# define sunif	unif_rand
# define sexp	exp_rand
#endif

#if !defined(MATHLIB_PRIVATE_H) /* defined by nmath.h */
/* If isnan is a macro, as C99 specifies, the C++
   math header will undefine it. This happens on OS X */
# ifdef __cplusplus
  int JR_isnancpp(double); /* in mlutils.c */
#  define ISNAN(x)     JR_isnancpp(x)
# else
#  define ISNAN(x)     (isnan(x)!=0)
# endif

# define R_FINITE(x)    JR_finite(x)
int JR_finite(double);

# ifdef _WIN32  /* not Win32 as no config information */
#  ifdef RMATH_DLL
#   define R_EXTERN extern __declspec(dllimport)
#  else
#   define R_EXTERN extern
#  endif
R_EXTERN double NA_REAL;
R_EXTERN double R_PosInf;
R_EXTERN double R_NegInf;
R_EXTERN int N01_kind;
#  undef R_EXTERN
#else
extern int N01_kind;
# endif

#endif /* MATHLIB_PRIVATE_H */

#ifdef  __cplusplus
}
#endif

#include "rbeta.hh"
#include "rbinom.hh"
#include "rcauchy.hh"
#include "rchisq.hh"
#include "rexp.hh"
#include "rf.hh"
#include "rgamma.hh"
#include "rgeom.hh"
#include "rhyper.hh"
#include "rlnorm.hh"
#include "rlogis.hh"
#include "rnbinom.hh"
#include "rnchisq.hh"
#include "rnorm.hh"
#include "rpois.hh"
#include "rt.hh"
#include "runif.hh"
#include "rweibull.hh"
#include "rweibull2.hh"

#endif /* RMATH_H */
