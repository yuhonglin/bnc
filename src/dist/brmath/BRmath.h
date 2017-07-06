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
#ifndef BRMATH_H
#define BRMATH_H

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

#define bessel_i	br_bessel_i
#define bessel_j	br_bessel_j
#define bessel_k	br_bessel_k
#define bessel_y	br_bessel_y
#define bessel_i_ex	br_bessel_i_ex
#define bessel_j_ex	br_bessel_j_ex
#define bessel_k_ex	br_bessel_k_ex
#define bessel_y_ex	br_bessel_y_ex
#define beta		br_beta
#define choose		br_choose
#define dbeta		br_dbeta
#define dbinom		br_dbinom
#define dcauchy		br_dcauchy
#define dchisq		br_dchisq
#define dexp		br_dexp
#define dF		br_dF
#define dgamma		br_dgamma
#define dgeom		br_dgeom
#define dhyper		br_dhyper
#define digamma		br_digamma
#define dlnorm		br_dlnorm
#define dlogis		br_dlogis
#define dnbeta		br_dnbeta
#define dnbinom		br_dnbinom
#define dnbinom_mu	br_dnbinom_mu
#define dnchisq		br_dnchisq
#define dnf		br_dnf
#define dnorm4		br_dnorm4
#define dnt		br_dnt
#define dpois		br_dpois
#define dpsifn		br_dpsifn
#define dsignrank	br_dsignrank
#define dt		br_dt
#define dtukey		br_dtukey
#define dunif		br_dunif
#define dweibull	br_dweibull
#define dweibull2	br_dweibull2
#define dwilcox		br_dwilcox
#define fmax2		br_fmax2
#define fmin2		br_fmin2
#define fprec		br_fprec
#define fround		br_fround
#define ftrunc		br_ftrunc
#define fsign		br_fsign
#define gammafn		br_gammafn
#define imax2		br_imax2
#define imin2		br_imin2
#define lbeta		br_lbeta
#define lchoose		br_lchoose
#define lgammafn	br_lgammafn
#define lgammafn_sign	br_lgammafn_sign
#define lgamma1p	br_lgamma1p
#define log1pmx		br_log1pmx
#define logspace_add	br_logspace_add
#define logspace_sub	br_logspace_sub
#define pbeta		br_pbeta
#define pbeta_raw	br_pbeta_raw
#define pbinom		br_pbinom
#define pcauchy		br_pcauchy
#define pchisq		br_pchisq
#define pentagamma	br_pentagamma
#define pexp		br_pexp
#define pF		br_pF
#define pgamma		br_pgamma
#define pgeom		br_pgeom
#define phyper		br_phyper
#define plnorm		br_plnorm
#define plogis		br_plogis
#define pnbeta		br_pnbeta
#define pnbinom		br_pnbinom
#define pnbinom_mu     	br_pnbinom_mu
#define pnchisq		br_pnchisq
#define pnf		br_pnf
#define pnorm5		br_pnorm5
#define pnorm_both	br_pnorm_both
#define pnt		br_pnt
#define ppois		br_ppois
#define psignrank	br_psignrank
#define psigamma	br_psigamma
#define pt		br_pt
#define ptukey		br_ptukey
#define punif		br_punif
#define pythag		br_pythag
#define pweibull	br_pweibull
#define pweibull2	br_pweibull2
#define pwilcox		br_pwilcox
#define qbeta		br_qbeta
#define qbinom		br_qbinom
#define qcauchy		br_qcauchy
#define qchisq		br_qchisq
#define qchisq_appr	br_qchisq_appr
#define qexp		br_qexp
#define qF		br_qF
#define qgamma		br_qgamma
#define qgeom		br_qgeom
#define qhyper		br_qhyper
#define qlnorm		br_qlnorm
#define qlogis		br_qlogis
#define qnbeta		br_qnbeta
#define qnbinom		br_qnbinom
#define qnbinom_mu     	br_qnbinom_mu
#define qnchisq		br_qnchisq
#define qnf		br_qnf
#define qnorm5		br_qnorm5
#define qnt		br_qnt
#define qpois		br_qpois
#define qsignrank	br_qsignrank
#define qt		br_qt
#define qtukey		br_qtukey
#define qunif		br_qunif
#define qweibull	br_qweibull
#define qweibull2	br_qweibull2
#define qwilcox		br_qwilcox
#define rbeta		br_rbeta
#define rbinom		br_rbinom
#define rcauchy		br_rcauchy
#define rchisq		br_rchisq
#define rexp		br_rexp
#define rF		br_rF
#define rgamma		br_rgamma
#define rgeom		br_rgeom
#define rhyper		br_rhyper
#define rlnorm		br_rlnorm
#define rlogis		br_rlogis
#define rnbeta		br_rnbeta
#define rnbinom		br_rnbinom
#define rnchisq		br_rnchisq
#define rnf		br_rnf
#define rnorm		br_rnorm
#define rnt		br_rnt
#define rpois		br_rpois
#define rsignrank	br_rsignrank
#define rt		br_rt
#define rtukey		br_rtukey
#define runif		br_runif
#define rweibull	br_rweibull
#define rweibull2	br_rweibull2
#define rwilcox		br_rwilcox
#define sign		br_sign
#define tetragamma	br_tetragamma
#define trigamma	br_trigamma

#define dnorm dnorm4
#define pnorm pnorm5
#define qnorm qnorm5

#ifdef  __cplusplus
namespace bnc {
    struct RNG;
}
// typedef bnc::RNG BRNG;
extern "C" {
#else
typedef struct BRNG BRNG;
#endif
	/* R's versions with !R_FINITE checks */

double BR_pow(double x, double y);
double BR_pow_di(double, int);

	/* Normal Distribution */

double	dnorm(double, double, double, int);
double	pnorm(double, const double&, const double&, const int&, const int&);
double	qnorm(double, double, double, int, int);
void	pnorm_both(const double&, double *, double *, const int&, const int&);/* both tails */

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
  int BR_isnancpp(double); /* in mlutils.c */
#  define ISNAN(x)     BR_isnancpp(x)
# else
#  define ISNAN(x)     (isnan(x)!=0)
# endif

# define R_FINITE(x)    BR_finite(x)
int BR_finite(double);

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

#endif /* RMATH_H */
