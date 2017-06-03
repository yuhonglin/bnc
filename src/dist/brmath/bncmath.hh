namespace brmath {

	/* Normal Distribution */

double	br_dnorm4(double, double, double, int);
double	br_pnorm(double, double, double, int, int);
double	br_qnorm(double, double, double, int, int);
void	br_pnorm_both(double, double *, double *, int, int);/* both tails */

	/* Uniform Distribution */

double	br_dunif(double, double, double, int);
double	br_punif(double, double, double, int, int);
double	br_qunif(double, double, double, int, int);
	
	/* Gamma Distribution */

double	br_dgamma(double, double, double, int);
double	br_pgamma(double, double, double, int, int);
double	br_qgamma(double, double, double, int, int);
	
double  br_log1pmx(double);
double  br_log1pexp(double); // <-- ../nmath/plogis.c
double  br_lgamma1p(double);
double  br_logspace_add(double, double);
double  br_logspace_sub(double, double);

	/* Beta Distribution */

double	br_dbeta(double, double, double, int);
double	br_pbeta(double, double, double, int, int);
double	br_qbeta(double, double, double, int, int);

	/* Lognormal Distribution */

double	br_dlnorm(double, double, double, int);
double	br_plnorm(double, double, double, int, int);
double	br_qlnorm(double, double, double, int, int);

	/* Chi-squared Distribution */

double	br_dchisq(double, double, int);
double	br_pchisq(double, double, int, int);
double	br_qchisq(double, double, int, int);

	/* Non-central Chi-squared Distribution */

double	br_dnchisq(double, double, double, int);
double	br_pnchisq(double, double, double, int, int);
double	br_qnchisq(double, double, double, int, int);

	/* F Distibution */

double	br_dF(double, double, double, int);
double	br_pF(double, double, double, int, int);
double	br_qF(double, double, double, int, int);

	/* Student t Distibution */

double	br_dt(double, double, int);
double	br_pt(double, double, int, int);
double	br_qt(double, double, int, int);

	/* Binomial Distribution */

double	br_dbinom(double, double, double, int);
double	br_pbinom(double, double, double, int, int);
double	br_qbinom(double, double, double, int, int);

	/* Cauchy Distribution */

double	br_dcauchy(double, double, double, int);
double	br_pcauchy(double, double, double, int, int);
double	br_qcauchy(double, double, double, int, int);

	/* Exponential Distribution */

double	br_dexp(double, double, int);
double	br_pexp(double, double, int, int);
double	br_qexp(double, double, int, int);

	/* Geometric Distribution */

double	br_dgeom(double, double, int);
double	br_pgeom(double, double, int, int);
double	br_qgeom(double, double, int, int);

	/* Hypergeometric Distibution */

double	br_dhyper(double, double, double, double, int);
double	br_phyper(double, double, double, double, int, int);
double	br_qhyper(double, double, double, double, int, int);

	/* Negative Binomial Distribution */

double	br_dnbinom(double, double, double, int);
double	br_pnbinom(double, double, double, int, int);
double	br_qnbinom(double, double, double, int, int);

double	br_dnbinom_mu(double, double, double, int);
double	br_pnbinom_mu(double, double, double, int, int);
double	br_qnbinom_mu(double, double, double, int, int);

	/* Poisson Distribution */

double	br_dpois(double, double, int);
double	br_ppois(double, double, int, int);
double	br_qpois(double, double, int, int);

	/* Weibull Distribution with shape-scale parameterization */

double	br_dweibull(double, double, double, int);
double	br_pweibull(double, double, double, int, int);
double	br_qweibull(double, double, double, int, int);

    	/* Weibull Distribution with shape-rate parameterization */

double	br_dweibull2(double, double, double, int);
double	br_pweibull2(double, double, double, int, int);
double	br_qweibull2(double, double, double, int, int);

	/* Logistic Distribution */

double	br_dlogis(double, double, double, int);
double	br_plogis(double, double, double, int, int);
double	br_qlogis(double, double, double, int, int);

	/* Non-central Beta Distribution */

double	br_dnbeta(double, double, double, double, int);
double	br_pnbeta(double, double, double, double, int, int);
double	br_qnbeta(double, double, double, double, int, int);

	/* Non-central F Distribution */

double  br_dnf(double, double, double, double, int);
double	br_pnf(double, double, double, double, int, int);
double	br_qnf(double, double, double, double, int, int);

	/* Non-central Student t Distribution */

double	br_dnt(double, double, double, int);
double	br_pnt(double, double, double, int, int);
double	br_qnt(double, double, double, int, int);

	/* Studentized Range Distribution */

double	br_ptukey(double, double, double, double, int, int);
double	br_qtukey(double, double, double, double, int, int);
	
	/* Wilcoxon Rank Sum Distribution */

double br_dwilcox(double, double, double, int);
double br_pwilcox(double, double, double, int, int);
double br_qwilcox(double, double, double, int, int);

	/* Wilcoxon Signed Rank Distribution */

double br_dsignrank(double, double, int);
double br_psignrank(double, double, int, int);
double br_qsignrank(double, double, int, int);

	/* Gamma and Related Functions */
double	br_gammafn(double);
double	br_lgammafn(double);
double	br_lgammafn_sign(double, int*);
void    br_dpsifn(double, int, int, int, double*, int*, int*);
double	br_psigamma(double, double);
double	br_digamma(double);
double	br_trigamma(double);
double	br_tetragamma(double);
double	br_pentagamma(double);

double	br_beta(double, double);
double	br_lbeta(double, double);
	
double	br_choose(double, double);
double	br_lchoose(double, double);
	
	/* Bessel Functions */

double	br_bessel_i(double, double, double);
double	br_bessel_j(double, double);
double	br_bessel_k(double, double, double);
double	br_bessel_y(double, double);
double	br_bessel_i_ex(double, double, double, double *);
double	br_bessel_j_ex(double, double, double *);
double	br_bessel_k_ex(double, double, double, double *);
double	br_bessel_y_ex(double, double, double *);
	
}
