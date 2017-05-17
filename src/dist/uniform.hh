#ifndef UNIFORM_H
#define UNIFORM_H

#include <limits>
#include <cmath>

#include <dist/dist.hh>
#include <rng/rng.hh>
#include <util/logger.hh>

namespace bnc {
    
    /* Density function */
    template<SCALE_P scale_p=NOR_P>
    double dunif(const double& x, const double &min, const double &max) 
    {
	if (std::isnan(x) || std::isnan(min) || std::isnan(max)) {
	    LOG_WARNING("inputs contain NaN, return NaN.");
	    return std::numeric_limits<double>::quiet_NaN();
	}

	if (min > max) {
	    LOG_WARNING("min > max, return NaN");
	    return std::numeric_limits<double>::quiet_NaN();
	}

	if (min<=x && x<=max) {
	    return scale_p==NOR_P? 1./(max-min) : -std::log(max-min);
	} else {
	    return scale_p==NOR_P? -std::numeric_limits<double>::infinity() \
		: -std::log(max-min);
	}
    }

    double dunif(const double& x) 
    {
	if (std::isnan(x)) {
	    LOG_WARNING("x is NaN, return NaN.")
		return std::numeric_limits<double>::quiet_NaN();
	}
	    
	if (x>=0. && x<=1.) {
	    return 1.;
	} else {
	    return 0.;
	}
    }

    /* Probability function */
    template<TAIL_P tail_p=TAIL_LOWER, SCALE_P scale_p=NOR_P>
    double punif(const double& q, const double &min, const double &max) {
	if (std::isnan(q) || std::isnan(min) || std::isnan(max)) {
	    LOG_WARNING("inputs contain NaN, return NaN.")
		return std::numeric_limits<double>::quiet_NaN();
	}

	if (min > max) {
	    LOG_WARNING("min > max, return NaN");
	    return std::numeric_limits<double>::quiet_NaN();
	}
	if (tail_p==TAIL_LOWER) {
	    if (q < min) return 0.;
	    else if (q >= max) return 1.;
	    else return (q-min)/(max-min);
	} else {
	    if (q > max) return 0.;
	    else if (q <= min) return 1.;
	    else return (max-q)/(max-min);
	}
    }

    /* Quantile function */
    template<TAIL_P tail_p=TAIL_LOWER, SCALE_P scale_p=NOR_P>
    double qunif(const double& p, const double &min, const double &max) {
	if (std::isnan(p) || std::isnan(min) || std::isnan(max)) {
	    LOG_WARNING("invalid inputs, return NaN");
	    return std::numeric_limits<double>::quiet_NaN();
	}
	if (p > 1. || p < 0.) {
	    LOG_WARNING("invalid probability, return NaN");
	    return std::numeric_limits<double>::quiet_NaN();
	}

	if (tail_p == TAIL_LOWER) {
	    return min + (scale_p==NOR_P? p:std::exp(p))*(max-min);
	} else {
	    return max - (scale_p==NOR_P? p:std::exp(p))*(max-min);
	}
    }

    template<TAIL_P tail_p=TAIL_LOWER, SCALE_P scale_p=NOR_P>
    double qunif(const double& p) {
	if (p > 1. || p < 0. || std::isnan(p)) {
	    LOG_WARNING("invalid probability, return NaN");
	    return std::numeric_limits<double>::quiet_NaN();
	}
	if (tail_p == TAIL_LOWER) {
	    return 0. + (scale_p==NOR_P? p:std::exp(p));
	} else {
	    return 1. - (scale_p==NOR_P? p:std::exp(p));
	}
    }

    /* Sampling function */
    template<class RNGType>
    void runif(double *x, RNGType* rng, const double& min, const double& max) {
	if (rng==nullptr || std::isnan(min) || std::isnan(max)) {
	    LOG_WARNING("invalid inputs, return NaN.");
	    (*x) = std::numeric_limits<double>::quiet_NaN();
	    return;
	}

	if (min > max) {
	    LOG_WARNING("min > max, return NaN");
	    (*x) = std::numeric_limits<double>::quiet_NaN();
	    return;
	}
	(*x) = rng->uniform() * (max-min) + min;
    }

    template<class RNGType>
    void runif(double *x, RNGType* rng) {
	if (rng==nullptr) {
	    LOG_WARNING("rng is NULL, return NaN.");
	    (*x) = std::numeric_limits<double>::quiet_NaN();
	    return;
	}
	(*x) = rng->uniform();
    }
	
}  // namespace bnc

#endif /* UNIFORM_H */
