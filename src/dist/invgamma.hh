#ifndef INVGAMMA_H
#define INVGAMMA_H

#include <cmath>

#include <util/logger.hh>
#include <util/constant.hh>
#include <dist/brutil.hh>

extern "C" {
    double lgamma(double x);
}

namespace bnc {
    /* D function */
    double dinvgamma(const double& x, const double &shape,
		     const double&scale, const SCALE& give_log) {
	if (shape <=0 or scale <= 0) {
	    LOG_WARNING("Invalid inputs, return NaN");
	    return NaN;
	}

	const double alpha = shape;
	const double beta = scale;
	
	double log_density = alpha * std::log(beta)
	    - lgamma(alpha) - (alpha+1)*std::log(x) - (beta/x);

	if (give_log==NORMAL) {
	    return std::exp(log_density);
	} else {
	    return log_density;
	}
    }
}  // namespace bnc

#endif /* INVGAMMA_H */
