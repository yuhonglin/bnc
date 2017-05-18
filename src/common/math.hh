#ifndef MATH_H
#define MATH_H

#include <cmath>

namespace bnc {

    inline double fmin2(const double& x, const double& y) {
	if (isnan(x) || isnan(y))
	    return x + y;
	return (x < y) ? x : y;
    }

    inline double fmax2(const double& x, const double& y)
    {
	if (isnan(x) || isnan(y))
	    return x + y;
	return (x < y) ? y : x;
    }

}

#endif /* MATH_H */
