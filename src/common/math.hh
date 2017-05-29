#ifndef MATH_H
#define MATH_H

#include <cmath>

namespace bnc {

    // matrix decomposition methods
    enum MAT_DECOMP{
	EIGEN_DECOMP,
	CHOL_DECOMP
    };

    
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

    inline int imin2(const int &x, const int &y)
    {
	return (x < y) ? x : y;
    }

    inline int imax2(const int &x, const int &y)
    {
	return (x < y) ? y : x;
    }

    inline double fsign(const double &x, const double &y)
    {
	if (isnan(x) || isnan(y))
	    return x + y;
	return ((y >= 0) ? fabs(x) : -fabs(x));
    }

}

#endif /* MATH_H */
