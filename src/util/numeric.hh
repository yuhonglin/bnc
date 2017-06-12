#ifndef NUMERIC_H
#define NUMERIC_H

#include <cmath>



namespace bnc {
    inline bool eq(const double&a, const double&b, const double& tol) {
	return std::abs(a-b)<=tol ? true : false;
    }
    inline bool le(const double&a, const double&b, const double& tol) {
	return a-b<=tol ? true : false;
    }    
}  // namespace bnc

#endif /* NUMERIC_H */
