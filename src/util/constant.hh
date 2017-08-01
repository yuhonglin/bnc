#ifndef CONSTANT_H
#define CONSTANT_H

#include <limits>
#include <cmath>

namespace bnc {
    const double NaN = std::numeric_limits<double>::quiet_NaN();
    const double MISS = NaN;  // missing value
    const double INF = std::numeric_limits<double>::infinity();
    const double NEGINF = -std::numeric_limits<double>::infinity();
    const double XTOL = 1e-20;  // min difference between x
    const double DNTOL = 1e-20;  // diff normal tolerance

    const double MCHBASE = 2; // machine base, normally 2 (binary)
    const int    DBL_MANTISSA = std::numeric_limits<double>::digits;
    const int    DBL_EPS = std::numeric_limits<double>::epsilon()/2;
    // The minimum exponent before (gradual) underflow occurs.
    const int    DBL_EMIN = int(std::log(std::numeric_limits<double>::min())/log(2)+1);
    const int    DBL_RMIN = std::numeric_limits<double>::min();
    const double DBL_BIG  = std::numeric_limits<double>::max();
    const double DBL_RTRMIN = std::sqrt(std::numeric_limits<double>::min());
    const double DBL_RTBIG = std::sqrt(std::numeric_limits<double>::max());
    const double DBL_RTEPS = std::sqrt(std::numeric_limits<double>::epsilon());
}  // namespace bnc

#endif /* CONSTANT_H */
