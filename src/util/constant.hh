#ifndef CONSTANT_H
#define CONSTANT_H

#include <limits>

namespace bnc {
    const double NaN = std::numeric_limits<double>::quiet_NaN();
    const double MISS = NaN;  // missing value
    const double INF = std::numeric_limits<double>::infinity();
    const double NEGINF = -std::numeric_limits<double>::infinity();
    const double XTOL = 1e-7;
    const double DNTOL = 1e-7;  // diff normal tolerance
}  // namespace bnc

#endif /* CONSTANT_H */
