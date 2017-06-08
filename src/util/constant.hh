#ifndef CONSTANT_H
#define CONSTANT_H

#include <limits>

namespace bnc {
    const double NaN = std::numeric_limits<double>::quiet_NaN();
    const double MISS = NaN; // missing value
}  // namespace bnc

#endif /* CONSTANT_H */
