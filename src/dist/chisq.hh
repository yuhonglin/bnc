#ifndef CHISQ_H
#define CHISQ_H

#include <dist/brutil.hh>
#include <dist/gamma.hh>

namespace bnc {
    template<class RNGType>
    double rchisq(const int &df, RNGType *rng)
    {
	if (!R_FINITE(df) || df < 0.0) ML_ERR_return_NAN;

	return rgamma(df / 2.0, 2.0, rng);
    }

    // imitate R interfaces
    R_RFUNC_INTERFACE_1ARG(rchisq);
    
}  // namespace bnc

#endif /* CHISQ_H */
