#ifndef CALLBACKS_H
#define CALLBACKS_H

#include <config.h>
#include "nmath.h"
#include <rng/BasicRNG.hh>

template <class RNGType>
inline double unif_rand(RNGType *rng)
{
    return rng->uniform();
}

template <class RNGType>
inline double exp_rand(RNGType *rng)
{
    return rng->exponential();
}

template <class RNGType>
inline double norm_rand(RNGType *rng)
{
    return rng->normal();
}


#endif /* CALLBACKS_H */
