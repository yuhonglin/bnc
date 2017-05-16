#ifndef NORMAL_H
#define NORMAL_H

#include <rng/>
#include <dist/JRmath.h>


namespace bnc {

    double	dnorm(double, double, double, int) {
        return jags::dnorm(double, double, double, int);
    };
    double	pnorm(double, double, double, int, int);
    double	qnorm(double, double, double, int, int);

    /** Since some samples are arrays, to make the interface 
        identical, we do not use return values.
     */
    template<class JRNG>
    void rnorm(double&, double, JRNG*);
    
}


#endif /* NORMAL_H */
