#ifndef DIST_H
#define DIST_H

namespace bnc {
    enum SCALE { NORMAL=0, LOG=1 };
    enum TAIL { UPTAIL=0, LOWTAIL=1 };
    enum DECOM { EIGEN_DECOM, CHOL_DECOM};
}  // namespace bnc

#include "matrix/matrix.hh"
#include "dist/unif.hh"
#include "dist/norm.hh"
#include "dist/gamma.hh"
#include "dist/beta.hh"
#include "dist/binom.hh"
#include "dist/nbinom.hh"
#include "dist/pois.hh"
#include "dist/geom.hh"
#include "dist/hyper.hh"
#include "dist/cauchy.hh"
#include "dist/exp.hh"
#include "dist/f.hh"
#include "dist/weibull.hh"
#include "dist/weibull2.hh"
#include "dist/t.hh"
#include "dist/chisq.hh"

// multivariate
//#include "mvnorm.hh"
//#include "tmvnorm.hh"
//#include "wishart.hh"

#endif /* DIST_H */
