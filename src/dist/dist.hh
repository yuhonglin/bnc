#ifndef DIST_H
#define DIST_H

#include "matrix/matrix.hh"
#include "dist/unif.hh"
#include "dist/norm.hh"
#include "dist/gamma.hh"
#include "dist/beta.hh"
#include "dist/binom.hh"
#include "dist/pois.hh"
#include "dist/geom.hh"
#include "dist/hyper.hh"
#include "dist/cauchy.hh"
#include "dist/exp.hh"
#include "dist/f.hh"

namespace bnc {
    enum SCALE { NORMAL=0, LOG=1 };
    enum TAIL { UPTAIL=0, LOWTAIL=1 };
}  // namespace bnc


#endif /* DIST_H */
