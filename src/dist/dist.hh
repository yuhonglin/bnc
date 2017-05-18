/**
   The common header file of distributions, it includes
   - Common types definition
   - Error handling
*/

#ifndef DIST_H
#define DIST_H
namespace bnc {

/* prability types */
    enum SCALE_P {
	NOR_P, // normal scale
	LOG_P // log scale
    };

/* tail type, used in ```pdist``` functions */
    enum TAIL_P {
	TAIL_LOWER = 0,  // lower
	TAIL_UPPER = 1,  // upper
	TAIL_BOTH  = 2   // both
    };

}  // namespace bnc

#endif /* DIST_H */
