#ifndef RNG_H_
#define RNG_H_

#include <string>
#include <vector>

namespace bnc {

    /**
     * Abstract class for a psuedo random number generator
     *
     * @short Random Number Generator
     */
    struct RNG
    {
	const std::string _name;

	RNG(std::string const &name);
	virtual ~RNG();
	/**
	 * This static utility function may be used by an RNG object to coerce
	 * values in the range [0,1] to the open range (0,1)
	 */
	static double fixup(double x);
	/**
	 * Returns the name of the RNG
	 */
	std::string const &name() const;

	static double uniform(void *);
	static double normal(void *);
	static double exponential(void *);
    };

} /* namespace bnc */

#endif /* RNG_H_ */
