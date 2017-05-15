#ifndef _MARSAGLIA_RNG_H_
#define _MARSAGLIA_RNG_H_

#include <config.h>

#include <rng/RmathRNG.h>

using std::vector;

#define i2_32m1 2.328306437080797e-10/* = 1/(2^32 - 1) */

namespace jags {
namespace base {
    template<NormKind norm_kind>
    class MarsagliaRNG : public RmathRNG<norm_kind, MarsagliaRNG<norm_kind> >
    {
	unsigned int I[2];
	void fixupSeeds();
    public:
	MarsagliaRNG(unsigned int seed);
	void init(unsigned int seed);
	bool setState(std::vector<int> const &state);
	void getState(std::vector<int> &state) const;
	double uniform();
    };

}}

namespace jags {
namespace base {

    template<NormKind norm_kind>
    MarsagliaRNG<norm_kind>::MarsagliaRNG(unsigned int seed)
	: RmathRNG<norm_kind, MarsagliaRNG<norm_kind> >("base::Marsaglia-Multicarry")
    {
	init(seed);
    }

    template<NormKind norm_kind>
    void MarsagliaRNG<norm_kind>::fixupSeeds()
    {
	if (I[0] == 0) I[0] = 1;
	if (I[1] == 0) I[1] = 1;
    }

    template<NormKind norm_kind>
    void MarsagliaRNG<norm_kind>::init(unsigned int seed)
    {
	/* Initial scrambling */
	for(unsigned int j = 0; j < 50; j++)
	    seed = (69069 * seed + 1);

	for (unsigned int j = 0; j < 2; j++) {
	    seed = (69069 * seed + 1);
	    I[j] = seed;
	}
	fixupSeeds();
    }

    template<NormKind norm_kind>
    bool MarsagliaRNG<norm_kind>::setState(vector<int> const &state)
    {
	if (state.size() != 2)
	    return false;

	I[0] = static_cast<unsigned int>(state[0]);
	I[1] = static_cast<unsigned int>(state[1]);
	fixupSeeds();
	return true;
    }

    template<NormKind norm_kind>
    void MarsagliaRNG<norm_kind>::getState(vector<int> &state) const
    {
	state.clear();
	state.push_back(static_cast<int>(I[0]));
	state.push_back(static_cast<int>(I[1]));
    }

    template<NormKind norm_kind>
    double MarsagliaRNG<norm_kind>::uniform()
    {
	I[0]= 36969*(I[0] & 0177777) + (I[0]>>16);
	I[1]= 18000*(I[1] & 0177777) + (I[1]>>16);
	/* result in in [0,1) */
	return this->fixup(((I[0] << 16)^(I[1] & 0177777)) * i2_32m1); 
    }

}}


#endif /* _MARSAGLIA_RNG_H_ */
