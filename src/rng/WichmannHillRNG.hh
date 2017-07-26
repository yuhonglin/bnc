#ifndef _WICHMANN_HILL_RNG_H_
#define _WICHMANN_HILL_RNG_H_

#include <rng/RmathRNG.hh>

namespace bnc {
    template<NormKind norm_kind>
    class WichmannHillRNG : public RmathRNG<norm_kind, WichmannHillRNG<norm_kind> >
    {
	unsigned int I[3];
	void fixupSeeds();
    public:
	WichmannHillRNG(unsigned int seed);
	WichmannHillRNG(const WichmannHillRNG<norm_kind>& rng);	
	double uniform();
	void init(unsigned int seed);
	bool setState(std::vector<int> const &state);
	void getState(std::vector<int> &state) const;
    };

}

#endif /* _WICHMANN_HILL_RNG_H_ */

using std::vector;

namespace bnc {
        
    template<NormKind norm_kind>
    void WichmannHillRNG<norm_kind>::fixupSeeds()
    {
	I[0] = I[0] % 30269; I[1] = I[1] % 30307; I[2] = I[2] % 30323;

	/* map values equal to 0 mod modulus to 1. */
	if(I[0] == 0) I[0] = 1;
	if(I[1] == 0) I[1] = 1;
	if(I[2] == 0) I[2] = 1;
    }

    template<NormKind norm_kind>
    WichmannHillRNG<norm_kind>::WichmannHillRNG(unsigned int seed)
	: RmathRNG<norm_kind, WichmannHillRNG<norm_kind> >("base::Wichmann-Hill")
    {
	init(seed);
    }

    template<NormKind norm_kind>
    WichmannHillRNG<norm_kind>::WichmannHillRNG(const WichmannHillRNG<norm_kind>& rng)
	: RmathRNG<norm_kind, WichmannHillRNG<norm_kind> >(rng._name)
    {
	I[0] = rng.I[0];
	I[1] = rng.I[1];
	I[2] = rng.I[2];		
    }
    
    template<NormKind norm_kind>
    double WichmannHillRNG<norm_kind>::uniform()
    {
	I[0] = I[0] * 171 % 30269;
	I[1] = I[1] * 172 % 30307;
	I[2] = I[2] * 170 % 30323;
	double value = I[0] / 30269.0 + I[1] / 30307.0 + I[2] / 30323.0;
	return this->fixup(value - (int) value); /* in [0,1) */
    }

    template<NormKind norm_kind>
    void WichmannHillRNG<norm_kind>::init(unsigned int seed)
    {
	/* Initial scrambling */
	for(unsigned int j = 0; j < 50; j++)
	    seed = (69069 * seed + 1);
  
  
	for(unsigned int j = 0; j < 3; ++j) {
	    seed = (69069 * seed + 1);
	    I[j] = seed;
	}
	this->fixupSeeds();
    }

    template<NormKind norm_kind>
    bool WichmannHillRNG<norm_kind>::setState(vector<int> const &state)
    {
	if (state.size() != 3)
	    return false;

	for (unsigned int i = 0; i < 3; ++i) {
	    I[i] = static_cast<unsigned int>(state[i]);
	}
	this->fixupSeeds();
	return true;
    }

    template<NormKind norm_kind>
    void WichmannHillRNG<norm_kind>::getState(vector<int> &state) const
    {
	state.clear();
	for (unsigned int i = 0; i < 3; ++i) {
	    state.push_back(static_cast<int>(I[i]));
	}
    }

}
