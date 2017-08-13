#ifndef _SUPER_DUPER_RNG_H_
#define _SUPER_DUPER_RNG_H_

#include <vector>

#include <rng/RmathRNG.hh>

namespace bnc {

    template<NormKind norm_kind>
    class SuperDuperRNG : public RmathRNG<norm_kind, SuperDuperRNG<norm_kind> >
    {
	unsigned int I[2];
	void fixupSeeds();
    public:
	SuperDuperRNG(unsigned int seed);
	SuperDuperRNG(const SuperDuperRNG<norm_kind>& rng);
	double uniform();
	void init(unsigned int seed);
	bool setState(std::vector<int> const &state);
	void getState(std::vector<int> &state) const;

	template <class Archive>
	void serialize( Archive & ar ) {
	    ar( I[0] );
	    ar( I[1] );
	}	

    };

}

#endif /* _SUPER_DUPER_RNG_H_ */


#define I1 I[0]
#define I2 I[1]

#define i2_32m1 2.328306437080797e-10/* = 1/(2^32 - 1) */

using std::vector;

namespace bnc {

    template<NormKind norm_kind>
    SuperDuperRNG<norm_kind>::SuperDuperRNG(unsigned int seed)
	: RmathRNG<norm_kind, SuperDuperRNG<norm_kind> >("base::Super-Duper")
    {
	init(seed);    
    }

    
    template<NormKind norm_kind>
    SuperDuperRNG<norm_kind>::SuperDuperRNG(const SuperDuperRNG<norm_kind>& rng)
	: RmathRNG<norm_kind, SuperDuperRNG<norm_kind> >(rng._name)
    {
	I[0] = rng.I[0];
	I[1] = rng.I[1];
    }
    
    
    template<NormKind norm_kind>
    void SuperDuperRNG<norm_kind>::fixupSeeds()
    {
	if(I1 == 0) I1 = 1;
	/* I2 = Congruential: must be ODD */
	I2 |= 1;
    }

    template<NormKind norm_kind>
    double SuperDuperRNG<norm_kind>::uniform()
    {
	/* This is Reeds et al (1984) implementation;
	 * modified using __unsigned__	seeds instead of signed ones
	 */
	I1 ^= ((I1 >> 15) & 0377777); /* Tausworthe */
	I1 ^= I1 << 17;
	I2 *= 69069;		/* Congruential */
	return this->fixup((I1^I2) * i2_32m1); /* in [0,1) */
    }


    template<NormKind norm_kind>
    void SuperDuperRNG<norm_kind>::init(unsigned int seed)
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
    bool SuperDuperRNG<norm_kind>::setState(vector<int> const &state)
    {
	if (state.size() != 2) 
	    return false;
  
	for (unsigned int j = 0; j < 2; j++) {
	    I[j] = static_cast<unsigned int>(state[j]);
	}
	fixupSeeds();
	return true;
    }


    template<NormKind norm_kind>
    void SuperDuperRNG<norm_kind>::getState(vector<int> &state) const
    {
	state.clear();
	for (unsigned int j = 0; j < 2; j++) {
	    state.push_back(static_cast<int>(I[j]));
	}
    }

}
