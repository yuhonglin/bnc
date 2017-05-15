#include "MarsagliaRNG.hh"

using std::vector;

#define i2_32m1 2.328306437080797e-10/* = 1/(2^32 - 1) */

namespace bnc {
    namespace base {
        template<NormKind _N01_kind>
        MarsagliaRNG<_N01_kind>::MarsagliaRNG(unsigned int seed)
            : RmathRNG<_N01_kind>("base::Marsaglia-Multicarry")
        {
            init(seed);
        }

        template<NormKind _N01_kind>
        void MarsagliaRNG<_N01_kind>::fixupSeeds()
        {
            if (I[0] == 0) I[0] = 1;
            if (I[1] == 0) I[1] = 1;
        }

        template<NormKind _N01_kind>
        void MarsagliaRNG<_N01_kind>::init(unsigned int seed)
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

        template<NormKind _N01_kind>
        bool MarsagliaRNG<_N01_kind>::setState(vector<int> const &state)
        {
            if (state.size() != 2)
                return false;

            I[0] = static_cast<unsigned int>(state[0]);
            I[1] = static_cast<unsigned int>(state[1]);
            fixupSeeds();
            return true;
        }

        template<NormKind _N01_kind>
        void MarsagliaRNG<_N01_kind>::getState(vector<int> &state) const
        {
            state.clear();
            state.push_back(static_cast<int>(I[0]));
            state.push_back(static_cast<int>(I[1]));
        }

        template<NormKind _N01_kind>
        double MarsagliaRNG<_N01_kind>::uniform()
        {
            I[0]= 36969*(I[0] & 0177777) + (I[0]>>16);
            I[1]= 18000*(I[1] & 0177777) + (I[1]>>16);
            /* result in in [0,1) */
            return fixup(((I[0] << 16)^(I[1] & 0177777)) * i2_32m1); 
        }

    }
}
