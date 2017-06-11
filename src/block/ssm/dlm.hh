#ifndef DLM_H
#define DLM_H

#include <matrix/matrix.hh>

namespace bnc {
    class StaticDLM
    {
	// observed_(t+1) ~ N(FF*observed_t, V)
	Matrix FF;
	Matrix V;
	// latent_(t+1) ~ N(GG*latent_t, W)
	Matrix GG;
	Matrix W;
	// theta0 ~ N(m0, C0)
	Vector m0;
	Matrix C0;
	
    public:
	void filter();
	void smooth();
	void predict();

	void sample();
	
	// mle
	void mle();
    };

    
    
}  // namespace bnc

#endif /* DLM_H */
