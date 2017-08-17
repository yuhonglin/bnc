#include <iostream>
#include <vector>

#include <rng/rng.hh>
#include <block/ssm/dlm.hh>
#include <util/rembed.hh>
#include <util/profiling.hh>

using namespace std;
using namespace bnc;

int main(int argc, char *argv[])
{

    auto R = REmbed::get_instance(argc, argv);

    int N = 500000;

    int T = 10;
    int nObs = 5;
    int nLatent = 3;
    
    bnc::MersenneTwisterRNG<bnc::AHRENS_DIETER> rng(100);
    
    Vector tmp;

    std::srand(100);

    // transform matrices
    Matrix G = Matrix::Random(nLatent,nLatent);
    Matrix F = Matrix::Random(nObs,nLatent);
    
    // variance-covariance matrices
    tmp = Vector::Random(nLatent);
    Matrix W = tmp*tmp.transpose();
    W.diagonal().array() += 0.1;
    tmp = Vector::Random(nObs);
    Matrix V = tmp*tmp.transpose();
    V.diagonal().array() += 0.1;    

    // m0 and C0
    Vector m0 = Vector::Random(3);
    tmp = Vector::Random(3);    
    Matrix C0 = tmp*tmp.transpose();
    C0.diagonal().array() += 0.1;    

    // generate some data
    Matrix ct = Matrix::Random(nObs, T);
    ct.array() += 1;

    // sample
    vector<Matrix> vecmat;
    vecmat.reserve(N);
    tic();
    DLM dlm;
    for (int i=0; i<N; i++) {
	auto sample = dlm.sample(ct, G, F, W, V, m0, C0, &rng);
	vecmat.push_back(sample);
//	if (i%100==0)
//	    cout << i << endl;
    }
    toc();
    for (int i = 0; i<nLatent; i++) {
	for (int j = 0; j<T+1; j++) {
	    std::ostringstream nms;
	    nms << "sample." << i+1 << "." << j+1;
	    Vector tmp(N);
	    for (int k = 0; k<N; k++) {
		tmp[k] = vecmat[k](i,j);
	    }
	    R.define_var(nms.str(), tmp);
	}
    }

    // expose parameters to R
#define EXPOSE(x)				\
    R.define_var(#x, x);

    EXPOSE(N);
    EXPOSE(T);
    EXPOSE(nObs);
    EXPOSE(nLatent);
    EXPOSE(ct);
    EXPOSE(G);    
    EXPOSE(F);    
    EXPOSE(W);    
    EXPOSE(V);    
    EXPOSE(m0);
    EXPOSE(C0);

#undef EXPOSE    


    // expose smoothing means and variances
    auto means = dlm.getFilterMean();
    Matrix filter_means(T+1,nLatent);
    for (int i=0; i<T+1; i++) {
	filter_means.row(i) = means[i].transpose();
    }

    R.define_var("fltmean", filter_means);
    
    auto var = dlm.getFilterCov();
    for (int i = 0; i<T+1; i++) {
	std::ostringstream nms;
	nms << "cov." << i+1;
	R.define_var(nms.str(), var[i]);
    }

    R.eval("source('testdlm.R')");
    
    R.repl();
    
    return 0;
}


