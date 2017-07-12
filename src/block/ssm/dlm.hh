#ifndef DLM_H
#define DLM_H

#include <type_traits>
#include <vector>
#include <iostream>

#include <matrix/matrix.hh>
#include <dist/mvnorm.hh>


namespace bnc {
    class StaticDLM {
    private:
	std::vector<Matrix> S;  // \Sigma
	std::vector<Matrix> hS; // \hat\Sigma
	std::vector<Vector> U;  // \mu
	std::vector<Vector> hU; // \hat\mu
	std::vector<Matrix> L;  // L
	Matrix         A;
	Matrix         C;
	Matrix         Sw;  // dynamic noise
	Matrix         Sv;  // observation noise
	int           len;

    public:
	StaticDLM(const int& length, const Vector m0, const Matrix& C0,
		  const Matrix& AA, const Matrix& CC,
		  const Matrix& W, const Matrix& V)
	    : S(length+1), hS(length+1), U(length+1), hU(length+1), L(length+1) {
	    U[0] = m0;
	    S[0] = C0;
	    A  = AA;
	    C  = CC;
	    Sw = W;
	    Sv = V;
	    len = length;
	}
	
	void filter(const Matrix & y) {	    
	    Matrix K(A.rows(), Sv.cols());
	    for (int i=0; i<len; i++) {		
		hU[i] = A*U[i];		
		hS[i] = A*S[i]*A.transpose() + Sw;		
		    K = hS[i] * C.transpose() *
			(C*hS[i]*C.transpose()+Sv).inverse(); // FIXME: use solve and Cholesky
     	       U[i+1] = hU[i] + K*(y.col(i)-C*hU[i]);	       
	       S[i+1] = hS[i] - K*C*hS[i];
	    }
	}

	template<class RNGType>
	Matrix sample(const Matrix& y, RNGType* rng) {
	    filter(y);
	    Matrix ret(U[0].size(),len);
	    Matrix L(S[0].rows(),hS[0].rows());
	    Matrix Var(S[0].rows(),S[0].cols());
	    Vector E(U[0].rows());
	    ret.col(len-1) = rmvnorm(U[len], S[len], rng);
	    for (int i = len-2; i>=0; i--) {
		L   = S[i+1]*A.transpose()*hS[i+1].inverse(); // FIXME: use solve and Cholesky
		E   = U[i+1] + L*(ret.col(i+1)-hU[i+1]);
		Var = S[i+1] - L*hS[i+1].inverse()*A*S[i+1];
		ret.col(i) = rmvnorm(E, Var, rng);		
	    }

	    return ret;
	}

	std::vector<Vector> getFilterMean() {
	    return std::vector<Vector>(U.begin()+1,U.end());
	}

	std::vector<Matrix> getFilterCov() {
	    return std::vector<Matrix>(S.begin()+1,S.end());
	}
	
    };  // class StaticDLM
}  // namespace bnc

#endif /* DLM_H */
