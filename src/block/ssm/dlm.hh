/*
  Dynamic Linear Model based on Eigen
  Author: Honglin Yu
 */
#ifndef DLM_H
#define DLM_H

#include <type_traits>
#include <vector>
#include <iostream>

#include <util/traits.hh>
#include <matrix/matrix.hh>
#include <dist/mvnorm.hh>


namespace bnc {
    class DLM {
    private:
	std::vector<Matrix> S;  // \Sigma
	std::vector<Matrix> hS; // \hat\Sigma
	std::vector<Vector> U;  // \mu
	std::vector<Vector> hU; // \hat\mu
	std::vector<Vector> sU; // smoothing marginal mean
	std::vector<Matrix> sS; // smoothing marginal variance
	int                len;

	// for vector
	template<class T>
	inline typename T::value_type&
	nth(T& v,
	    const typename
	    std::enable_if<is_special<T, std::vector>::value,int>::type & n) {
	    return v[n];
	}
	// for element
	template<class T>
	inline T& nth(T& v, const int& n) {
	    return v;
	}
	
    public:
	DLM() {
	}

	template<class DynMatType, class ObsMatType, class DynCovType,
		 class ObsCovType>
	void filter(const Matrix & y, const DynMatType& A, const ObsMatType& C,
		    const DynCovType& Sw, const ObsCovType& Sv,
		    const Vector& m0, const Matrix& C0) {
	    const int length = y.size()+1;
	    len = y.size();
	    S.resize(length);
	    hS.resize(length);
	    U.resize(length);
	    hU.resize(length);
	    U[0] = m0;
	    S[0] = C0;
	    
	    Matrix K(nth(A,0).rows(), Sv.cols());
	    for (int i=0; i<len; i++) {		
		hU[i] = nth(A,i)*U[i];
		hS[i] = nth(A,i)*S[i]*nth(A,i).transpose() + nth(Sw,i);
		// /*
		//  * Use naive inverse. slower by may be safer
		// */
		// K = hS[i] * nth(C,i).transpose() *
		//    (nth(C,i)*hS[i]*nth(C,i).transpose() + nth(Sv,i)).inverse();
		K = (nth(C,i)*hS[i]*nth(C,i).transpose() + nth(Sv,i)).llt()
		    .solve(nth(C,i) * hS[i]).transpose();
		U[i+1] = hU[i] + K*(y.col(i)-nth(C,i)*hU[i]);
		S[i+1] = hS[i] - (K*nth(C,i)*hS[i]);
		// /*
		//  * Force S[i+1] to be symmetric (due to numerical error)
		//  * seems not needed
		// */
		// S[i+1].triangularView<Eigen::Lower>()
		// = S[i+1].transpose().triangularView<Eigen::Lower>();
	    }
	}

	template<class RNGType, class DynMatType, class ObsMatType,
		 class DynCovType, class ObsCovType>
	Matrix sample(const Matrix& y, const DynMatType& A, const ObsMatType& C,
		      const DynCovType& Sw, const ObsCovType& Sv,
		      const Vector& m0, const Matrix& C0,
		      RNGType* rng) {
	    filter(y, A, C, Sw, Sv, m0, C0);
	    Matrix ret(U[0].size(),len);
	    Matrix L(S[0].rows(),hS[0].rows());
	    Matrix Var(S[0].rows(),S[0].cols());
	    Vector E(U[0].rows());
	    ret.col(len-1) = rmvnorm(U[len], S[len], rng);
	    for (int i = len-2; i>=0; i--) {
//		L   = S[i+1]*nth(A,i).transpose()*hS[i+1].inverse(); // FIXME: use solve and Cholesky
		L   = hS[i+1].llt().solve(nth(A,i)*S[i+1]).transpose();
		E   = U[i+1] + L*(ret.col(i+1)-hU[i+1]);
		Var = S[i+1] - L*nth(A,i)*S[i+1];
                ret.col(i) = rmvnorm(E, Var, rng);
            }

	    return ret;
	}

	template<class DynMatType, class ObsMatType, class DynCovType,
		 class ObsCovType>
	void smooth(const Matrix & y, const DynMatType& A, const ObsMatType& C,
		    const DynCovType& Sw, const ObsCovType& Sv,
		    const Vector& m0, const Matrix& C0) {
	    filter(y, A, C, Sw, Sv, m0, C0);
	    sU.resize(y.size());
	    sS.resize(y.size());

	    Matrix L(S[0].rows(),hS[0].rows());
	    sU[len-1] = U[len];
	    sS[len-1] = S[len];
	    for (int i = len-2; i>=0; i--) {
		L      = S[i+1]*nth(A,i).transpose()*hS[i+1].inverse(); // FIXME: use solve and Cholesky
		sU[i]  = U[i+1] + L*(sU[i+1] - hU[i+1]);
		sS[i]  = S[i+1] + L*(sS[i+1] - hS[i+1])*L.transpose();
            }
	}


	std::vector<Vector> getFilterMean() {
	    return std::vector<Vector>(U.begin()+1,U.end());
	}

	std::vector<Matrix> getFilterCov() {
	    return std::vector<Matrix>(S.begin()+1,S.end());
	}

	std::vector<Vector> getSmoothMean() {
	    return sU;
	}

	std::vector<Matrix> getSmoothCov() {
	    return sS;
	}
	
    };  // class DLM
}  // namespace bnc

#endif /* DLM_H */
