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

	/**
	 * Pick the nth element if first input is a vector
	 * otherwise return the first input.
	 */
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

	/**
	 * The Kalman filter function
	 * 
	 * \param y : the input data
	 * \param A : the dynamic transform matrix
	 * \param C : the observation matrix
	 * \param Sw : the covariance of noise in latent dynamics
	 * \param Sv : the covariance of noise in observation
	 * \param m0 : the mean of latent at time 0
	 * \param C0 : the variance of latent at time 0
	 */
	template<class DynMatType, class ObsMatType, class DynCovType,
		 class ObsCovType>
	void filter(const Matrix & y, const DynMatType& A, const ObsMatType& C,
		    const DynCovType& Sw, const ObsCovType& Sv,
		    const Vector& m0, const Matrix& C0) {
	    const int length = y.cols()+1;
	    len = y.cols();
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
		K = (nth(C,i)*hS[i]*nth(C,i).transpose() + nth(Sv,i)).llt()
		    .solve(nth(C,i) * hS[i]).transpose();
		U[i+1] = hU[i] + K*(y.col(i)-nth(C,i)*hU[i]);
		S[i+1] = hS[i] - (K*nth(C,i)*hS[i]);
	    }
	}

	/**
	 * Sample a single sample from the posterior distribution of
	 * the latent variable.
	 * NOT use it when multiple samples are needed for two reasons
	 *    1. Unnecessary filterings will be called.
	 *    2. Multiple Cholesky decomposition will be done in rmvnorm.
	 * \param y : the input data
	 * \param A : the dynamic transform matrix
	 * \param C : the observation matrix
	 * \param Sw : the covariance of noise in latent dynamics
	 * \param Sv : the covariance of noise in observation
	 * \param m0 : the mean of latent at time 0
	 * \param C0 : the variance of latent at time 0
	 * \param rng : rng
	 * \return Matrix: a matrix of latent samples. The ith column
	 *                 is a sample of latent variable at time i.
	 */
	template<class RNGType, class DynMatType, class ObsMatType,
		 class DynCovType, class ObsCovType>
	Matrix sample(const Matrix& y, const DynMatType& A, const ObsMatType& C,
		      const DynCovType& Sw, const ObsCovType& Sv,
		      const Vector& m0, const Matrix& C0,
		      RNGType* rng) {
	    filter(y, A, C, Sw, Sv, m0, C0);
	    Matrix ret(U[0].rows(),len+1);
	    Matrix L(S[0].rows(),hS[0].rows());
	    Matrix Var(S[0].rows(),S[0].cols());
	    Vector E(U[0].rows());
	    ret.col(len) = rmvnorm(U[len], S[len], rng);
	    for (int i = len-1; i>=0; i--) {
		L   = hS[i].llt().solve(nth(A,i)*S[i]).transpose();
		E   = U[i] + L*(ret.col(i+1)-hU[i]);
		Var = S[i] - L*nth(A,i)*S[i];
                ret.col(i) = rmvnorm(E, Var, rng);
            }

	    return ret;
	}

	/**
	 * The Kalman smoothing function
	 * \param y : the input data
	 * \param A : the dynamic transform matrix
	 * \param C : the observation matrix
	 * \param Sw : the covariance of noise in latent dynamics
	 * \param Sv : the covariance of noise in observation
	 * \param m0 : the mean of latent at time 0
	 * \param C0 : the variance of latent at time 0
	 * 
	 */
	template<class DynMatType, class ObsMatType, class DynCovType,
		 class ObsCovType>
	void smooth(const Matrix & y, const DynMatType& A, const ObsMatType& C,
		    const DynCovType& Sw, const ObsCovType& Sv,
		    const Vector& m0, const Matrix& C0) {
	    filter(y, A, C, Sw, Sv, m0, C0);
	    sU.resize(y.cols()+1);
	    sS.resize(y.cols()+1);

	    Matrix L(S[0].rows(),hS[0].rows());
	    sU[len] = U[len];
	    sS[len] = S[len];
	    for (int i = len-1; i>=0; i--) {
		// L      = S[i+1]*nth(A,i).transpose()*hS[i+1].inverse(); // FIXME: use solve and Cholesky
		L      = (hS[i].llt().solve(nth(A,i)*S[i])).transpose();
		sU[i]  = U[i] + L*(sU[i+1] - hU[i]);
		sS[i]  = S[i] + L*(sS[i+1] - hS[i])*L.transpose();
            }
	}

	std::vector<Vector> getFilterMean() {
	    return std::vector<Vector>(U.begin(),U.end());
	}

	std::vector<Matrix> getFilterCov() {
	    return std::vector<Matrix>(S.begin(),S.end());
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
