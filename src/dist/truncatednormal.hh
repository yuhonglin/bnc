#ifndef TRUNCATEDNORMAL_H
#define TRUNCATEDNORMAL_H

#include <cmath>
#include <vector>
#include <memory>
#include <numeric>
#include <algorithm>

#include <util/logger.hh>
#include <util/constant.hh>
#include <common/math.hh>
#include <matrix/matrix.hh>
#include <dist/norm.hh>
#include <dist/mvnorm.hh>
#include <dist/chisq.hh>
#include <dist/mvnorm.hh>

#define TRUNCATEDNORMAL_IDENTICAL

namespace bnc {
    namespace truncatednormal {
	// compute ln(P(a<X<b)), X~N(0,1)
	// I do not know why the old one computes
	// this way?
	template <class T>
	typename enable_if<is_class<T>::value,Vector>::type
	lnNpr(T &a, T &b) {
	    double mu = 0;
	    double s = 1;
	    int t = 1;
	    int f = 0;
	    return (pnorm(b, mu, s, t, f) - pnorm(a, mu, s, t, f))
		.array().log();
	}

	double lnNpr(double& a, double& b) {
	    return std::log(pnorm(b, 0, 1, true, false) - pnorm(a, 0, 1, true, false));
	}

	// samples a column vector of length=length(l)=length(u)
	// from the standard multivariate normal distribution,
	// truncated over the region [l,u], where l>0 and
	// l and u are column vectors;
	// uses acceptance-rejection from Rayleigh distr;
	// method due to Marsaglia (1964);

	template<class T>
	inline void forward_shift_iter(const T& i1, const T& i2, const int& n) {
	    auto s = i1 - n;
	    auto i = i1;
	    while (i!=i2) {
		*s = *i;
		s++;
		i++;
	    }
	}

	template<class T>
	void remove_by_index(T& data, const std::vector<int>& deleteIndices)
	{
	    if (deleteIndices.size() == 0)
		return;
	    int i = 0;
	    for (; i<deleteIndices.size()-1; i++) {
		forward_shift_iter(data.begin()+deleteIndices[i]+1, data.begin()+deleteIndices[i+1], i+1);
	    }
	    forward_shift_iter(data.begin()+deleteIndices[i]+1, data.end(), deleteIndices.size());
	    data.erase(data.end()-deleteIndices.size(), data.end());
	}

	template <class RNGType>
	Vector ntail(Vector& l, Vector& u, RNGType* rng) {
	    Vector x(l.size());
	    int n = l.size();
	    Vector c = l.array().square()/2.;
	    Vector f = c.array() - u.array().square()/2.;
	    for (int i = 0; i<f.size(); i++)
		f(i) = std::expm1(f(i));
	    // sample using Rayleigh
	    x = c.array() - (runif(n, 0., 1., rng).array()*f.array()+1).log();
	    // keep list of rejected
	    std::vector<int> I; I.reserve(l.size());
	    Vector tmp = runif(n, 0., 1., rng);
	    for (int i=0; i<n; i++)
		if (tmp(i)*tmp(i)*x(i)>c(i)) I.push_back(i);
	    int d = I.size();
	    while(d>0) { // while there are rejections
		Vector cy(I.size());
		for (int i=0; i<I.size(); i++) cy(i) = c(I[i]);
		Vector uniftmp = runif(d,0.,1.,rng);
		Vector y(I.size());
		for (int i=0; i<I.size(); i++) y(i) = cy(i) - std::log(1+uniftmp(i)*f(I[i]));
		std::vector<int> idx; idx.reserve(I.size());
		Vector uniftmp1 = runif(d,0.,1., rng);
		for (int i=0; i<I.size(); i++)
		    if (uniftmp1(i)*uniftmp1(i)*y(i) < cy(i))
			idx.push_back(i);
		for (int i=0; i<idx.size(); i++)
		    x(I[idx[i]]) = y(idx[i]);
		// remove accepted from list
		remove_by_index(I,idx);
		d = I.size();
	    }
	    return (2*x.array()).sqrt();
	}	

	template<class RNGType>	
	Vector trnd(const Vector& l, const Vector& u, RNGType *rng) {
	    int n = l.size();
	    Vector x = rnorm(n, 0., 1., rng);
	    // keep list of rejected
	    std::vector<int> I; I.reserve(l.size());
	    for (int i=0; i<l.size(); i++)
		if (x(i)<l(i) || x(i)>u(i)) I.push_back(i);
	    int d = I.size();
	    while (d>0) {
		Vector ly(I.size());
		Vector uy(I.size());
		for (int i=0; i<I.size(); i++) {
		    ly(i) = l(I[i]);
		    uy(i) = u(I[i]);
		}
		int tmpfixme = ly.size();
		Vector y = rnorm(tmpfixme, 0., 1., rng);
		std::vector<int> idx; idx.reserve(I.size());
		for(int i=0; i<ly.size(); i++)
		    if (y(i) > ly(i) && y(i) < uy(i))
			idx.push_back(i);
		for(int i=0; i<idx.size(); i++)
		    x(I[idx[i]]) = y[idx[i]];
		// remove accepted from list
		remove_by_index(I, idx);
		d = I.size();
	    }
	    return x;
	}

	// samples a column vector of length=length(l)=length(u)
	// from the standard multivariate normal distribution,
	// truncated over the region [l,u], where -a<l<u<a for some
	// 'a' and l and u are column vectors;
	// uses acceptance rejection and inverse-transform method;
	template<class RNGType>	
	Vector tn(const Vector& l, const Vector& u, RNGType *rng) {
	    const double tol = 2.05; // controls switch between methods
	    // threshold can be tuned for maximum
	    // speed for each platform
	    Vector x(l.size());
	    std::vector<int> I; I.reserve(l.size());
	    std::vector<int> notI; notI.reserve(l.size());
	    for (int i=0; i<l.size(); i++)
		if (std::abs(u(i)-l(i))>tol)
		    I.push_back(i);
		else
		    notI.push_back(i);

	    // case abs(u-l) > tol
	    if (I.size()!=0) {
		Vector tl(I.size());
		Vector tu(I.size());
		for (int i=0;i<I.size();i++) {
		    tl(i) = l(I[i]);
		    tu(i) = u(I[i]);
		}
		auto tmp = trnd(tl,tu,rng);
		for (int i=0;i<I.size();i++) {
		    x(I[i]) = tmp(i);
		}
	    }
	    // case abs(u-l) < tol
	    if (notI.size()!=0) {
		Vector tl(notI.size());
		Vector tu(notI.size());
		for (int i=0;i<notI.size();i++) {
		    tl(i) = l(notI[i]);
		    tu(i) = u(notI[i]);
		}
		Vector pl = pnorm(tl, 0., 1., 1, 0);
		Vector pu = pnorm(tu, 0., 1., 1, 0);
		int tmplenfixme = tl.size();
		Vector tmpfixme = pl.array()
		    + (pu-pl).array() * runif(tmplenfixme, 0., 1.,rng).array();
		auto tmp = qnorm(tmpfixme, 0., 1., 1, 0);
		for (int i=0;i<notI.size();i++) {
		    x(notI[i]) = tmp(i);
		}
	    }
	    return x;
	}

	
	//  truncated normal generator
	// * efficient generator of a vector of length(l)=length(u)
	// from the standard multivariate normal distribution,
	// truncated over the region [l,u];
	// infinite values for 'u' and 'l' are accepted;
	// * Remark:
	// If you wish to simulate a random variable
	// 'Z' from the non-standard Gaussian N(m,s^2)
	// conditional on l<Z<u, then first simulate
	// X=trandn((l-m)/s,(u-m)/s) and set Z=m+s*X;
	//
	// Reference:
	// Z. I. Botev (2015),
	// "The Normal Law Under Linear Restrictions:
	//  Simulation and Estimation via Minimax Tilting", submitted to JRSS(B)
	template <class RNGType>
	Vector trandn(Vector& l, Vector& u, RNGType* rng) {
	    ASSERT_MSG((l.array()<=u.array()).all(),
		       "Truncation limits have to be l<=u");

	    Vector x = Vector::Zero(l.size());
	    // treshold for switching between methods
	    // threshold can be tuned for maximum speed for each Matlab version
	    const double a = .4; 
	    // three cases to consider:
	    // case 1: a<l<u
	    std::vector<int> I; I.reserve(l.size());
	    for (int i=0; i<l.size(); i++)
		if (l(i) > a)
		    I.push_back(i);
	    if (I.size()!=0) {
		Vector tl(I.size());
		Vector tu(I.size());
		for (int i=0;i<I.size();i++) {
		    tl(i) = l(I[i]);
		    tu(i) = u(I[i]);
		}
		auto tmp = ntail(tl,tu,rng);
		for (int i=0;i<I.size();i++) {
		    x[I[i]] = tmp[i];
		}
	    }
	    // case 2: l<u<-a
	    I.clear();
	    for (int i=0; i<l.size(); i++)
		if (u(i) < -a)
		    I.push_back(i);
	    if (I.size()!=0) {
		Vector tl(I.size());
		Vector tu(I.size());
		for (int i=0;i<I.size();i++) {
		    tl(i) = -u(I[i]);
		    tu(i) = -l(I[i]);
		}
		auto tmp = -ntail(tl,tu,rng);
		for (int i=0;i<I.size();i++) {
		    x[I[i]] = tmp[i];
		}
	    }
	    // case 3: otherwise
	    I.clear();
	    for (int i=0; i<l.size(); i++)
		if (!(u(i)<-a || l(i)>a))
		    I.push_back(i);
	    if (I.size()!=0) {
		Vector tl(I.size());
		Vector tu(I.size());
		for (int i=0;i<I.size();i++) {
		    tl(i) = l(I[i]);
		    tu(i) = u(I[i]);
		}
		auto tmp = tn(tl,tu,rng);
		for (int i=0;i<I.size();i++) {
		    x[I[i]] = tmp[i];
		}
	    }
	    return x;
	}
	

      //  Computes permuted lower Cholesky factor L for Sig
      //  by permuting integration limit vectors l and u.
      //  Outputs perm, such that Sig(perm,perm)=L%*%t(L).
      //
      // Reference: 
      //  Gibson G. J., Glasbey C. A., Elston D. A. (1994),
      //  "Monte Carlo evaluation of multivariate normal integrals and
      //  sensitivity to variate ordering", 
      //  In: Advances in Numerical Methods and Applications, pages 120--126
      void cholperm(/*in&output*/ Matrix& Sig, Vector& l, Vector& u,
		    /*output*/ Matrix& L, Vector& perm) {

	const double eps = 1e-10; // round-off error tolerance
	int d = l.size(); Vector z = Vector::Zero(d);
	L.array() = 0.; for (int i=0;i<d;i++) perm(i) = i;
	Vector D = Sig.diagonal();

	std::unique_ptr<double> uniqdata(new double[d*5]);
	double * data = uniqdata.get(); //	    
	double * datapr = data;
	double * datas = datapr + d;
	double * datacols = datas + d;
	double * datatl = datacols + d;
	double * datatu = datatl + d;
   
	for (int j=0; j<d; j++) {
	  Eigen::Map<Vector> pr(datapr, d); pr.array() = INF;
	  Eigen::Map<Vector> s(datas,d-j);
	  Eigen::Map<Vector> cols(datacols,d-j);
	  Eigen::Map<Vector> tl(datatl,d-j);
	  Eigen::Map<Vector> tu(datatu,d-j);

	  if (j > 1) {
	    s = D.segment(j,d-j).array() -
		L.block(j,0,d-j,j).array().square().rowwise().sum();
	    cols = L.block(j,0,d-j,j) * z.head(j);
	  } else if (j==1) {
	    s = D.tail(d-j).array() -
		L.block(j,0,d-j,1).array().square();
	    cols = L.block(j,0,d-j,1)*z(0);
	  } else {
	    s = D.tail(d-j);
	    cols.array() = 0.;
	  }
	
	  for (int i=0; i<s.size(); i++) {
	    if (s(i) < 0)
	      s(i) = std::sqrt(eps);
	    else
	      s(i) = std::sqrt(s(i));
	  }
	
	  tl = (l.tail(d-j)-cols).array() / s.array();
	  tu = (u.tail(d-j)-cols).array() / s.array();
	  pr.tail(d-j) = lnNpr(tl, tu);
	  // find smallest marginal dimension
	  double k; pr.minCoeff(&k);
	  // flip dimensions k-->j
	  Sig.row(j).swap(Sig.row(k));
	  Sig.col(j).swap(Sig.col(k));
	  std::swap(D(j),D(k));
	  L.row(j).swap(L.row(k)); // update only rows of L
	  std::swap(l(j), l(k)); std::swap(u(j), u(k)); // update integration limit
	  std::swap(perm(j), perm(k)); // keep track of permutation
	  // construct L sequentially via Cholesky computation
	  double tmp = Sig(j,j) - L.block(j,0,1,j).array().square().sum();
	  if (tmp < -0.001) {
	    LOG_WARNING("Sigma is not positive semi-definite");
	  }
	
	  if (tmp < 0) tmp = eps;
	  L(j,j) = std::sqrt(tmp);
	  if (j < d-1) {
	    if (j > 1) {
	      L.block(j+1,j,d-j-1,1)
		= (Sig.block(j+1,j,d-j-1,1) -
		   L.block(j+1,0,d-j-1,j)*L.block(j,0,1,j).transpose()).array()
		/ L(j,j);
	    } else if (j==1) {
	      L.block(j+1,j,d-j-1,1)
		= (Sig.block(j+1,j,d-j-1,1) -
		   L.block(j+1,0,d-j-1,1)*L(j,0)).array()
		/ L(j,j);
	    } else if (j==0) {
	      L.block(j+1,j,d-j-1,1) = Sig.block(j+1,j,d-j-1,1).array() /
		L(j,j);
	    }
	  }
	  // find mean value, z(j), of truncated normal:
	  double tmpl =
	    (l(j) - (L.block(j,0,1,j+1).transpose().array() *
		     z.head(j+1).array()).sum()) / L(j,j);
	  double tmpu =
	    (u(j) - (L.block(j,0,1,j+1).transpose().array() *
		     z.head(j+1).array()).sum()) / L(j,j);
	
	  double w = lnNpr(tmpl, tmpu); // aids in computing expected value
	  // of trunc. normal
	  z(j) = (std::exp(-.5*tmpl*tmpl-w)-std::exp(-.5*tmpu*tmpu-w)) /
	      std::sqrt(2*3.1415926535897932384626433);
	}
      }

      // implements grad_psi(x) to find optimal exponential twisting;
      // assume scaled 'L' with zero diagonal;
      void gradpsi(const Vector& y, const Matrix& L, const Vector& l,
		   const Vector& u, /*output:*/ Vector& grad, Matrix& Jac) {
	int d = u.size();
	Vector c(d);
	Vector x(d);
	Vector mu(d);
	x.head(d-1) = y.head(d-1); x(d-1) = 0.;
	mu.head(d-1) = y.tail(d-1); mu(d-1) = 0.;
	// compute now ~l and ~u
	c.tail(d-1) = L.bottomRightCorner(d-1,d)*x; c(0) = 0.;

	Vector lt = l - mu - c;
	Vector ut = u - mu - c;
	// compute gradients avoiding catastrophic cancellation
	Vector w = lnNpr(lt, ut);
	Vector pl = (-lt.array().square()*0.5-w.array()).exp() /
	  std::sqrt(2*3.1415926535897932384626433);
	Vector pu = (-ut.array().square()*0.5-w.array()).exp() /
	  std::sqrt(2*3.1415926535897932384626433);
	Vector P = pl-pu;

	// output the gradient
	Vector dfdx = -mu.head(d-1) + (P.transpose()*L.topLeftCorner(d,d-1)).transpose();
	Vector dfdm = mu - x + P;
	grad.head(d-1) = dfdx;
	grad.tail(d-1) = dfdm.head(d-1);
	// here compute Jacobian matrix
	for (int i=0; i<d; i++) {
	  if (std::isinf(lt(i)))
	    lt(i) = 0.;
	  if (std::isinf(ut(i)))
	    ut(i) = 0.;
	}

	Vector dP = -P.array().square()+lt.array()*pl.array()-ut.array()*pu.array();
	Matrix DL = L.array().colwise()*dP.array();
	Matrix mx = DL; mx.diagonal().array() -= 1;
	Matrix xx = L.transpose()*DL;

	Jac.topLeftCorner(d-1,d-1) = xx.topLeftCorner(d-1,d-1);
	Jac.bottomLeftCorner(d-1,d-1) = mx.topLeftCorner(d-1,d-1);
	Jac.topRightCorner(d-1,d-1) = mx.topLeftCorner(d-1,d-1).transpose();	
	Jac.bottomRightCorner(d-1,d-1) =
	  (dP.head(d-1).array() + 1.).matrix().asDiagonal();
      }

      Vector nleq(const Vector &l, const Vector &u, const Matrix &L) {
	const int d = l.size();
	Vector x = Vector::Zero(2*d-2); // initial point for Newton iteration
	double err = INF; int iter = 0;
	Vector grad(2*d-2);
	Matrix Jac(2*d-2, 2*d-2);
	// Newton correction
	while(err>1e-10) { 
	  gradpsi(x,L,l,u,grad,Jac);
	  x -= Jac.lu().solve(grad);  // FIXME: Jac is symmetry, can it be faster?
	                              // Ref: Bunch, James R., and Beresford N. Parlett.
                        	      // "Direct methods for solving symmetric indefinite
	                              // systems of linear equations." SIAM Journal on
	                              // Numerical Analysis 8.4 (1971): 639-655.
	  err = grad.array().square().sum();
	  iter++;
	  if (iter>100) {
	    LOG_WARNING("Covariance matrix is ill-conditioned and method failed, return rubbish.");
	    break;
	  }
	}
	return x;
      }

      // implements psi(x,mu); assume scaled 'L' without diagonal;
	inline double psy(const Vector& in_x, const Matrix& L, const Vector& l,
			  const Vector& u, const Vector& in_mu) {
	    const int d = u.size();
	    Vector x(d);
	    x.head(d-1) = in_x; x(d-1) = 0.;
	    Vector mu(d);
	    mu.head(d-1) = in_mu; mu(d-1) = 0.;

	    // compute now ~l and ~u
	    Vector c = L*x;
	    Vector tmp1 = l-mu-c;
	    Vector tmp2 = u-mu-c;
	    return (lnNpr(tmp1, tmp2).array() +
		    .5*mu.array().square() - x.array()*mu.array()).sum();
	}

      // generates the proposals from the exponentially tilted 
      // sequential importance sampling pdf;
      // output:    'logpr', log-likelihood of sample
      //             Z, random sample 
	template <class RNGType>
	void mvnrnd(const int& n, Matrix& L, Vector& l, Vector& u, const Vector& in_mu,
		    /*output:*/ Vector &p, Matrix& Z, RNGType * rng) {
	    int d = l.size();
	    Vector mu(d);
	    mu.segment(0,d-1) = in_mu; mu(d-1) = 0;
	    Z = Matrix::Zero(d,n);
	    p.array() = 0;
	    for (int k=0; k<d; k++) {
		// compute matrix multiplication L*Z
		Vector col = Z.topLeftCorner(k+1,n).transpose() *
		    L.block(k,0,1,k+1).transpose();
		// compute limits of truncation
		Vector tl  = -col.array() + (l(k) - mu(k));
		Vector tu  = -col.array() + (u(k) - mu(k));
		// simulate N(mu,1) conditional on [tl,tu]
		Z.row(k) = mu(k) + trandn(tl,tu,rng).array();
		// update likelihood ratio
		p.array() += lnNpr(tl,tu).array() +
		    .5*mu(k)*mu(k) - mu(k)*Z.row(k).transpose().array();
	    }
	}
      
      //// truncated multivariate normal generator
      // simulates 'n' random vectors exactly/perfectly distributed
      // from the d-dimensional N(0,Sig) distribution (zero-mean normal
      // with covariance 'Sig') conditional on l<X<u;
      // infinite values for 'l' and 'u' are accepted;
      // output:   'd' times 'n' array 'rv' storing random vectors;
      //
      // * Example:
      //  d=60;n=10^3;Sig=0.9*matrix(1,d,d)+0.1*diag(d);l=(1:d)/d*4;u=l+2;
      //  X=mvrandn(l,u,Sig,n);boxplot(t(X)) // plot marginals
      //
      // * Notes: Algorithm may not work if 'Sig' is close to being rank deficient.
      // Reference:
      // Z. I. Botev (2015), "The Normal Law Under Linear Restrictions:
      //  Simulation and Estimation via Minimax Tilting", submitted to JRSS(B)
	template <class RNGType>
	Matrix rtmvnorm(const int& n, Vector l, Vector u, Matrix Sig,
			RNGType* rng) {
	    int d = l.size();
	    ASSERT_MSG(u.size()==d && d==Sig.rows() && d==Sig.cols()
		       && (l.array()<=u.array()).all(), "invalid inputs");
	    // Cholesky decomposition of matrix
	    Matrix Lfull(d,d);
	    Vector perm(d);
	    cholperm(Sig, l, u, Lfull, perm);
	    Vector D = Lfull.diagonal();
	    if ((D.array()<1e-10).any()) {
		LOG_WARNING("Method may fail as covariance matrix is singular!");
	    }
	    // rescale
	    Matrix L = Lfull.array().colwise() / D.array();
	    u.array() /= D.array();
	    l.array() /= D.array();
	    // remove diagonal
	    L.diagonal().array() = 0.;
	    // find optimal tilting parameter via non-linear equation solver
	    Vector xmu = nleq(l,u,L);
	    //assign saddlepoint x* and mu*
	    Vector x = xmu.head(d-1);
	    Vector mu = xmu.tail(d-1);
	    // compute psi star
	    double psistar = psy(x,L,l,u,mu);
	    // start acceptance rejection sampling
	    int iter=0; Matrix rv(d,n);
	    Vector logpr(n);
	    Matrix Z(d,n);
	    int accepted = 0;
	    while(true) {
		mvnrnd(n,L,l,u,mu,logpr,Z,rng);
		Vector unsam = runif(n, 0., 1., rng);
		for (int i=0; i<n; i++) {
		    if (accepted >=n )
			break;
		    if (-std::log(unsam(i))> psistar-logpr(i)) {
			// accepted
			rv.col(accepted) = Z.col(i);
			accepted++;
		    }
		}
		iter++;
		if (iter==1e3) {
		    LOG_WARNING("Acceptance prob. smaller than 0.001");
		}  else if (iter>1e4) {
		    for (int i=0; i<n-accepted; i++) {
			rv.col(accepted+i) = Z.col(i);
		    }
		    accepted = n;
		}
		if (accepted >= n) break;
	    }

	    // finish sampling; postprocessing
	    rv = Lfull*rv; // reverse scaling of L
	    VectorTemplate<int> index(d);
	    std::iota(index.data(), index.data()+d, 0);
	    std::sort(index.data(), index.data()+d, [&perm](const int& i1, const int& i2)
		      { return perm(i1) < perm(i2); } );
	    return Eigen::PermutationWrapper<VectorTemplate<int>>(index).transpose()*rv;
	}

    } // namespace truncatednormal    
}  // namespace bnc

#endif /* TRUNCATEDNORMAL_H */
