#ifndef LOCALLEVEL_H
#define LOCALLEVEL_H

#include <utility>

#include <util/logger.hh>
#include <matrix/matrix.hh>
#include <dist/norm.hh>

namespace bnc {

    class LocalLevelModel
    {
    public:
	/**
	 * univariate filter
	 * Reference: [1] (formula 2.15)
	 */
	static std::pair<Vector, Vector>
	filter(const Vector& obs,  const double& a1,
	       const double& P1,   const double& eta2,
	       const double& eps2)
	    {
		Vector retM(obs.size());
		Vector retV(obs.size());		
		double vt;
		double Kt;
		double at = a1;
		double att;
		double at1;
		double Ft;
		double Pt = P1;
		double Ptt;
		double Pt1;
		for (int i=0; i<obs.size(); i++) {
		    vt = obs(i) - at;
		    Ft = Pt + eps2;
		    Kt = Pt/Ft;
		    att = at + Kt*vt;
		    Ptt = Pt*(1-Kt);
		    at1 = at+Kt*vt;
		    Pt1 = Pt*(1-Kt) + eta2;

		    retM(i) = att;
		    retV(i) = Ptt;
		    
		    Pt = Pt1;
		    at = at1;
		}
		
		return std::make_pair(retM, retV);
	    }

	static std::pair<Vector, Vector>
	smooth(const Vector& obs,  const double& a1,
	       const double& P1,   const double& eta2,
	       const double& eps2)
	    {
		// first do forward filtering and store
		// a_t, P_t, F_t, Lt (=1-Kt)
		Vector va(obs.size()), vP(obs.size()),
		    vF(obs.size()), vL(obs.size()),
		    vv(obs.size());

		double vt;
		double Kt;
		double at = a1;
		double att;
		double at1;
		double Ft;
		double Pt = P1;
		double Ptt;
		double Pt1;
		for (int i=0; i<obs.size(); i++) {
		    vt = obs(i) - at;
		    Ft = Pt + eps2;
		    Kt = Pt/Ft;
		    att = at + Kt*vt;
		    Ptt = Pt*(1-Kt);
		    at1 = at+Kt*vt;
		    Pt1 = Pt*(1-Kt) + eta2;

		    va(i) = at;
		    vP(i) = Pt;
		    vF(i) = Ft;
		    vL(i) = 1-Kt;
		    vv(i) = vt;
		    
		    Pt = Pt1;
		    at = at1;
		}

		// do the backward smoothing
		Vector retM(obs.size());
		Vector retV(obs.size());
		double rt = 0., Nt = 0.;
		for (int i=obs.size()-1; i>=0; i--) {
		    // for mean
		    rt = vv(i)/vF(i) + vL(i)*rt;
		    retM(i) = va(i) + vP(i)*rt;
		    // for variance
		    Nt = 1/vF(i) + pow(vL(i),2)*Nt;
		    retV(i) = vP(i)*(1 - vP(i)*Nt);
		}

		return std::make_pair(retM, retV);
	    }

	/*
	 *  Sampling functions
	 */

	// sample one posterior state series
	template<class RNGType>
	static Vector sample_state(const Vector& obs, const double& a1,
				   const double& P1, const double& eta2,
				   const double& eps2, RNGType *rng)
	    {
		auto postparam = smooth(obs, a1,
					P1, eta2, eps2);
		return rnorm(obs.size(), postparam.first,
			     postparam.second, rng);
	    }

	// sample multiple posterior state series
	template<class RNGType>
	static Matrix sample_state(const int &n, const Vector& obs,
		      const double& a1, const double& P1,
		      const double& eta2, const double& eps2,
		      RNGType *rng)
	    {
		Matrix ret(n, obs.size());
		auto postparam = smooth(obs, a1,
					P1, eta2, eps2);
		for (int i=0; i<n; i++)
		{
		    ret.row(i) = rnorm(obs.size(), postparam.first,
				       postparam.second.array().sqrt(), rng)
			.transpose();
		}
		return ret;
	    }


	template<class RNGType>
	static Matrix sample_obs(const int &n, const Vector& obs,
		      const double& a1, const double& P1,
		      const double& eta2, const double& eps2,
		      RNGType *rng)
	    {
		Matrix ret(n, obs.size());
		auto postparam = smooth(obs, a1,
					P1, eta2, eps2);
		for (int i=0; i<n; i++)
		{
		    ret.row(i) = rnorm(obs.size(), postparam.first,
				       (postparam.second.array() + eps2).sqrt(),
				       rng)
			.transpose();
		}
		return ret;
	    }
	
	class Missing
	{
	public:
	    Missing();
	    virtual ~Missing();
	};


	LocalLevelModel() {};
	~LocalLevelModel() {};
    };
    
}  // namespace bnc



#endif /* LOCALLEVEL_H */

// [1] J.Durbin and S.J. Koopman, Time series analysis by state space methods
