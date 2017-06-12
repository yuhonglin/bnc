#ifndef GRADIENTDESCENT_H
#define GRADIENTDESCENT_H

#include <cmath>

#include <util/logger.hh>
#include <util/constant.hh>
#include <optim/LineSearch/MoreThuente.hh>
#include <optim/optim.hh>

namespace bnc {
    namespace optim {
	
	template<typename LS>
	class GradientDescentTemplate
	{
	public:
	    template<typename F, typename G>
	    static Result min(F f, G g, const Vector& x0,
			      const Vector& l, const Vector& u,
			      const double& tol=DNTOL,
			      const int& maxIter=1000) {
		ASSERT_MSG( (l.array()<=x0.array()).all() &&
			    (x0.array()<=u.array()).all(),
			    "Infeasible inputs" );
		
		Result res;
		res.x = x0;
		res.nIter = 0;
		res.info = SUCCESS;

		Vector dx(x0.size());
		Vector direct(x0.size());
		while (true) {
		    res.nIter++;
		    // compute direct and bound
		    direct = -g(res.x);
		    double uu  = INF;
		    double tmp = 0.;
		    for (int i=0; i<x0.size(); i++) {
			if (direct(i) > 0) {
			    tmp = (u(i)-res.x(i))/direct(i);
			    if (uu>tmp) uu = tmp;
			} else {
			    tmp = (l(i)-res.x(i))/direct(i);
			    if (uu>tmp) uu = tmp;
			}
		    }
		    uu = std::max(std::min(uu,1e15),0.);
		    // search along direct
		    double step = LS::search(f, g, res.x, direct,
					     0., uu, tol);
			
		    dx = step*direct;
		    res.x += step*direct;
		    // check convergence
		    if (dx.norm() <= tol) {
			res.f = f(res.x);
			break;
		    }
		    // check for maxIter
		    if (res.nIter>=maxIter) {
			res.info = MAXITERREACHED;
			break;
		    }
		}

		return res;
	    }
	};
	
	using GradientDescent = GradientDescentTemplate<optim::lsrch::MoreThuente>;
    }  // namespace Namespace
}  // namespace bnc

#endif /* GRADIENTDESCENT_H */
