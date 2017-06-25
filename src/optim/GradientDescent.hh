#ifndef GRADIENTDESCENT_H
#define GRADIENTDESCENT_H

#include <cmath>
#include <iostream>
#include <type_traits>

#include <util/logger.hh>
#include <util/constant.hh>
#include <optim/LineSearch/MoreThuente.hh>
#include <optim/optim.hh>
#include <optim/bounded.hh>

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
		    f(res.x);
		    direct = -g(res.x);
		    direct.normalize();
		    double uu  = INF;
		    double tmp = 0.;
		    // Find the common limit (uu)
		    // and handle the active contraints.
		    // The active contraints only
		    // exist and are only useful when uu=0.
		    for (int i=0; i<x0.size(); i++) {
			if (direct(i) > 0) {
			    tmp = (u(i)-res.x(i))/direct(i);
			} else {
			    tmp = (l(i)-res.x(i))/direct(i);
			}
			if (eq(tmp,0,tol)) {
			    // current contraint is active
			    // set this direction to 0.
			    direct(i) = 0.;
			    continue;
			} else {
			    if (uu>tmp) uu = tmp;
			}
		    }
		    // rescale direct


		    uu = std::max(std::min(uu,1e15),0.);

		    // search along direct
		    double step = LS::search(f, g, res.x, direct,
					     0., uu, (uu+0.)*0.5);
		    dx = step*direct;
		    res.x += step*direct;

		    if (std::is_base_of<Bounded, LS>::value) {
			// Safeguard res.x to avoid numerical errors
			// This is important because if the bounds
			// are violated, in the next loop, uu may be 0.
			for (int i=0; i<x0.size(); i++) {
			    if (res.x(i)>u(i)) {
				res.x(i) = u(i);
				continue;
			    }
			    if (res.x(i)<l(i)) {
				res.x(i) = l(i);
				continue;
			    }
			}
		    }
		    // check convergence
		    if (dx.norm() <= tol) {
			g(res.x);
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
