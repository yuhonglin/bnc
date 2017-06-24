/*
 * Basic facts of this algorithm:
 *
 * 1. In general, the returned value may not satisfy Wolfe
 *    condition.
 * 2. But when the init gradient is decreasing, alpha_min=0
 *    and alpha_max satisfies some condition (easily satisfied
 *    in unconstrained optimisation), the Wolfe condition is
 *    guaranteed to be satisfied by the returned step.
 * 3. When one of the bounds is returned, the Wolfe condition
 *    may NOT be satisfied.
 * 4. If the returned step is not a bound, the Wolfe condition
 *    is satisfied. (it will be a local minimiser of \psi or \phi)
 * 5. Normally, the returned step will satisfy Wolfe condition
 *    with \eta=\mu. But in some cases, the step will satisfy
 *    Wolfe condition with any \eta>0. This is done by let
 *    the \phi'=0.
 * 6. As long as the function is "normal", the algrithm will
 *    terminate in finite steps.
 * 7. Even if there exist local minimisers of \psi or \phi,
 *    the algorithm may not find them (they may not even
 *    satisfy the Wolfe condition in general). Whether one of
 *    them will be found relates to the properties of bounds.
 * 8. In all, for unconstraint optimisation problem and
 *    feasible optimisation method, this result can be guaranteed
 *    to satisfy the sufficient decrease condition (but
 *    not for the curvature condition).
 */

#ifndef MORETHUENTE_H
#define MORETHUENTE_H

#include <iostream>
#include <cmath>
#include <limits>

#include <util/constant.hh>
#include <util/numeric.hh>
#include <util/logger.hh>
#include <matrix/matrix.hh>


namespace bnc { namespace optim { namespace lsrch {
	
/* 
 * Line search method
 * Reference:
 *     Moré, Jorge J., and David J. Thuente. "Line search algorithms 
 *     with guaranteed sufficient decrease." ACM Transactions on 
 *     Mathematical Software (TOMS) 20.3 (1994): 286-307.
 * 
 * For a continuous function f:[0,inf]->R with f'(0)<0, the goal of 
 * this function is to find a alpha>=0 such that
 *
 *                     f(alpha) <= f(0) + mu*f'(0)*alpha
 *                  |f'(alpha)| <= mu*|f'(0)|
 *                            l <= alpha <= u
 * with some parameter mu and eta in (0,1).
 * 
 * when l=0 and u is huge enough (e.g. depends on the lower bound of 
 * f or \psi'(u)>=0), the Wolfe condition is guaranteed. But in 
 * general, even the sufficient decrease condition is actually 
 * not always guaranteed (this can be easily seen by examples).
 * The nonguarantee cases all only happen when l or u is returned.
 * So the algorithm will return a status showing that l or u is 
 * returned.
 *
 */
	    class MoreThuente
	    {
	    public:
		// the interface
		/**
		 * search
		 * 
		 * \param fFunc 
		 * \param gFunc 
		 * \param x
		 * \param direct 
		 * \param stpmin 
		 * \param stpmax 
		 * \return double
		 */
		template<typename T, typename G>
		static double search(T fFunc, G gFunc,
				     const Vector& x, const Vector& direct,
				     const double& stpmin, 
				     const double& stpmax,
				     const double& init = 1.0) {
		    // check inputs
		    if (stpmin > stpmax || stpmin < 0) {
			LOG_WARNING("invalid stpmax or stpmin");
			return -1;
		    }
		    
		    // finit = phi(0)
		    // ginit = phi'(0)
		    const double finit = fFunc(x);
		    const double ginit = gFunc(x).dot(direct);

		    // check inputs
		    if (ginit>=0.) {
			return 0.;
		    }
		    
		    // constants
		    const double xtrapl = 1.1;
		    const double xtrapu = 4.0;
		    const double p5 = .5;
		    const double ftol = 1e-3;
		    const double gtol = .9;
		    const double xtol = XTOL;
		    const double p66 = .66;

		    // // Check if 1.0 satisfies the Wolfe conditions
		    // // This can make it a newton step when it is possible
		    // if (stpmax > 1.0 &&
		    // 	stpmin < 1.0 &&
		    // 	fFunc(x+direct) <= finit+ftol*gFunc(x).dot(direct) &&
		    // 	std::abs(gFunc(x+direct).dot(direct)) <=
		    // 	gtol*std::abs(gFunc(x).dot(direct)))
		    // 	return 1.0;
		    
		    bool brackt = false;
		    bool stage1 = true;
		    
		    double gtest  = ftol*ginit;
		    double width  = stpmax - stpmin;
		    double width1 = width/p5;

		    // The initial value of stp depends on the direct
		    // finding algorithm and the problem. For example,
		    // For gradient descent and quadratic function,
		    // "(stpmax+stpmin)*0.5" is better, but for Newton's
		    // method, 1.0 is better, so it should be set by
		    // the algorithm
		    double stp    = init;
		    
		    auto tmp = x + stp*direct;
		    double f = fFunc(tmp);
		    double g = gFunc(tmp).dot(direct);
		    
		    double stx    = 0.;
		    double fx     = finit;
		    double gx     = ginit;
		    double sty    = 0.;
		    double fy     = finit;
		    double gy     = ginit;
		    double stmin  = 0.;
		    // According to paper, stmax should be initialised
		    // as INF. But according to the lbfgsb.f, it should
		    // be "stp + xtrapu*stp", I think the code is wrong,
		    // or at least not suitable for here. If INF causes
		    // error, try "stpmax*xtrapu" (I think initial stmax
		    // is OK as long as stmax>stpmax).
		    double stmax  = INF; //

		    double fm, fxm, fym, gm, gxm, gym;
		    double ftest;

		    while(true) {
			// If psi(stp) <= 0 and f'(stp) >= 0 for some step
			// then the algorithm enters the second stage.
			// This stage is not necessary but can generate
			// better result
			ftest = finit + stp*gtest;
			if (stage1 & (f <= ftest) & (g >= 0.))
			    stage1 = false;

			// test for warnings
			if (brackt & ((stp <= stmin) | (stp >= stmax))) {
			    LOG_DEBUG("Rounding errors prevent progress");
			    return stp;
			}
			if (brackt & ((stmax-stmin) <= xtol*stmax)) {
			    LOG_DEBUG("xtol test satisfied");
			    return stp;
			}
			if ((stp == stpmax) & (f <= ftest) & (g <= gtest)) {
			    LOG_DEBUG("stp == stpmax");
			    return stpmax;
			}
			if ((stp == stpmin) & (f > ftest | g <= gtest)) {
			    LOG_DEBUG("stp == stpmin");
			    return stpmin;
			}
			if (stp == stx) {
			    LOG_DEBUG("stp == stx");
			    return stp;
			}

			// test for convergence
			if ((f <= ftest) & (std::abs(g) <= gtol*(-ginit))) {
			    LOG_DEBUG("Converged");
			    return stp;
			}
		    
                        // A modified function is used to predict the step during the
                        // first stage if a lower function value has been obtained but 
                        // the decrease is not sufficient.
			// COMMENT: this is where I do not understand, according
			// to the paper (it says, the alg will use stage1 from beginning
			// until stage two then never return back), actually the
			// condition is only stage1==true. And based on limited test,
			// I cannot found big difference between this two. This may be
			// a further modification of the authors in implementation.
			if (stage1 && f <= fx && f > ftest) {
			    // Define the modified function and derivative values.
			    fm  =  f - stp*gtest;
			    fxm = fx - stx*gtest;
			    fym = fy - sty*gtest;
			    gm  =  g - gtest;
			    gxm = gx - gtest;
			    gym = gy - gtest;

			    // Call dcstep to update stx, sty, and to compute the new step.
			    dcstep(stx,fxm,gxm,sty,fym,gym,stp,fm,gm,brackt,stmin,stmax);

			    // Reset the function and derivative values for f.
			    fx = fxm + stx*gtest;
			    fy = fym + sty*gtest;
			    gx = gxm + gtest;
			    gy = gym + gtest;
			} else {
			    dcstep(stx,fx,gx,sty,fy,gy,stp,f,g,brackt,stmin,stmax);
			}

                        // Decide if a bisection step is needed.

			if (brackt) {
			    if (std::abs(sty-stx) >= p66*width1) stp = stx + p5*(sty - stx);
			    width1 = width;
			    width = std::abs(sty-stx);
			}

                        // Set the minimum and maximum steps allowed for stp.
			if (brackt) {
			    stmin = std::min(stx,sty);
			    stmax = std::max(stx,sty);
			} else {
			    stmin = stp + xtrapl*(stp - stx);
			    stmax = stp + xtrapu*(stp - stx);
			}
 
                        // Force the step to be within the bounds stpmax and stpmin.
			stp = std::max(stp,stpmin);
			stp = std::min(stp,stpmax);

                        // If further progress is not possible, let stp be the best
                        // point obtained during the search.
			if ((brackt && ((stp <= stmin) || (stp >= stmax)))
			    || (brackt && ((stmax-stmin) <= (xtol*stmax)))) stp = stx;
			
			// Obtain another function and derivative.
			auto tmp = x + stp*direct;
			f = fFunc(tmp);
			g = gFunc(tmp).dot(direct);
		    }
		}

		static void dcstep(double &stx, double &fx, double &dx, double &sty,
				   double &fy, double &dy, double &stp, double &fp,
				   double &dp, bool &brackt, double &stpmin, double &stpmax) {
		    const double p66 = 0.66;

		    // local variables
		    double gamma,p,q,r,s,sgnd,stpc,stpf,stpq,theta;

		    // FORMER: sgnd = dp*(dx/std::abs(dx));
		    if (dx < 0)
			sgnd = -dp;
		    else
			sgnd = dp;

		    stpf = 0;
		    stpc = 0;
		    stpq = 0;
		    
                    // First case: A higher function value. The minimum is bracketed. 
                    // If the cubic step is closer to stx than the quadratic step, the 
                    // cubic step is taken, otherwise the average of the cubic and 
                    // quadratic steps is taken.

		    if (fp > fx) {
			theta = 3.*(fx - fp)/(stp - stx) + dx + dp;
			s = std::max(std::max(std::abs(theta),std::abs(dx)),std::abs(dp));
			gamma = s*std::sqrt(std::pow(theta/s,2) - (dx/s)*(dp/s));
			if (stp < stx) gamma = -gamma;
			p = (gamma - dx) + theta;
			q = ((gamma - dx) + gamma) + dp;
			r = p/q;
			stpc = stx + r*(stp - stx);
			stpq = stx + ((dx/((fx - fp)/(stp - stx) + dx))/2.)*
	                                                            (stp - stx);
			if (std::abs(stpc-stx) < std::abs(stpq-stx))
			    stpf = stpc;
			else
			    stpf = stpc + (stpq - stpc)/2.;

			brackt = true;
			
                    // Second case: A lower function value and derivatives of opposite 
                    // sign. The minimum is bracketed. If the cubic step is farther from 
                    // stp than the secant step, the cubic step is taken, otherwise the 
                    // secant step is taken.
		    } else if (sgnd < 0.) {
			theta = 3.*(fx - fp)/(stp - stx) + dx + dp;
			s = std::max(std::max(std::abs(theta),std::abs(dx)),std::abs(dp));
			gamma = s*std::sqrt(pow(theta/s,2) - (dx/s)*(dp/s));
			if (stp > stx) gamma = -gamma;
			p = (gamma - dp) + theta;
			q = ((gamma - dp) + gamma) + dx;
			r = p/q;
			stpc = stp + r*(stx - stp);
			stpq = stp + (dp/(dp - dx))*(stx - stp);
			if (std::abs(stpc-stp) > std::abs(stpq-stp)) 
			    stpf = stpc;
			else
			    stpf = stpq;

			brackt = true;

                    // Third case: A lower function value, derivatives of the same sign,
                    // and the magnitude of the derivative decreases.
		    } else if (std::abs(dp) < std::abs(dx)) {

                    // The cubic step is computed only if the cubic tends to infinity 
                    // in the direction of the step or if the minimum of the cubic
                    // is beyond stp. Otherwise the cubic step is defined to be the 
                    // secant step.

			theta = 3.*(fx - fp)/(stp - stx) + dx + dp;
			s = std::max(std::max(std::abs(theta),std::abs(dx)),std::abs(dp));

                    // The case gamma = 0 only arises if the cubic does not tend
                    // to infinity in the direction of the step.

			gamma = s*std::sqrt(std::max(0.,pow(theta/s,2)-(dx/s)*(dp/s)));
			if (stp > stx) gamma = -gamma;
			p = (gamma - dp) + theta;
			q = (gamma + (dx - dp)) + gamma;
			r = p/q;
			if (r < 0. && gamma != 0.)
			    stpc = stp + r*(stx - stp);
			else if (stp > stx)
			    stpc = stpmax;
			else
			    stpc = stpmin;

			stpq = stp + (dp/(dp - dx))*(stx - stp);

			if (brackt) {

                    // A minimizer has been bracketed. If the cubic step is 
                    // closer to stp than the secant step, the cubic step is 
                    // taken, otherwise the secant step is taken.

			    if (std::abs(stpc-stp) < std::abs(stpq-stp))
				stpf = stpc;
			    else
				stpf = stpq;

			    if (stp > stx)
				stpf = std::min(stp+p66*(sty-stp),stpf);
			    else
				stpf = std::max(stp+p66*(sty-stp),stpf);
			} else {
			
                            // A minimizer has not been bracketed. If the cubic step is 
                            // farther from stp than the secant step, the cubic step is 
                            // taken, otherwise the secant step is taken.

			    if (std::abs(stpc-stp) > std::abs(stpq-stp))
				stpf = stpc;
			    else
				stpf = stpq;

			    stpf = std::min(stpmax,stpf);
			    stpf = std::max(stpmin,stpf);
			}

                    // Fourth case: A lower function value, derivatives of the
                    // same sign, and the magnitude of the derivative does not
                    // decrease. If the minimum is not bracketed, the step is either
                    // stpmin or stpmax, otherwise the cubic step is taken.

		    } else {
			if (brackt) {
			    theta = 3.*(fp - fy)/(sty - stp) + dy + dp;
			    s = std::max(std::max(std::abs(theta),std::abs(dy)),std::abs(dp));
			    gamma = s*std::sqrt(pow(theta/s,2) - (dy/s)*(dp/s));
			    if (stp > sty) gamma = -gamma;
			    p = (gamma - dp) + theta;
			    q = ((gamma - dp) + gamma) + dy;
			    r = p/q;
			    stpc = stp + r*(sty - stp);
			    stpf = stpc;
			} else if (stp > stx) 
			    stpf = stpmax;
			else
			    stpf = stpmin;
		    }
		
                    // Update the interval which contains a minimizer.

		    if (fp > fx) {
			sty = stp;
			fy = fp;
			dy = dp;
		    } else {
			if (sgnd < 0.) {
			    sty = stx;
			    fy = fx;
			    dy = dx;
			}
			stx = stp;
			fx = fp;
			dx = dp;
		    }

                    // Compute the new step.

		    stp = stpf;
		    
		}
		
	    };
	
	} // namespace linesearch
    } // namespace optim
} // namespace bnc

#endif 
