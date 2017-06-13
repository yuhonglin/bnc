/*
 * FIXME: because in many functions, computing f and g are
 * similar, so make sure that f and g are always computed 
 * one after the other. But there is one case: 
 * to test BNC_LS_U1, we need both ft and fl. But the 
 * following code may or may not need gt and gl
 * so one needs to either compute both gt gl before BNC_LS_U1
 * or recompute ft or fl when gt gl is needed. This may 
 * be fixed by recorder U1,U2,U3.
 *
 */


#ifndef MORETHUENTE_H
#define MORETHUENTE_H

#include <iostream>
#include <cmath>

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
 *                  |f'(alpha)| <= eta*|f'(0)|
 *                            l <= alpha <= u
 * with some parameter mu and eta in (0,1).
 * 
 */
	    class MoreThuente
	    {

#define BNC_LS_phi(s)   f(x+direct*s)
#define BNC_LS_dphi(s)  (g(x+direct*s).dot(direct))
#define BNC_LS_U1       (ft > fl)
#define BNC_LS_U2       (ft <= fl && gt*(al-at)>0)
#define BNC_LS_U3       (ft <= fl && gt*(al-at)<0)

	
		static inline double safeguard_U2(const double& at, const double& al,
						  const double& au) {
		    return std::min(at+1.1*(at-al), au);
		}

		static inline double safeguard_refine(const double& at, const double& al,
						      const double& au) {
		    return (al+std::max(7/12.*at, al))/2.;
		}

		/**
		 * safeguard_bisect
		 * 
		 * This function is rather complex but actually its purpose is just to
		 * speed up convergence when I+ is already in [l,u]
		 */

		static double compute_ac (const double& al, const double& at,
					  const double& fl, const double& ft,
					  const double& gl, const double& gt,
					  const double& lpt, const double& lmt) {
		    const double a = (fl-ft + (gl*at-gt*al) - 0.5*(gl-gt)*lpt) /	
			(lmt * (al*(al+4*at) - (0.5*at+3)*at - 1.5)); 
		    const double b = 0.5*(gl-gt)/lmt - 1.5*a*lpt;		
		    const double c = 3*a*al*at - (gl*at-gt*al)/lmt;		
		    return (-2*b + std::sqrt(4*b*b-12*a*c)) / (6*a);
		}

		static double compute_aq (const double& al, const double& at,
					  const double& fl, const double& ft,
					  const double& gl, const double& gt,
					  const double& lpt, const double& lmt) {
		    const double a = ((ft-fl) + gl*lmt) / ((2*at-lpt)/lmt);
		    const double b = gl - 2*a*at;
		    return -b/(2*a);
		}

		static double compute_as (const double& al, const double& at,
					  const double& fl, const double& ft,
					  const double& gl, const double& gt,
					  const double& lpt, const double& lmt) {
		    double a = 0.5*(gl-gt)/lmt;			
		    return al - 0.5*gl/a;
		}

		template <typename F, typename G>
		static inline double safeguard_bisect(const double& al, const double& at,
						      const double& au, const double& fl,
						      const double& ft, const double& gl,
						      const double& gt, F f, G g,
						      const Vector& x,
						      const Vector& direct) {
		    // when this function is called, I+ (\in [l,u]) should have
		    // been already found.
		    // And no worries for this because this function only aims
		    // to speedup the near-solution (local) convergences. When
		    // I+ is not inside [l,u], it usually implies al and au
		    // are far from a*.

		    const double lpt = al+at;
		    const double lmt = al-at;
		    double a, b, c;
		    double ac, aq, as;
		    if (ft > fl) {
			// case 1
			// First, compute ac, the minimiser of cubic
			// interpolation of fl,ft,gl,gt
			// assuming the cubic is ax^3+bx^2+cx+d
			ac = compute_ac(al,at,fl,ft,gl,gt,lpt,lmt);
			// Then, compute aq, the minimiser of quadratic
			// interpolation of fl, ft, gl
			aq = compute_aq(al,at,fl,ft,gl,gt,lpt,lmt);
			// choose
			return (std::abs(ac-al) < std::abs(aq-al)) ?  ac : (0.5*(aq+ac));
			
		    } else if (gt*gl < 0) {
			// First, compute ac, the minimiser of cubic
			// interpolation of fl,ft,gl,gt
			// assuming the cubic is ax^3+bx^2+cx+d
			ac = compute_ac(al,at,fl,ft,gl,gt,lpt,lmt);			
			// Then, compute as, the minimiser of quadratic
			// interpolation of fl, gl and gt
			as = compute_as(al,at,fl,ft,gl,gt,lpt,lmt);
			return (std::abs(ac-at)>=std::abs(as-at)) ? ac : as;
		    } else if (std::abs(gt) <= std::abs(gl)) {
			// First, compute ac, the minimiser of cubic
			// interpolation of fl,ft,gl,gt
			// assuming the cubic is ax^3+bx^2+cx+d
			ac = compute_ac(al,at,fl,ft,gl,gt,lpt,lmt);
			// Then, compute as, the minimiser of quadratic
			// interpolation of fl, gl and gt
			as = compute_as(al,at,fl,ft,gl,gt,lpt,lmt);
			double ret = (std::abs(ac-at)<std::abs(as-at)) ? ac : as;
			return (at>al) ? std::min(at+0.66*(au-at), ret) :
			    std::max(at + 0.66*(au-at), ret);
		    } else {
			// finally, return the minimiser of cubic interpolation
			// of fu, ft, gu, gt
			const double fu = BNC_LS_phi(au);
			const double gu = BNC_LS_dphi(au); // dphi call must be preceded by
			                                   // a phi call with same input
			const double upt = au+at;
			const double umt = au-at;
			double a = (fu-ft + (gu*at-gt*au) - 0.5*(gu-gt)*upt) /
			    (umt * (au*(au+4*at) - (0.5*at+3)*at - 1.5));
			double b = 0.5*(gu-gt)/umt - 1.5*a*upt;
			double c = 3*a*au*at - (gu*at-gt*au)/umt;
			double ac = (-2*b + std::sqrt(4*b*b-12*a*c)) / (6*a);
			// Then, compute aq, the minimiser of quadratic
			// interpouation of fu, ft, gu
			a = ((ft-fu) + gu*umt) / ((2*at-au-at)/umt);
			b = gu - 2*a*at;
			return -b/(2*a);
		    }
		}
		

	    public:
		template<typename T, typename G>
		static double search(T f, G g, const Vector& x, const Vector& direct,
				     const double& l=1e-15, const double& u=1e15,
				     const double& tol=XTOL) {
		    const double f0 = f(x);
		    const double g0 = g(x).dot(direct);
		    if (g0 > 0) {
			// Error: direction is not a decreasing one.
			// For smooth function, we can use -direct
			// but this case happens often due to errors
			// in higher level code, so better return
			// something to notice them.
			return -1.;
		    }

		    double al = 0., au = INF;
		    double at = (l+u)/2.; // any initial value
		    double prevI = au-al;

		    bool   insided = false; // whether [al,au] \in [l,u]

		    double fl, ft, gt, gl;
		    
		    while (true) {
			// test convergence
			prevI = std::abs(au-al);
			if (le(prevI, tol, tol))
			    break;

			// compute ft, fl, gt
			ft = BNC_LS_phi(at);
			gt = BNC_LS_dphi(at); // FIXME: checking BNC_LS_U1 does not need
			                      //        does not need gt
			fl = BNC_LS_phi(al);
			
			if (BNC_LS_U1) {
			    // [al,at] will contain a*
			    au = at;
			    if (insided && le(prevI*0.66,std::abs(au-al),tol)) { 
				// interval decrease is not enough
				// try bisect safeguards
				// only do this when insided,
				// otherwise it the interval may not
				// converge inside [l,u]
				gl = BNC_LS_dphi(al); // last phi call is on al
				at = safeguard_bisect(al, at, au,
						      fl, ft, gl, gt, f, g, x, direct);
				continue;
			    }
			    at = safeguard_refine(at, al, au);
			    if (eq(at, l, tol)) {
				return l;
			    }

			} else if (BNC_LS_U2) {
			    al = at;
			    at = safeguard_U2(at, al, au);
			    if (eq(at, u, tol)) {
				// I+ not found, return upper bound
				// the {at} are increasing
				return u;
			    }
			} else {
			    // [at,al] will contain a*
			    al = at;
			    au = al;
			    if (insided && le(prevI*0.66,std::abs(au-al),tol)) {
				// interval decrease is not enough
				// try bisect safeguards
				// only do this when insided,
				// otherwise it the interval may not
				// converge inside [l,u]
				gl = BNC_LS_dphi(al); // last phi call is on al
				at = safeguard_bisect(al, at, au,
						      fl, ft, gl, gt, f, g, x, direct);
				continue;
			    }
			    at = safeguard_refine(at, al, au);
			    if (eq(at, l, tol)) {
				return l;
			    }
			}
			insided = (au<=u) && (al>=l); // FIXME: I feel this can be
			                               // faster because when one
			                               // pass of case 1,2,3 is processed
			                               // we may be able to know that
			                               // au <= u or al >= l already. So
			                               // maybe no need to check both
			                               // everytime
		    }  // while

		    return al;
		}
		
#undef BNC_LS_phi
#undef BNC_LS_psi
	
	    public:
		MoreThuente();
		~MoreThuente();
	    };
	
	} // namespace linesearch
    } // namespace optim
} // namespace bnc

#endif 
