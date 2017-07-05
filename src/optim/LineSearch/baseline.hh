#ifndef MORETHUENTE_H_
#define MORETHUENTE_H_

#include <cmath>

#include <matrix/matrix.hh>

namespace bnc {

    namespace baseline {
    
	class MoreThuente {
	public:
	    /*
	     * Line Search baseline, modified from project: https://github.com/PatWie/CppNumericalSolvers/tree/master/include/cppoptlib
	     */
	    template<typename F, typename G>
	    static double search(F fFunc, G gFunc, const Vector &x, const Vector &searchDir, const double& l, const double& u, const double& tol=1e-10) {
		return linesearch(fFunc, gFunc, x, searchDir);
	    }
	    
	    template<typename F, typename G>
	    static double linesearch(F fFunc, G gFunc, const Vector &x, const Vector &searchDir, const  double alpha_init = 1.0) {
		// assume step width
		double ak = alpha_init;

		double fval = fFunc(x);
		Vector g = gFunc(x);

		Vector s = searchDir;
		Vector xx = x;

		cvsrch(fFunc, gFunc, xx, fval, g, ak, s);

		return ak;
	    }

	    template<typename F, typename G>    
	    static int cvsrch(F fFunc, G gFunc, Vector &x, double f, Vector &g, double &stp, Vector &s) {
		// we rewrite this from MIN-LAPACK and some MATLAB code
		int info           = 0;
		int infoc          = 1;
		const double xtol   = 1e-15;
		const double ftol   = 1e-4;
		const double gtol   = 1e-2;
		const double stpmin = 1e-15;
		const double stpmax = 1e15;
		const double xtrapf = 4;
		const int maxfev   = 20;
		int nfev           = 0;

		double dginit = g.dot(s);
		if (dginit >= 0.0) {
		    // no descent direction
		    // TODO: handle this case
		    return -1;
		}

		bool brackt      = false;
		bool stage1      = true;

		double finit      = f;
		double dgtest     = ftol * dginit;
		double width      = stpmax - stpmin;
		double width1     = 2 * width;
		Vector wa = x.eval();

		double stx        = 0.0;
		double fx         = finit;
		double dgx        = dginit;
		double sty        = 0.0;
		double fy         = finit;
		double dgy        = dginit;

		double stmin;
		double stmax;

		while (true) {
		    // make sure we stay in the interval when setting min/max-step-width
		    if (brackt) {
			stmin = std::min(stx, sty);
			stmax = std::max(stx, sty);
		    } else {
			stmin = stx;
			stmax = stp + xtrapf * (stp - stx);
		    }

		    // Force the step to be within the bounds stpmax and stpmin.
		    stp = std::max(stp, stpmin);
		    stp = std::min(stp, stpmax);

		    // Oops, let us return the last reliable values
		    if (
			(brackt && ((stp <= stmin) || (stp >= stmax)))
			|| (nfev >= maxfev - 1 ) || (infoc == 0)
			|| (brackt && ((stmax - stmin) <= (xtol * stmax)))) {
			stp = stx;
		    }

		    // test new point
		    x = wa + stp * s;
		    f = fFunc(x);
		    g = gFunc(x);
		    nfev++;
		    double dg = g.dot(s);
		    double ftest1 = finit + stp * dgtest;

		    // all possible convergence tests
		    if ((brackt & ((stp <= stmin) | (stp >= stmax))) | (infoc == 0))
			info = 6;

		    if ((stp == stpmax) & (f <= ftest1) & (dg <= dgtest))
			info = 5;

		    if ((stp == stpmin) & ((f > ftest1) | (dg >= dgtest)))
			info = 4;

		    if (nfev >= maxfev)
			info = 3;

		    if (brackt & (stmax - stmin <= xtol * stmax))
			info = 2;

		    if ((f <= ftest1) & (fabs(dg) <= gtol * (-dginit)))
			info = 1;

		    // terminate when convergence reached
		    if (info != 0)
			return -1;

		    if (stage1 & (f <= ftest1) & (dg >= std::min(ftol, gtol)*dginit))
			stage1 = false;

		    if (stage1 & (f <= fx) & (f > ftest1)) {
			double fm = f - stp * dgtest;
			double fxm = fx - stx * dgtest;
			double fym = fy - sty * dgtest;
			double dgm = dg - dgtest;
			double dgxm = dgx - dgtest;
			double dgym = dgy - dgtest;

			cstep( stx, fxm, dgxm, sty, fym, dgym, stp, fm, dgm, brackt, stmin, stmax, infoc);

			fx = fxm + stx * dgtest;
			fy = fym + sty * dgtest;
			dgx = dgxm + dgtest;
			dgy = dgym + dgtest;
		    } else {
			// this is ugly and some variables should be moved to the class scope
			cstep( stx, fx, dgx, sty, fy, dgy, stp, f, dg, brackt, stmin, stmax, infoc);
		    }

		    if (brackt) {
			if (fabs(sty - stx) >= 0.66 * width1)
			    stp = stx + 0.5 * (sty - stx);
			width1 = width;
			width = fabs(sty - stx);
		    }
		}
		return 0;
	    }

	    static int cstep(double& stx, double& fx, double& dx, double& sty, double& fy, double& dy, double& stp,
			     double& fp, double& dp, bool& brackt, double& stpmin, double& stpmax, int& info) {
		info = 0;
		bool bound = false;

		// Check the input parameters for errors.
		if ((brackt & ((stp <= std::min(stx, sty) ) | (stp >= std::max(stx, sty)))) | (dx * (stp - stx) >= 0.0)
		    | (stpmax < stpmin)) {
		    return -1;
		}

		double sgnd = dp * (dx / fabs(dx));

		double stpf = 0;
		double stpc = 0;
		double stpq = 0;

		if (fp > fx) {
		    info = 1;
		    bound = true;
		    double theta = 3. * (fx - fp) / (stp - stx) + dx + dp;
		    double s = std::max(theta, std::max(dx, dp));
		    double gamma = s * sqrt((theta / s) * (theta / s) - (dx / s) * (dp / s));
		    if (stp < stx)
			gamma = -gamma;
		    double p = (gamma - dx) + theta;
		    double q = ((gamma - dx) + gamma) + dp;
		    double r = p / q;
		    stpc = stx + r * (stp - stx);
		    stpq = stx + ((dx / ((fx - fp) / (stp - stx) + dx)) / 2.) * (stp - stx);
		    if (fabs(stpc - stx) < fabs(stpq - stx))
			stpf = stpc;
		    else
			stpf = stpc + (stpq - stpc) / 2;
		    brackt = true;
		} else if (sgnd < 0.0) {
		    info = 2;
		    bound = false;
		    double theta = 3 * (fx - fp) / (stp - stx) + dx + dp;
		    double s = std::max(theta, std::max(dx, dp));
		    double gamma = s * sqrt((theta / s) * (theta / s)  - (dx / s) * (dp / s));
		    if (stp > stx)
			gamma = -gamma;

		    double p = (gamma - dp) + theta;
		    double q = ((gamma - dp) + gamma) + dx;
		    double r = p / q;
		    stpc = stp + r * (stx - stp);
		    stpq = stp + (dp / (dp - dx)) * (stx - stp);
		    if (fabs(stpc - stp) > fabs(stpq - stp))
			stpf = stpc;
		    else
			stpf = stpq;
		    brackt = true;
		} else if (fabs(dp) < fabs(dx)) {
		    info = 3;
		    bound = 1;
		    double theta = 3 * (fx - fp) / (stp - stx) + dx + dp;
		    double s = std::max(theta, std::max( dx, dp));
		    double gamma = s * sqrt(std::max(static_cast<double>(0.), (theta / s) * (theta / s) - (dx / s) * (dp / s)));
		    if (stp > stx)
			gamma = -gamma;
		    double p = (gamma - dp) + theta;
		    double q = (gamma + (dx - dp)) + gamma;
		    double r = p / q;
		    if ((r < 0.0) & (gamma != 0.0)) {
			stpc = stp + r * (stx - stp);
		    } else if (stp > stx) {
			stpc = stpmax;
		    } else {
			stpc = stpmin;
		    }
		    stpq = stp + (dp / (dp - dx)) * (stx - stp);
		    if (brackt) {
			if (fabs(stp - stpc) < fabs(stp - stpq)) {
			    stpf = stpc;
			} else {
			    stpf = stpq;
			}
		    } else {
			if (fabs(stp - stpc) > fabs(stp - stpq)) {
			    stpf = stpc;
			} else {
			    stpf = stpq;
			}

		    }
		} else {
		    info = 4;
		    bound = false;
		    if (brackt) {
			double theta = 3 * (fp - fy) / (sty - stp) + dy + dp;
			double s = std::max(theta, std::max(dy, dp));
			double gamma = s * sqrt((theta / s) * (theta / s) - (dy / s) * (dp / s));
			if (stp > sty)
			    gamma = -gamma;

			double p = (gamma - dp) + theta;
			double q = ((gamma - dp) + gamma) + dy;
			double r = p / q;
			stpc = stp + r * (sty - stp);
			stpf = stpc;
		    } else if (stp > stx)
			stpf = stpmax;
		    else {
			stpf = stpmin;
		    }
		}

		if (fp > fx) {
		    sty = stp;
		    fy = fp;
		    dy = dp;
		} else {
		    if (sgnd < 0.0) {
			sty = stx;
			fy = fx;
			dy = dx;
		    }

		    stx = stp;
		    fx = fp;
		    dx = dp;
		}

		stpf = std::min(stpmax, stpf);
		stpf = std::max(stpmin, stpf);
		stp = stpf;

		if (brackt & bound) {
		    if (sty > stx) {
			stp = std::min(stx + static_cast<double>(0.66) * (sty - stx), stp);
		    } else {
			stp = std::max(stx + static_cast<double>(0.66) * (sty - stx), stp);
		    }
		}

		return 0;

	    }

	};

    }
}

#endif /* MORETHUENTE_H_ */
