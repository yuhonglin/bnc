/**
 * Agorithms adapted from Rstats function
 * Scale matrix must be definite
 * see 
 * 1) https://github.com/wch/r-source/blob/af7f52f70101960861e5d995d3a4bec010bc89e6/src/library/stats/src/rWishart.c
 */

#ifndef WISHART_H
#define WISHART_H

#include <Eigen/Eigenvalues>

#include <util/logger.hh>
#include <util/checker.hh>
#include <common/math.hh>
#include <matrix/matrix.hh>
#include <dist/norm.hh>
#include <dist/mvnorm.hh>
#include <dist/chisq.hh>
#include <dist/mvnorm.hh>

extern "C" {
    double br_lgammafn(double);
    double br_gammafn(double);
}

namespace bnc {

    enum RWISHART_ALG { RWISHART_DEFINITION,   // by definition (slow)
			RWISHART_RCORE         // adapted from R-core
    };

    /**
     * std_rWishart_factor, used by rwishart when 
     * "RWISHART_ALG==RWISHART_RCORE"
     * 
     */
    template<MAT_STRUCTURE ms, class RNGType>
    Matrix std_rWishart_factor(const int& nu, const int& p, RNGType* rng) {
	ASSERT_MSG(!(nu < p || p <= 0), "Invalid inputs");

	Matrix ret(p,p);
	for (int j=0; j<p; j++) {      // jth column
	    ret(j*(p+1)) = sqrt( rchisq(nu-j,rng) );
	    for (int i=0; i<j; i++) {
		if (ms == UPPER_TRIANGLE) {
		    ret(i+j*p) = rng->normal();
		    ret(j+i*p) = 0.;
		} else {
		    // ms=LOWER_TRIANGLE
		    ret(j+i*p) = rng->normal();
		    ret(i+j*p) = 0.;
		}
	    }
	}
	return ret;
    }

    template<RWISHART_ALG alg=RWISHART_RCORE, class RNGType>
    Matrix rwishart(const int& df, const Matrix& scale, RNGType* rng)
    {
	ASSERT_MSG(!(scale.rows()!=scale.cols() || df<1),
		   "Invalid inputs");

	if (alg == RWISHART_DEFINITION) {
	    // by definition
	    Matrix ret = Matrix::Zero(scale.rows(), scale.cols());
	    Vector v;
	    for (int i=0; i<df; i++) {
		v = rmvnorm(Vector::Zero(scale.rows()),scale,rng);
		ret += mdot<UPPER_TRIANGLE>(v.transpose());
	    }
	    tri2sym_inplace<UPPER_TRIANGLE>(ret);
	    return ret;
	}

	if (alg == RWISHART_RCORE) {
	    Matrix U = scale.llt().matrixU();
	    ASSERT_MSG(!U.hasNaN(), "Cholesky error");
	    Matrix tmp = std_rWishart_factor<UPPER_TRIANGLE>(df, scale.rows(), rng);
	    tmp = tmp.triangularView<Eigen::Upper>()*U;
	    return mdot<MAT_FULL>(tmp);
	}

	LOG_ERROR("RWishart algorithm not supported, return NaN.");
	return Matrix::Constant(scale.rows(),scale.cols(),
				std::numeric_limits<double>::quiet_NaN());
    }


    double dwishart(const Matrix& x, const int& df,
		    const Matrix& scale, const SCALE& s=NORMAL)
    {
	if (s == NORMAL)
	{
	    // adapt from MCMCpack's dwish
	    double gammapart = 1.;
	    for (int i=0; i<x.rows(); i++) {
		gammapart *= br_gammafn((df-i)/2.);
	    }
	    double denom = gammapart * pow(2.,df*x.rows()/2.) *
		pow(3.141592653589793238462643383280, x.rows()*(x.rows()-1)/4.);
	    double detS = scale.determinant();
	    double detW = x.determinant();
	    double tracehold = scale.llt().solve(x).trace();
	    double num  = pow(detS, -df/2.)*pow(detW, (df-x.rows()-1)/2.)*
		exp(-.5*tracehold);
	    return num/denom;

	} else
	{
	    // modified from MCMCpack's dwish
	    double gammapart = 0.;
	    for (int i=0; i<x.rows(); i++) {
		gammapart += br_lgammafn((df-i)/2.);
	    }
	    double denom = gammapart + df*x.rows()/2.*log(2.) +
		log(3.141592653589793238462643383280)/4.*x.rows()*(x.rows()-1);
	    double detS = scale.determinant();
	    double detW = x.determinant();
	    double tracehold = scale.llt().solve(x).trace();
	    double num  =  -df/2.*log(detS) + log(detW)*(df-x.rows()-1)/2. - .5*tracehold;
	    return num - denom;
	}
    }
    
}  // namespace bnc

#endif /* WISHART_H */
