#ifndef BRUTIL_H
#define BRUTIL_H

#include <type_traits>
#include <limits>
#include <cmath>
#include <cfloat>

#include <util/constant.hh>
#include <util/logger.hh>
#include <common/math.hh>
#include <dist/brconst.hh>

namespace bnc {
    
#define ML_POSINF	INF
#define ML_NEGINF	-INF
#define ML_NAN		NaN
#define ML_ERR_return_NAN { LOG_WARNING("return NaN"); return ML_NAN; }

#define dup_no_copy(x) x
    
typedef enum { FALSE = 0, TRUE } Rboolean;    

    template <class RNGType>
    inline double unif_rand(RNGType *rng)
    {
	return rng->uniform();
    }

    template <class RNGType>
    inline double exp_rand(RNGType *rng)
    {
	return rng->exponential();
    }

    template <class RNGType>
    inline double norm_rand(RNGType *rng)
    {
	return rng->normal();
    }

    inline bool R_FINITE(const double &x) {
      return std::isfinite(x);
    }

    inline bool ISNAN(const double &x) {
      return std::isnan(x);
    }

    template<class T>
    inline T R_forceint(const T &x) {
	return nearbyint(x);
    }

    static double myfmod(double x1, double x2)
    {
	double q = x1 / x2;
	return x1 - floor(q) * x2;
    }

    double BR_pow(double x, double y) /* = x ^ y */
    {
	if(x == 1. || y == 0.)
	    return(1.);
	if(x == 0.) {
	    if(y > 0.) return(0.);
	    /* y < 0 */return(ML_POSINF);
	}
	if (R_FINITE(x) && R_FINITE(y))
	    return(pow(x,y));
	if (ISNAN(x) || ISNAN(y)) {
	    return(x + y);
	}
	if(!R_FINITE(x)) {
	    if(x > 0)		/* Inf ^ y */
		return((y < 0.)? 0. : ML_POSINF);
	    else {			/* (-Inf) ^ y */
		if(R_FINITE(y) && y == floor(y)) /* (-Inf) ^ n */
		    return((y < 0.) ? 0. : (myfmod(y,2.) ? x  : -x));
	    }
	}
	if(!R_FINITE(y)) {
	    if(x >= 0) {
		if(y > 0)		/* y == +Inf */
		    return((x >= 1)? ML_POSINF : 0.);
		else		/* y == -Inf */
		    return((x < 1) ? ML_POSINF : 0.);
	    }
	}
	return(ML_NAN);		/* all other cases: (-Inf)^{+-Inf,
				   non-int}; (neg)^{+-Inf} */
    }

    double BR_pow_di(double x, int n)
    {
	double pow = 1.0;

	if (ISNAN(x)) return x;
	if (n != 0) {
	    if (!R_FINITE(x)) return BR_pow(x, (double)n);
	    if (n < 0) { n = -n; x = 1/x; }
	    for(;;) {
		if(n & 01) pow *= x;
		if(n >>= 1) x *= x; else break;
	    }
	}
	return pow;
    }

#define R_RFUNC_INTERFACE_2ARG(FUNC)					\
    template<class RNGType>						\
    Vector FUNC(const int&n, const double& a, const double& b,		\
		 RNGType *rng)						\
    {									\
        Vector ret(n);							\
									\
	for (int i=0; i<n; i++) {					\
	    ret(i) = FUNC(a, b, rng);					\
	}								\
	return ret;							\
									\
    }									\
									\
    template<class RNGType, class AType, class BType>			\
    typename std::enable_if<std::is_class<AType>::value&&std::is_class<BType>::value, \
	      Vector>::type						\
    FUNC(const int&n, const AType& a, const BType& b,			\
		RNGType *rng)						\
    {									\
        Vector ret(n);							\
	int ai=0, bi=0;							\
	for (int i=0; i<n; i++) {					\
	    ret(i) = FUNC(a(ai), b(bi),					\
		       rng);						\
	    if (++ai >= a.size()) ai = 0;				\
	    if (++bi >= b.size()) bi = 0;				\
	}								\
	return ret;							\
    }									\
									\
    template<class RNGType, class BType>				\
    typename std::enable_if<std::is_class<BType>::value,Vector>::type	\
    FUNC(const int&n, const double &a, const BType &b,			\
		 RNGType *rng)						\
    {									\
	Vector ret(n);							\
	int bi=0;							\
	for (int i=0; i<n; i++) {					\
	    ret(i) = FUNC(a, b(bi), rng);				\
	    if (++bi >= b.size()) bi = 0;				\
	}								\
	return ret;							\
    }									\
									\
    template<class RNGType, class AType>				\
    typename std::enable_if<std::is_class<AType>::value,Vector>::type	\
    FUNC(const int&n, const AType &a, const double &b,			\
		 RNGType *rng)						\
    {									\
	Vector ret(n);							\
	int ai=0;							\
	for (int i=0; i<n; i++) {					\
	    ret(i) = FUNC(a(ai), b, rng);				\
	    if (++ai >= a.size()) ai = 0;				\
	}								\
	return ret;							\
    }									\


#define R_RFUNC_INTERFACE_1ARG(FUNC)					\
    template<class RNGType>						\
    Vector FUNC(const int&n, const double& a,				\
		 RNGType *rng)						\
    {									\
        Vector ret(n);							\
									\
	for (int i=0; i<n; i++) {					\
	    ret(i) = FUNC(a, rng);					\
	}								\
	return ret;							\
									\
    }									\
									\
    template<class RNGType, class AType>				\
    Vector FUNC(const int&n, const AType& a,				\
		RNGType *rng)						\
    {									\
        Vector ret(n);							\
	int ai=0;							\
	for (int i=0; i<n; i++) {					\
	    ret(i) = FUNC(a(ai),					\
		       rng);						\
	    if (++ai >= a.size()) ai = 0;				\
	}								\
	return ret;							\
    }									\


#define R_RFUNC_INTERFACE_3ARG(FUNC)					\
    template<class RNGType>						\
    Vector FUNC(const int&n, const double& a, const double& b,		\
		const double& c, RNGType *rng)				\
    {									\
        Vector ret(n);							\
									\
	for (int i=0; i<n; i++) {					\
	    ret(i) = FUNC(a, b, c, rng);				\
	}								\
	return ret;							\
									\
    }									\
									\
    template<class RNGType, class AType, class BType, class CType>	\
    Vector FUNC(const int&n, const AType& a, const BType& b,		\
		const CType &c, RNGType *rng)				\
    {									\
        Vector ret(n);							\
	int ai=0, bi=0, ci=0;						\
	for (int i=0; i<n; i++) {					\
	    ret(i) = FUNC(a(ai), b(bi),					\
			  c(ci), rng);					\
	    if (++ai >= a.size()) ai = 0;				\
	    if (++bi >= b.size()) bi = 0;				\
	    if (++ci >= c.size()) ci = 0;				\
	}								\
	return ret;							\
    }									\
									\
    template<class RNGType, class BType, class CType>			\
    Vector FUNC(const int&n, const double& a, const BType& b,		\
		const CType &c, RNGType *rng)				\
    {									\
        Vector ret(n);							\
	int ai=0, bi=0, ci=0;						\
	for (int i=0; i<n; i++) {					\
	    ret(i) = FUNC(a, b(bi),					\
			  c(ci), rng);					\
	    if (++bi >= b.size()) bi = 0;				\
	    if (++ci >= c.size()) ci = 0;				\
	}								\
	return ret;							\
    }									\
									\
    template<class RNGType, class AType, class CType>			\
    Vector FUNC(const int &n, const AType &a, const double &b,		\
		const CType &c, RNGType *rng)				\
    {									\
        Vector ret(n);							\
	int ai=0, ci=0;							\
	for (int i=0; i<n; i++) {					\
	    ret(i) = FUNC(a(ai), b,					\
			  c(ci), rng);					\
	    if (++ai >= a.size()) ai = 0;				\
	    if (++ci >= c.size()) ci = 0;				\
	}								\
	return ret;							\
    }									\
									\
    template<class RNGType, class AType, class BType>			\
    Vector FUNC(const int&n, const AType& a, const BType& b,		\
		const double &c, RNGType *rng)				\
    {									\
        Vector ret(n);							\
	int ai=0, bi=0;							\
	for (int i=0; i<n; i++) {					\
	    ret(i) = FUNC(a(ai), b(bi),					\
			  c, rng);					\
	    if (++ai >= a.size()) ai = 0;				\
	    if (++bi >= b.size()) bi = 0;				\
	}								\
	return ret;							\
    }									\
									\
    template<class RNGType, class AType>				\
    Vector FUNC(const int&n, const AType &a, const double &b,		\
		const double &c, RNGType *rng)				\
    {									\
        Vector ret(n);							\
	int ai=0;							\
	for (int i=0; i<n; i++) {					\
	    ret(i) = FUNC(a(ai), b,					\
			  c, rng);					\
	    if (++ai >= a.size()) ai = 0;				\
	}								\
	return ret;							\
    }									\
									\
    template<class RNGType, class BType>				\
    Vector FUNC(const int&n, const double &a, const BType &b,		\
		const double &c, RNGType *rng)				\
    {									\
        Vector ret(n);							\
	int bi=0;							\
	for (int i=0; i<n; i++) {					\
	    ret(i) = FUNC(a, b(bi),					\
			  c, rng);					\
	    if (++bi >= b.size()) bi = 0;				\
	}								\
	return ret;							\
    }									\
									\
    template<class RNGType, class CType>				\
    Vector FUNC(const int&n, const double& a, const double& b,		\
		const CType &c, RNGType *rng)				\
    {									\
        Vector ret(n);							\
	int ci=0;							\
	for (int i=0; i<n; i++) {					\
	    ret(i) = FUNC(a, b,						\
			  c(ci), rng);					\
	    if (++ci >= c.size()) ci = 0;				\
	}								\
	return ret;							\
    }									\
									\

    // Get first value of a scalar or vector/matrix class safely
    template<typename T>
    typename std::enable_if<std::is_scalar<T>::value,T>::type& first(T &x) {
	return x;
    }
    template<typename T>
    const typename std::enable_if<std::is_class<T>::value,T>::type::Scalar& first(T& x) {
	return x(0);
    }

    // Get nth value of a scalar or vector/matrix class safely
    template<typename T>
    const typename std::enable_if<std::is_scalar<T>::value,T>::type& nth(const T &x,
							       const int &n) {
	return x;
    }
    template<typename T>
    const typename std::enable_if<std::is_class<T>::value,T>::type::Scalar& nth(const T& x,
								const int &n) {
	return x(n);
    }

    // Check if length is beyond array's length
    template<typename T>
    bool isOut(const T &x, const typename std::enable_if<std::is_scalar<T>::value,int>::type &n) {
	return false;
    }
    template<typename T>
    bool isOut(const T& x, const typename std::enable_if<std::is_class<T>::value,int>::type &n) {
	return n>=x.size();
    }

    // Check if a matrix is dynamic
    template<class T>
    struct matrix_is_dynamic {
	static const bool value = false;
    };

    template<class D>
    struct matrix_is_dynamic<Eigen::Matrix<D, Eigen::Dynamic, Eigen::Dynamic>> {
	static const bool value = true;
    };

    template<class D, int M>
    struct matrix_is_dynamic<Eigen::Matrix<D, M, Eigen::Dynamic>> {
	static const bool value = true;
    };

    template<class D, int N>
    struct matrix_is_dynamic<Eigen::Matrix<D, Eigen::Dynamic, N>> {
	static const bool value = true;
    };

    // Duplicate a matrix without copy values
    // TODO: may need to add Tensor support
    template<class D>
    Eigen::Matrix<D,Eigen::Dynamic,Eigen::Dynamic>
    clone_no_copy(const Eigen::Matrix<D,Eigen::Dynamic,Eigen::Dynamic>& x) {
	return Eigen::Matrix<D,Eigen::Dynamic,Eigen::Dynamic>(x.rows(),x.cols());
    }

    template<class D, int M>
    Eigen::Matrix<D,M,Eigen::Dynamic>
    clone_no_copy(const Eigen::Matrix<D,M,Eigen::Dynamic>& x) {
	return Eigen::Matrix<D,M,Eigen::Dynamic>(x.cols());
    }

    template<class D, int M>
    Eigen::Matrix<D,Eigen::Dynamic,M>
    clone_no_copy(const Eigen::Matrix<D,Eigen::Dynamic,M>& x) {
	return Eigen::Matrix<D,Eigen::Dynamic,M>(x.rows());
    }

    
#define R_DFUNC_INTERFACE_3ARG(TF, BRFUNC)	      	        \
    template<class X, class M, class L>				\
    typename std::enable_if<std::is_scalar<X>::value,X>::type		\
    TF(const X &x, const M &mu,					\
       const L &give_log)					\
    {								\
	    double ret;						\
	    ret = BRFUNC(x, first(mu),				\
			 static_cast<int>(first(give_log)));	\
	    return ret;						\
    }								\
								\
    template<class X, class M, class L>				\
    typename std::enable_if<!std::is_scalar<X>::value,X>::type		\
    TF(const X &x, const M &mu,	const L &give_log)		\
    {								\
	auto ret = dup_no_copy(x);				\
	int i1=0,i3=0;						\
	for (int i=0; i<x.size(); i++)				\
	{							\
	    ret(i) = BRFUNC(nth(x,i), nth(mu,i1),		\
			    nth(give_log,i3));			\
	    if (isOut(mu,++i1)) i1=0;				\
	    if (isOut(give_log,i3++)) i3=0;			\
	}							\
	return ret;						\
    }								\

    
#define R_DFUNC_INTERFACE_4ARG(TF, BRFUNC)	      	        \
    template<class X, class M, class S, class L>		\
    typename std::enable_if<std::is_scalar<X>::value,X>::type		\
    TF(const X &x, const M &mu,					\
       const S &sigma, const L &give_log)			\
    {								\
	    double ret;						\
	    ret = BRFUNC(x, first(mu), first(sigma),		\
			 static_cast<int>(first(give_log)));	\
	    return ret;						\
    }								\
								\
    template<class X, class M, class S, class L>		\
    typename std::enable_if<!std::is_scalar<X>::value,X>::type		\
    TF(const X &x, const M &mu,					\
       const S &sigma, const L &give_log)			\
    {								\
	auto ret = dup_no_copy(x);				\
	int i1=0,i2=0,i3=0;					\
	for (int i=0; i<x.size(); i++)				\
	{							\
	    ret(i) = BRFUNC(nth(x,i), nth(mu,i1),		\
			    nth(sigma,i2),			\
			    nth(give_log,i3));			\
	    if (isOut(mu,++i1)) i1=0;				\
	    if (isOut(sigma,++i2)) i2=0;			\
	    if (isOut(give_log,++i3)) i3=0;			\
	}							\
	return ret;						\
    }								\


#define R_DFUNC_INTERFACE_5ARG(TF, BRFUNC)	      	        \
    template<class X, class M, class S, class K, class L>	\
    typename std::enable_if<std::is_scalar<X>::value,X>::type		\
    TF(const X &x, const M &mu,					\
       const S &sigma, const K &k, const L &give_log)		\
    {								\
	    double ret;						\
	    ret = BRFUNC(x, first(mu), first(sigma),		\
			 first(k),				\
			 static_cast<int>(first(give_log)));	\
	    return ret;						\
    }								\
								\
    template<class X, class M, class S, class K, class L>	\
    typename std::enable_if<!std::is_scalar<X>::value,X>::type		\
    TF(const X &x, const M &mu,					\
       const S &sigma, const K& k, const L &give_log)		\
    {								\
	auto ret = dup_no_copy(x);				\
	int i1=0,i2=0,i3=0,i4=0;				\
	for (int i=0; i<x.size(); i++)				\
	{							\
	    ret(i) = BRFUNC(nth(x,i), nth(mu,i1),		\
			    nth(sigma,i2),			\
			    nth(k,i3),				\
			    nth(give_log,i4));			\
	    if (isOut(mu,++i1)) i1=0;				\
	    if (isOut(sigma,++i2)) i2=0;			\
	    if (isOut(k,++i3)) i3=0;	       		        \
            if (isOut(give_log,++i4)) i4=0;			\
	}							\
	return ret;						\
    }								\


    
#define R_PFUNC_INTERFACE_4ARG(TF, BRFUNC)				\
    template<class X, class M, class T, class L>			\
    typename std::enable_if<std::is_scalar<X>::value,X>::type			\
    TF(const X &x, const M &mu,						\
       const T& lower_tail,						\
       const L &give_log)						\
    {									\
	X ret;								\
	ret = BRFUNC(x, first(mu),					\
		     static_cast<int>(first(lower_tail)),		\
		     static_cast<int>(first(give_log)));		\
	return ret;							\
    }									\
									\
    template<class X, class M, class T, class L>			\
    typename std::enable_if<!std::is_scalar<X>::value,X>::type			\
    TF(const X &x, const M &mu,						\
       const T& lower_tail,						\
       const L &give_log)						\
    {									\
	X ret = dup_no_copy(x);						\
	int i1=0,i3=0,i4=0;						\
	for (int i=0; i<x.size(); i++)					\
	{								\
	    ret(i) = BRFUNC(nth(x,i),					\
			    nth(mu,i1),					\
			    static_cast<int>(nth(lower_tail,i3)),	\
			    static_cast<int>(nth(give_log,i4)));	\
	    if (isOut(mu,++i1)) i1=0;					\
	    if (isOut(lower_tail,++i3)) i3=0;				\
	    if (isOut(give_log,++i4)) i4=0;				\
	}								\
	return ret;							\
    }									\


#define R_PFUNC_INTERFACE_5ARG(TF, BRFUNC)				\
    template<class X, class M, class S, class T, class L>		\
    typename std::enable_if<std::is_scalar<X>::value,X>::type			\
    TF(const X &x, const M &mu,						\
       const S &sigma, const T& lower_tail,				\
       const L &give_log)						\
    {									\
	X ret;								\
	ret = BRFUNC(x, first(mu), first(sigma),			\
		     static_cast<int>(first(lower_tail)),		\
		     static_cast<int>(first(give_log)));		\
	return ret;							\
    }									\
									\
    template<class X, class M, class S, class T, class L>		\
    typename std::enable_if<!std::is_scalar<X>::value,X>::type			\
    TF(const X &x, const M &mu,						\
       const S &sigma, const T& lower_tail,				\
       const L &give_log)						\
    {									\
	X ret = dup_no_copy(x);						\
	int i1=0,i2=0,i3=0,i4=0;					\
	for (int i=0; i<x.size(); i++)					\
	{								\
	    ret(i) = BRFUNC(nth(x,i),					\
			    nth(mu,i1),					\
			    nth(sigma,i2),				\
			    static_cast<int>(nth(lower_tail,i3)),	\
			    static_cast<int>(nth(give_log,i4)));	\
	    if (isOut(mu,++i1)) i1=0;					\
	    if (isOut(sigma,++i2)) i2=0;				\
	    if (isOut(lower_tail,++i3)) i3=0;				\
	    if (isOut(give_log,++i4)) i4=0;				\
	}								\
	return ret;							\
    }									\


#define R_PFUNC_INTERFACE_6ARG(TF, BRFUNC)				\
    template<class X, class M, class S, class K, class T, class L>	\
    typename std::enable_if<std::is_scalar<X>::value,X>::type			\
    TF(const X &x, const M &mu,						\
       const S &sigma, const K& k, const T& lower_tail,			\
       const L &give_log)						\
    {									\
	X ret;								\
	ret = BRFUNC(x, first(mu), first(sigma), first(k),		\
		     static_cast<int>(first(lower_tail)),		\
		     static_cast<int>(first(give_log)));		\
	return ret;							\
    }									\
									\
    template<class X, class M, class S, class K, class T, class L>	\
    typename std::enable_if<!std::is_scalar<X>::value,X>::type			\
    TF(const X &x, const M &mu,						\
       const S &sigma, const K &k, const T &lower_tail,			\
       const L &give_log)						\
    {									\
	X ret = dup_no_copy(x);						\
	int i1=0,i2=0,i3=0,i4=0,i5=0;					\
	for (int i=0; i<x.size(); i++)					\
	{								\
	    ret(i) = BRFUNC(nth(x,i),					\
			    nth(mu,i1),					\
			    nth(sigma,i2),				\
                            nth(k,i3),					\
			    static_cast<int>(nth(lower_tail,i4)),	\
			    static_cast<int>(nth(give_log,i5)));	\
	    if (isOut(mu,++i1)) i1=0;					\
	    if (isOut(sigma,++i2)) i2=0;				\
	    if (isOut(k,++i3)) i3=0;					\
	    if (isOut(lower_tail,++i4)) i4=0;				\
	    if (isOut(give_log,++i5)) i5=0;				\
	}								\
	return ret;							\
    }									\



#define R_QFUNC_INTERFACE_4ARG(TF, BRFUNC)				\
    template<class X, class M, class T, class L>			\
    typename std::enable_if<std::is_scalar<X>::value,X>::type			\
    TF(const X &x, const M &mu,						\
       const T& lower_tail,						\
       const L &give_log)						\
    {									\
	X ret;								\
	ret = BRFUNC(x, first(mu),					\
		     static_cast<int>(first(lower_tail)),		\
		     static_cast<int>(first(give_log)));		\
	return ret;							\
    }									\
									\
    template<class X, class M, class T, class L>			\
    typename std::enable_if<!std::is_scalar<X>::value,X>::type			\
    TF(const X &x, const M &mu,						\
       const T& lower_tail,						\
       const L &give_log)						\
    {									\
	X ret = dup_no_copy(x);						\
	int i1=0,i3=0,i4=0;						\
	for (int i=0; i<x.size(); i++)					\
	{								\
	    ret(i) = BRFUNC(nth(x,i),					\
			    nth(mu,i1),					\
			    static_cast<int>(nth(lower_tail,i3)),	\
			    static_cast<int>(nth(give_log,i4)));	\
	    if (isOut(mu,++i1)) i1=0;					\
	    if (isOut(lower_tail,++i3)) i3=0;				\
	    if (isOut(give_log,++i4)) i4=0;				\
	}								\
	return ret;							\
    }									\


#define R_QFUNC_INTERFACE_5ARG(TF, BRFUNC)				\
    template<class X, class M, class S, class T, class L>		\
    typename std::enable_if<std::is_scalar<X>::value,X>::type			\
    TF(const X &x, const M &mu,						\
       const S &sigma, const T& lower_tail,				\
       const L &give_log)						\
    {									\
	X ret;								\
	ret = BRFUNC(x, first(mu), first(sigma),			\
		     static_cast<int>(first(lower_tail)),		\
		     static_cast<int>(first(give_log)));		\
	return ret;							\
    }									\
									\
    template<class X, class M, class S, class T, class L>		\
    typename std::enable_if<!std::is_scalar<X>::value,X>::type			\
    TF(const X &x, const M &mu,						\
       const S &sigma, const T& lower_tail,				\
       const L &give_log)						\
    {									\
	X ret = dup_no_copy(x);						\
	int i1=0,i2=0,i3=0,i4=0;					\
	for (int i=0; i<x.size(); i++)					\
	{								\
	    ret(i) = BRFUNC(nth(x,i),					\
			    nth(mu,i1),					\
			    nth(sigma,i2),				\
			    static_cast<int>(nth(lower_tail,i3)),	\
			    static_cast<int>(nth(give_log,i4)));	\
	    if (isOut(mu,++i1)) i1=0;					\
	    if (isOut(sigma,++i2)) i2=0;				\
	    if (isOut(lower_tail,++i3)) i3=0;				\
	    if (isOut(give_log,++i4)) i4=0;				\
	}								\
	return ret;							\
    }									\


#define R_QFUNC_INTERFACE_6ARG(TF, BRFUNC)				\
    template<class X, class M, class S, class T, class K, class L>	\
    typename std::enable_if<std::is_scalar<X>::value,X>::type			\
    TF(const X &x, const M &mu,						\
       const S &sigma, const K &k, const T& lower_tail,			\
       const L &give_log)						\
    {									\
	X ret;								\
	ret = BRFUNC(x, first(mu), first(sigma), first(k),		\
		     static_cast<int>(first(lower_tail)),		\
		     static_cast<int>(first(give_log)));		\
	return ret;							\
    }									\
									\
    template<class X, class M, class S, class K, class T, class L>	\
    typename std::enable_if<!std::is_scalar<X>::value,X>::type			\
    TF(const X &x, const M &mu,						\
       const S &sigma, const K& k, const T& lower_tail,			\
       const L &give_log)						\
    {									\
	X ret = dup_no_copy(x);						\
	int i1=0,i2=0,i3=0,i4=0,i5=0;					\
	for (int i=0; i<x.size(); i++)					\
	{								\
	    ret(i) = BRFUNC(nth(x,i),					\
			    nth(mu,i1),					\
			    nth(sigma,i2),				\
			    nth(k,i3),					\
			    static_cast<int>(nth(lower_tail,i4)),	\
			    static_cast<int>(nth(give_log,i5)));	\
	    if (isOut(mu,++i1)) i1=0;					\
	    if (isOut(sigma,++i2)) i2=0;				\
	    if (isOut(k,++i3)) i3=0;					\
	    if (isOut(lower_tail,++i4)) i4=0;				\
	    if (isOut(give_log,++i5)) i5=0;				\
	}								\
	return ret;							\
    }									\



}  // namespace bnc


#endif /* BRUTIL_H */
