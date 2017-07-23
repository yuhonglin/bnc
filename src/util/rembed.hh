/*
 *  REmbed: embed an r in c++ program.
 *  
 *  This class provides "RInside"-like interfaces.
 *    
 *  Motivation: everytime some simulation is done,
 *  we would like to do some post-analysis on the
 *  results (plotting, mcmc diagnostics etc.). This
 *  class enable the program to enter a R repl at 
 *  anytime. And it provide interfaces to define data
 *  to the R global environment before entering repl.
 *
 */

#ifndef REMBED_H
#define REMBED_H

#include <iostream>
#include <vector>
#include <type_traits>
#include <cstdio>

#include <matrix/matrix.hh>
#include <util/logger.hh>

extern "C" {
#include <Rembedded.h>
#include <Rinternals.h>
#include <R_ext/Parse.h>
}

#undef length

namespace bnc {

    class REmbed;

    struct REnvironment {

	std::string name;
	REmbed * ptr_r;
	
	REnvironment(REmbed* r, const std::string& n)
	    : ptr_r(r), name(n) {}
	
	REnvironment(const REnvironment& r) {
	    name = r.name;
	    ptr_r = r.ptr_r;
	}
	
	template<class T> T& operator= (T& t);

    };

  
    class REmbed
    {
    public:
	enum Status {
	    SUCCESS,
	    FILE_NOTFOUND,
	    LIBRARY_NOTFOUND,
	    EVAL_ERR,
	    PARSE_ERR,
	    PARSE_INCOM,
	    INVALID_OUTPUT
	};

    
	REmbed(int argc, char **argv) {
	    Rf_initEmbeddedR(argc, argv);
	};
	virtual ~REmbed() {
	    Rf_endEmbeddedR(0);
	};

	virtual void repl_callback() {
	}
    
	void repl() {
	    R_ReplDLLinit();
	    while(R_ReplDLLdo1() > 0) {
		/* add user actions here if desired */
		repl_callback();
	    }
	}

	// eval and return anwsers
	Status eval(const std::string & line, SEXP & ans) {
	    ParseStatus status;
	    SEXP cmdSexp, cmdexpr = R_NilValue;
	    int errorOccurred;

	    PROTECT(cmdSexp = Rf_allocVector(STRSXP, 1));
	    SET_STRING_ELT(cmdSexp, 0, Rf_mkCharLen(line.c_str(), line.size()));
	    
	    cmdexpr = PROTECT(R_ParseVector(cmdSexp, -1, &status, R_NilValue));

	    switch (status){
	    case PARSE_OK:
		for(int i=0; i<Rf_length(cmdexpr); i++){
		    ans = R_tryEval(VECTOR_ELT(cmdexpr, i),
				    R_GlobalEnv, &errorOccurred);
		    if (errorOccurred) {
			LOG_WARNING("Evaluating R code failed");
			UNPROTECT(2);
			return PARSE_ERR;
		    }
		}
		break;
	    case PARSE_INCOMPLETE:
		// need to read another line
		break;
	    case PARSE_NULL:
		LOG_WARNING("ParseStatus is null");
		UNPROTECT(2);
		return PARSE_ERR;
		break;
	    case PARSE_ERROR:
		LOG_WARNING("Parse Error.");
		UNPROTECT(2);
		return PARSE_ERR;
	    case PARSE_EOF:
		LOG_WARNING("ParseStatus is eof");
		break;
	    default:
		LOG_WARNING("Unknown ParseStatus");
		return PARSE_ERR;
		break;
	    }
	    UNPROTECT(2);
	    return SUCCESS;		
	}
	Status eval(const std::string & line) {
	    SEXP ans;
	    return eval(line, ans);
	}
	Status eval(const std::string & line, double & d) {
	    SEXP ans;
	    Status rc = eval(line, ans);
	    if (rc != SUCCESS)
		return rc;
	    
	    if (isNumber(ans) && !isComplex(ans)) {
		d = asReal(ans);
		return SUCCESS;
	    } else {
		return INVALID_OUTPUT;
	    }
	}
	Status eval(const std::string & line, Vector & v) {
	    SEXP ans;
	    Status rc = eval(line, ans);
	    if (rc != SUCCESS)
		return rc;

	    if (isNumeric(ans) && !isMatrix(ans)) {
		v.resize(XLENGTH(ans));
		for (int i=0; i<XLENGTH(ans); i++) {
		    v(i) = REAL(ans)[i];
		}
		return SUCCESS;
	    } else {
		return INVALID_OUTPUT;
	    }
	}
	Status eval(const std::string & line, Matrix & m) {	    
	    SEXP ans;
	    Status rc = eval(line, ans);
	    if (rc != SUCCESS)
		return rc;

	    if (isNumeric(ans) && isMatrix(ans) && !isComplex(ans)) {
		SEXP dims;
		PROTECT(dims = allocVector(INTSXP, 2));
		dims = getAttrib(ans, R_DimSymbol);
		m.resize(INTEGER(dims)[0], INTEGER(dims)[1]);
		for (int i=0; i<INTEGER(dims)[0]; i++) {
		    for (int j=0; j<INTEGER(dims)[1]; j++) {
			m(i,j) = REAL(ans)[j+i*INTEGER(dims)[1]];
		    }
		}
		UNPROTECT(1);
		return SUCCESS;
	    } else {
		return INVALID_OUTPUT;
	    }
            return SUCCESS;
	}
    
	// Similar to R's library command
	// Input: the library name
	Status library(const std::string& n) {
	    return eval(std::string("library(")
			+ n + std::string(")"));
	}

	// Similar to R's source command
	// Input: the R file path
	Status source(const std::string& path) {
	    return eval(std::string("source(\'")
			+ path + std::string("\')"));
	}

	// Define new variable in embeded R which can be explored
	// in REPL.
	// Input: variable name and its value
	void define_var(const std::string& n, const double& d) {
	    SEXP data, variableName;
	    PROTECT(data = Rf_allocVector(REALSXP, 1));
	    REAL(data)[0] = d;
	    PROTECT(variableName = Rf_install(n.c_str()));
	    Rf_defineVar(variableName, data, R_GlobalEnv);
	    UNPROTECT(2);	    
	}
	void define_var(const std::string& n, const Vector& v) {
	    SEXP data, variableName;
	    PROTECT(data = Rf_allocVector(REALSXP, v.size()));
	    std::memcpy(REAL(data), v.data(), sizeof(double)*v.size());
	    variableName = PROTECT(Rf_install(n.c_str()));
	    Rf_defineVar(variableName, data, R_GlobalEnv);
	    UNPROTECT(2);	    
	}
	void define_var(const std::string& n, const Matrix& m) {
	    SEXP data, variableName;
	    PROTECT(data = Rf_allocMatrix(REALSXP, m.rows(), m.cols()));
	    if (!m.IsRowMajor) {
		std::memcpy(REAL(data), m.data(), sizeof(double)* m.size());
	    } else {
		std::memcpy(REAL(data), m.transpose().data(), sizeof(double)*m.size());
	    }
	    PROTECT(variableName = Rf_install(n.c_str()));
	    Rf_defineVar(variableName, data, R_GlobalEnv);
	    UNPROTECT(2);
	}
	template<class T>
	typename std::enable_if<std::is_arithmetic<T>::value, void>::type
	define_var(const std::string& n, const std::vector<T>& v) {
	    SEXP data, variableName;
	    PROTECT(data = Rf_allocVector(REALSXP, v.size()));
	    for(int i=0; i<v.size(); i++)
		REAL(data)[i] = v[i];
	    variableName = PROTECT(Rf_install(n.c_str()));
	    Rf_defineVar(variableName, data, R_GlobalEnv);
	    UNPROTECT(2);	    
	}

	// bracket sugar
	// Usage
	// Rembed["varname"] = var;
	// but not support var = Rembed["varname"]
	REnvironment operator[] (const std::string& n) {
	    return REnvironment(this, n);
	}
    };

    template<class T>
    T& REnvironment::operator= (T& t) {
	ptr_r->define_var(name, t);
	return t;
    }

    
}  // namespace bnc

#endif /* REMBED_H */
