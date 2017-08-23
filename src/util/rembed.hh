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
 *  Usage:
 *     1. Add compilation include path
 *  [on Mac] for example,
 *     -I/opt/local/Library/Frameworks/R.framework/Headers/
 *     2. Add library path
 *  [on Mac] for example,
 *     -L/opt/local/Library/Frameworks/R.framework/Libraries/
 *     3. Add R library in compilation: -lR
 *     4. Set R_HOME environment variable
 *  [on Mac] for example,
 *     export R_HOME=/opt/local/Library/Frameworks/R.framework/Resources/
 *     5. get instance
 *     REmbed* r = REmbed->get_instance(argc, argv);
 *     
 *  Notice:
 *     One should place the code at the beginning of the main
 *     function because if the first call of get_instance failed,
 *     the process will quit immediately.
 *
 */

#ifndef REMBED_H
#define REMBED_H

#include <iostream>
#include <vector>
#include <map>
#include <type_traits>
#include <vector>
#include <string>
#include <cstdio>

#include <matrix/matrix.hh>
#include <util/logger.hh>

extern "C" {
#include <Rembedded.h>
#include <Rinternals.h>
#include <R_ext/Parse.h>
}

#undef length
#undef ncols
#undef nrows

namespace bnc {

    /* R data types */
    class R_Base
    {
    public:
	R_Base() { type = "unknown"; };
	R_Base(const std::string& s) : type(s) {};
	virtual ~R_Base() {};
	std::string type;
	virtual int size() = 0;
    };

    class R_Double : public R_Base
    {
    public:
	R_Double(const double & v) : R_Base("double"),
				     value(v) {};
	virtual ~R_Double() {};
	R_Double& operator= (const double & v) {
	    value = v;
	    return *this;
	}
	double value;
	virtual int size() { return 1; }
    };

    class R_String : public R_Base
    {
    public:
	R_String(const std::string & v) : R_Base("character"),
					  value(v) {};
	virtual ~R_String() {};
	R_String& operator= (const std::string & v) {
	    value = v;
	    return *this;
	}
	std::string value;
	virtual int size() { return value.size(); }
    };
    
    class R_Vector : public R_Base
    {
    public:
	R_Vector(const Vector& v) : R_Base("vector"),
				    value(v) {};
	virtual ~R_Vector() {};
	R_Vector& operator= (const Vector & v) {
	    value = v;
	    return *this;
	}
	Vector value;
	virtual int size() { return value.size(); }	
    };

    class R_Matrix : public R_Base
    {
    public:
	R_Matrix(const Matrix& v) : R_Base("matrix"),
				    value(v) {};
	virtual ~R_Matrix() {};
	R_Matrix& operator= (const Matrix & v) {
	    value = v;
	    return *this;
	}
	Matrix value;
	std::vector<std::string> rownames;
	std::vector<std::string> colnames;
	virtual int size() { return value.size(); }
    };

    class R_List : public R_Base
    {
    public:
	class R_List_Mapped {
	    R_List* rlist;
	    std::string name;
	    int idx;

	public:
	    R_List_Mapped(R_List* l, const std::string& n="", const int& i=0)
		: name(n), idx(i), rlist(l) {};
	    
	    template<class T>
	    typename std::enable_if<std::is_base_of<R_Base, T>::value, T>::type
	    operator= (const T& value) {
		if (name=="") {
		    rlist->intidx_value[idx] =
			std::make_shared<T>(value);
		} else {
		    rlist->stridx_value[name] =
			std::make_shared<T>(value);
		}
		return value;		
	    };

	    Vector operator= (const Vector& value) {
		if (name=="") {
		    rlist->intidx_value[idx] =
			std::make_shared<R_Vector>(value);
		} else {
		    rlist->stridx_value[name] =
			std::make_shared<R_Vector>(R_Vector(value));
		}
		return value;
	    };

	    double operator= (const double& value) {
		if (name=="") {
		    rlist->intidx_value[idx] =
			std::make_shared<R_Double>(R_Double(value));
		} else {
		    rlist->stridx_value[name] =
			std::make_shared<R_Double>(R_Double(value));
		}
		return value;		
	    };	    

	    Matrix operator= (const Matrix& value) {
		if (name=="") {
		    rlist->intidx_value[idx] =
			std::make_shared<R_Matrix>(R_Matrix(value));
		} else {
		    rlist->stridx_value[name] =
			std::make_shared<R_Matrix>(R_Matrix(value));
		}
		return value;
	    };
	    
	};
	
	R_List() : R_Base("list") {};
	virtual ~R_List() {};

	std::map<int, std::shared_ptr<R_Base>> intidx_value;
	std::map<std::string, std::shared_ptr<R_Base>> stridx_value;

	R_List_Mapped operator[](int idx) {
	    return R_List_Mapped(this, "", idx);
	}

	R_List_Mapped operator[](const std::string& idx) {
	    return R_List_Mapped(this, idx, 0);
	}
	
	virtual int size() { return stridx_value.size() + intidx_value.size(); }
    };

    
    /*
     * REmbed (including REnvironment) is the 
     * main interface to R process.
     *
     */
    
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

	/// use singleton
	//  because each process can only embed one R instance
	//  since it has lots of global variables to initialise.
    private:
	REmbed(int argc, char **argv) {
	    Rf_initEmbeddedR(argc, argv);
	};
    public:
	// default arguments
	// Usage:
	//   - one can use REmbed::dfltargs.clear() to clear it
	//   - or one can use REmbed::dfltargs.append() to extend it
	static std::vector<std::string> dfltargs;
	
	static REmbed& get_instance(int argc=0, char **argv=nullptr) {
	    int full_argc = argc + dfltargs.size();
	    std::vector<char*> full_argv;
	    for (int i = 0; i < argc; i++) {
		full_argv.push_back(argv[i]);
	    }
	    // must locate before push_back(argv) because
	    // the first argument is the exec's name
	    for (const auto& arg : dfltargs) {
		full_argv.push_back((char*)arg.data());
	    }
	    full_argv.push_back(nullptr);
	    
	    static REmbed r(full_argc, full_argv.data());
	    
	    return r;
	}
	    
    public:
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
	SEXP alloc_var(const double& d) {
	    SEXP data;
	    PROTECT(data = Rf_allocVector(REALSXP, 1));
	    REAL(data)[0] = d;
	    UNPROTECT(1);
	    return data;
	}
	SEXP define_var(const std::string& n, const double& d) {
	    SEXP data, variableName;
	    data = PROTECT(alloc_var(d));
	    PROTECT(variableName = Rf_install(n.c_str()));
	    Rf_defineVar(variableName, data, R_GlobalEnv);
	    UNPROTECT(2);
	    return data;
	}
	
	SEXP alloc_var(const Vector& v) {
	    SEXP data;
	    PROTECT(data = Rf_allocVector(REALSXP, v.size()));
	    std::memcpy(REAL(data), v.data(), sizeof(double)*v.size());
	    UNPROTECT(1);
	    return data;
	}
	
	SEXP define_var(const std::string& n, const Vector& v) {
	    SEXP data, variableName;
	    data = PROTECT(alloc_var(v));
	    variableName = PROTECT(Rf_install(n.c_str()));
	    Rf_defineVar(variableName, data, R_GlobalEnv);
	    UNPROTECT(2);
	    return data;
	}
	
	SEXP alloc_var(const Matrix& m) {
	    SEXP data;
	    PROTECT(data = Rf_allocMatrix(REALSXP, m.rows(), m.cols()));
	    if (!m.IsRowMajor) {
		std::memcpy(REAL(data), m.data(), sizeof(double)* m.size());
	    } else {
		std::memcpy(REAL(data), m.transpose().data(), sizeof(double)*m.size());
	    }
	    UNPROTECT(1);
	    return data;
	}
	
	SEXP define_var(const std::string& n, const Matrix& m) {
	    SEXP data, variableName;
	    data = PROTECT(alloc_var(m));
	    PROTECT(variableName = Rf_install(n.c_str()));
	    Rf_defineVar(variableName, data, R_GlobalEnv);
	    UNPROTECT(2);
	    return data;
	}
	
	// for std::vector<int>, std::vector<double> etc.
	template<class T>
	typename std::enable_if<std::is_arithmetic<T>::value, SEXP>::type
	define_var(const std::string& n, const std::vector<T>& v) {
	    SEXP data, variableName;
	    PROTECT(data = Rf_allocVector(REALSXP, v.size()));
	    for(int i=0; i<v.size(); i++)
		REAL(data)[i] = v[i];
	    variableName = PROTECT(Rf_install(n.c_str()));
	    Rf_defineVar(variableName, data, R_GlobalEnv);
	    UNPROTECT(2);
	    return data;
	}
	
	// define a name list of vectors
	SEXP define_var(const std::string& n, std::map<std::string, Vector> m) {
	    SEXP nms = PROTECT(allocVector(STRSXP, m.size()));
	    SEXP res = PROTECT(allocVector(VECSXP, m.size()));
	    int i = 0;
	    for (auto& kv : m) {
		// set name
		SET_STRING_ELT(nms, i, Rf_mkCharLen(kv.first.c_str(),kv.first.size()));
		// set value
		SEXP elem = PROTECT(allocVector(REALSXP, kv.second.size()));
		std::memcpy(REAL(elem), kv.second.data(),
			    sizeof(double)* kv.second.size());
		SET_VECTOR_ELT(res, i, elem);
		UNPROTECT(1);
		i++;
	    }
	    setAttrib(res, R_NamesSymbol, nms);
	    SEXP variableName = PROTECT(Rf_install(n.c_str()));
	    Rf_defineVar(variableName, res, R_GlobalEnv);
	    UNPROTECT(3);
	    return res;
	}
	
	// define a named list of vectors from map
	template<class T>
	SEXP define_var(const std::string& n, std::map<std::string, T> m) {

	    SEXP res = PROTECT(allocVector(VECSXP, m.size()));
	    SEXP nms = PROTECT(allocVector(STRSXP, m.size()));
	    int i = 0;
	    for (auto& kv : m) {
		// set name
		SET_STRING_ELT(nms, i, Rf_mkCharLen(kv.first.c_str(),kv.first.size()));
		// set value
		SEXP elem = PROTECT(allocVector(REALSXP, kv.second.size()));
		std::memcpy(REAL(elem), kv.second.data(),
			    sizeof(double)* kv.second.size());
		SET_VECTOR_ELT(res, i, elem);
		UNPROTECT(1);
		i++;
	    }
	    // set class name
	    setAttrib(res, R_NamesSymbol, nms);
	    // define the list in global environment
	    SEXP variableName = PROTECT(Rf_install(n.c_str()));
	    Rf_defineVar(variableName, res, R_GlobalEnv);

	    UNPROTECT(3);
	    return res;
	}
	
	// get a coda mcmc list
	template<class T>
	SEXP map_to_coda_mcmc(std::map<std::string, T> m, const int& start=1,
			      int end=0, const int& thin=1) {
	    if (m.size() <= 0) {
		return R_NilValue;
	    }

	    const int nrow = m.begin()->second.size();
	    const int ncol = m.size();
	    
	    SEXP ret = PROTECT(Rf_allocMatrix(REALSXP, nrow, ncol));
	    SEXP dimnames = PROTECT(allocVector(VECSXP, 2));
	    SEXP colnames = PROTECT(allocVector(STRSXP, ncol));
	    int i = 0;
	    for(auto& kv : m) {
		// set name
		SET_STRING_ELT(colnames, i, Rf_mkCharLen(kv.first.c_str(),kv.first.size()));
		// set value
		std::memcpy(REAL(ret)+i*nrow, kv.second.data(),
			    sizeof(double)*kv.second.size());
		i++;
	    }
	    SET_VECTOR_ELT(dimnames, 0, R_NilValue);
	    SET_VECTOR_ELT(dimnames, 1, colnames); 
	    // set column names
	    setAttrib(ret, R_DimNamesSymbol, dimnames);

	    // set class
	    SEXP cls = PROTECT(allocVector(STRSXP, 1));
	    SET_STRING_ELT(cls, 0, Rf_mkChar("mcmc"));
	    setAttrib(ret, R_ClassSymbol, cls);

	    // define the mcpar field
	    SEXP attname = PROTECT(allocVector(STRSXP, 1));
	    SET_STRING_ELT(attname, 0, Rf_mkChar("mcpar"));
	    SEXP mcpar = PROTECT(allocVector(REALSXP, 3));
	    REAL(mcpar)[0] = static_cast<double>(start);
	    REAL(mcpar)[1] = static_cast<double>(end);
	    REAL(mcpar)[2] = static_cast<double>(thin);
	    setAttrib(ret, attname, mcpar);

	    UNPROTECT(6);

	    return ret;
	}
	
	// define a coda mcmc object
	template<class T>
	SEXP define_coda_mcmc(const std::string& n, std::map<std::string, T> m,
			      const int& start=1, int end=0, const int& thin=1) {
	    if (library("coda")!=SUCCESS) {
		LOG_WARNING("Coda package not installed, do nothing");
		return nullptr;
	    }

	    if (m.size() != 0) {
		end = start + (m.begin()->second.size()-1)*thin;
	    }

	    SEXP res = PROTECT(map_to_coda_mcmc(m, start, end, thin));
	    // define the list in global environment
	    SEXP variableName = PROTECT(Rf_install(n.c_str()));
	    Rf_defineVar(variableName, res, R_GlobalEnv);
	    UNPROTECT(2);
	    return res;
	}
	
	// define a coda mcmc.list object
	template<class T>
	SEXP define_coda_mcmc_list(const std::string& n,
				   std::vector<std::map<std::string, T>> ml,
				   const int& start=1,
				   int end=0,
				   const int& thin=1) {
	    if (library("coda")!=SUCCESS) {
		LOG_WARNING("Coda package not installed, do nothing");
		return nullptr;
	    }

	    if (ml.size() <= 0 || ml[0].size() <= 0) {
		LOG_WARNING("Empty input, do nothing");
		return nullptr;
	    }
	    
	    end = start + (ml[0].begin()->second.size()-1)*thin;
	    
	    SEXP res = PROTECT(allocVector(VECSXP, ml.size()));
	    for (int i = 0; i < ml.size(); i++) {
		SEXP mcmc = PROTECT(map_to_coda_mcmc(ml[i], start, end, thin));
		SET_VECTOR_ELT(res, i, mcmc);
		UNPROTECT(1);
	    }

	    SEXP cls = PROTECT(allocVector(STRSXP, 1));
	    SET_STRING_ELT(cls, 0, Rf_mkChar("mcmc.list"));
	    setAttrib(res, R_ClassSymbol, cls);

	    // define the list in global environment
	    SEXP variableName = PROTECT(Rf_install(n.c_str()));
	    
	    Rf_defineVar(variableName, res, R_GlobalEnv);
	    UNPROTECT(3);

	    return res;
	}

	// for R specific types
	template<class T>
	typename std::enable_if<std::is_base_of<R_Base, T>::value, SEXP>::type
	alloc_var(std::shared_ptr<T> d) {
	    if (d == nullptr) return nullptr;

	    if (d->type == "vector")		
		return alloc_var(std::dynamic_pointer_cast<R_Vector>(d)->value);
	    else if (d->type == "matrix")
		return alloc_var(std::dynamic_pointer_cast<R_Matrix>(d)->value);
	    else if (d->type == "double")
		return alloc_var(std::dynamic_pointer_cast<R_Double>(d)->value);
	    else if (d->type == "list") {
		// count the length of the list
		int len = 0;
		for (auto& kv : std::static_pointer_cast<R_List>(d)->intidx_value) {
		    if (kv.first <= 0) {
			LOG_WARNING("invalid integer index, do nothing");
			return nullptr;
		    }
		    if (kv.first > len)
			len = kv.first;
		}
		int max_intidx = len;
		len += std::static_pointer_cast<R_List>(d)->stridx_value.size();
		
		SEXP res = PROTECT(allocVector(VECSXP, len));
		SEXP nms = PROTECT(allocVector(STRSXP, len));
		for (auto& kv : std::static_pointer_cast<R_List>(d)->intidx_value) {
		    SEXP elem = PROTECT(alloc_var(kv.second));
		    SET_VECTOR_ELT(res, kv.first-1, elem);
		    UNPROTECT(1);
		}

		int i = max_intidx;
		for (auto& kv : std::static_pointer_cast<R_List>(d)->stridx_value) {
		    SEXP elem = PROTECT(alloc_var(kv.second));
		    SET_VECTOR_ELT(res, i, elem);
		    SET_STRING_ELT(nms, i, Rf_mkCharLen(kv.first.c_str(),kv.first.size()));
		    i++;		    
		}
		setAttrib(res, R_NamesSymbol, nms);
		
		UNPROTECT(2);
		return res;
	    } else if (d->type == "unknown") {
		LOG_WARNING("unknown type, do nothing");
		return nullptr;
	    }
	    
	    TODO; // not implemented yet
	    return nullptr;
	}
	
	template<class T>
	typename std::enable_if<std::is_base_of<R_Base, T>::value, SEXP>::type
	define_var(const std::string& n, std::shared_ptr<T> d) {
	    SEXP data, variableName;
	    data = PROTECT(alloc_var(d));
	    PROTECT(variableName = Rf_install(n.c_str()));
	    Rf_defineVar(variableName, data, R_GlobalEnv);
	    UNPROTECT(2);
	    return data;
	}

	template<class T>
	typename std::enable_if<std::is_base_of<R_Base, T>::value, SEXP>::type
	define_var(const std::string& n, const T& d) {
	    SEXP data, variableName;
	    data = PROTECT(alloc_var(std::make_shared<T>(d)));
	    PROTECT(variableName = Rf_install(n.c_str()));
	    Rf_defineVar(variableName, data, R_GlobalEnv);
	    UNPROTECT(2);
	    return data;
	}
	
        // bracket sugar
	// Usage
	// Rembed["varname"] = var;
	// but not support var = Rembed["varname"]
	REnvironment operator[] (const std::string& n) {
	    return REnvironment(this, n);
	}
    };

    // default arguments
    std::vector<std::string> REmbed::dfltargs = {
	"-q" // no R header
    };
    
    template<class T>
    T& REnvironment::operator= (T& t) {
	ptr_r->define_var(name, t);
	return t;
    }

    
}  // namespace bnc

#endif /* REMBED_H */
