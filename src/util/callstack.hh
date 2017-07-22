// modified from the StackTrace project:
// https://sourceforge.net/projects/stacktrace/
// Usage:
//      1. Declare a CallStack variable, e.g., ```cs```
//      2. cout << cs.to_string() << endl;
//
// For example, if we want to show the call stack
// of LOG_WARNING in Logger, just do the following,
// before LOG_WARNING is called,
//
// ```C++
// Logger.showWarningCallback = []() {
//    CallStack cs;
//    std::cout << cs.to_string() << std::endll;
// }
// ```
//

#ifndef BACKTRACE_H
#define BACKTRACE_H

#include <string>
#include <vector>
#include <sstream>

#include <stdio.h>
#include <execinfo.h>
#include <cxxabi.h>
#include <dlfcn.h>
#include <stdlib.h>

#define MAX_CALLSTACK_DEPTH 32

#ifndef __GNUC__
#error "Only support Linux/GCC"
#endif

namespace bnc {
   
    struct StackEntry {
	/** Default constructor that clears all fields. */
	StackEntry () : line(0) {
	}

	std::string file;     ///< filename
	int         line;     ///< line number
	std::string function; ///< name of function or method

	/** Serialize entry into a text string. */
	std::string to_string () const {
	    std::ostringstream os;
	    os << file << " : " << function;
	    return os.str();
	}
    };
    
    class CallStack
    {
    public:
	CallStack(const int& num_discard = 0) {
	    using namespace abi;

	    // retrieve call-stack
	    void * trace[MAX_CALLSTACK_DEPTH];
	    int stack_depth = backtrace(trace, MAX_CALLSTACK_DEPTH);

	    for (int i = num_discard+1; i < stack_depth; i++) {
		Dl_info dlinfo;
		if(!dladdr(trace[i], &dlinfo))
		    break;

		const char * symname = dlinfo.dli_sname;

		int    status;
		char * demangled = abi::__cxa_demangle(symname, NULL, 0, &status);
		if(status == 0 && demangled)
		    symname = demangled;

		//printf("entry: %s, %s\n", dlinfo.dli_fname, symname);

		// store entry to stack
		if (dlinfo.dli_fname && symname) {
		    StackEntry e;
		    e.file     = dlinfo.dli_fname;
		    e.line     = 0; // unsupported
		    e.function = symname;
		    stack.push_back(e);
		} else {
		    break; // skip last entries below main
		}

		if (demangled)
		    free(demangled);
	    }	    
	}
	~CallStack () throw() {
	    // automatic cleanup
	}

	std::string to_string() const {
	    std::ostringstream os;
	    for (size_t i = 0; i < stack.size(); i++)
		os << stack[i].to_string() << std::endl;
	    return os.str();
	}

	std::vector<StackEntry> stack;
    };

}  // namespace bnc

#endif /* BACKTRACE_H */
