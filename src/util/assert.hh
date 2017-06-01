#ifndef ASSERT_H
#define ASSERT_H

#include <string>

using namespace std;

#include <util/logger.hh>

namespace bnc {

#ifndef FAST_MODE
    // in safe mode, we check and check and check...
    template<typename T>
    inline bool _assert_helper(T t)
    {
	if (t) return true;
	else return false;
    }
    
    template<typename T, typename... Args>
    inline bool _assert_helper(const T& t, Args... args)
    {
	if (t) return _assert_helper(args...);
	else return false;
    }
    
    template<typename... Args>
    inline void assert(const string& info, Args... args) {
	if(!_assert_helper(args...))
	{
	    LOG_ERROR(info);
	}
    }
#else
    // in fast mode, do nothing
    template<typename... Args>
    inline void assert(const string& info, Args... args) {
	return;
    }
#endif
    

}  // namespace bnc

#endif /* CHECKER_H */
