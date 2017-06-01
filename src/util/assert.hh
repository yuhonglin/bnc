#ifndef ASSERT_H
#define ASSERT_H

#include <string>

using namespace std;

#include <util/logger.hh>

namespace bnc {
    template<typename T>
    inline bool _assert_helper(T t)
    {
	if (t) return true;
	else return false;
    }
    
    template<typename T, typename... Args>
    inline bool _assert_helper(T t, Args... args)
    {
	if (t) return _assert_helper;
	else return false;
    }
    
    template<typename... Args>
    void assert(const string& info, Args... args) {
	if(!_assert_helper(args...))
	{
	    LOG_ERROR(info);
	}
    }

}  // namespace bnc

#endif /* CHECKER_H */
