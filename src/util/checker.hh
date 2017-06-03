#ifndef CHECKER_H
#define CHECKER_H

#include <string>

#include <util/logger.hh>

namespace bnc {

#ifndef BNC_FAST_MODE
    // in safe mode, we checkfor and checkfor and checkfor...
    template<typename T>
    inline bool _check_for_helper(T t)
    {
	if (t) return true;
	else return false;
    }
    
    template<typename T, typename... Args>
    inline bool _check_for_helper(const T& t, Args... args)
    {
	if (t) return _check_for_helper(args...);
	else return false;
    }
    
    template<typename... Args>
    inline void check_for(const std::string& info, Args... args) {
	if(!_check_for_helper(args...))
	{
	    LOG_ERROR(info);
	}
    }
#else
    // in fast mode, do nothing
    template<typename... Args>
    inline void check_for(const std::string& info, Args... args) {
	return;
    }
#endif
    

}  // namespace bnc

#endif /* CHECKER_H */
