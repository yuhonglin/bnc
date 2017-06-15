#ifndef PROFILING_H
#define PROFILING_H

#include <stack>
#include <ctime>

namespace bnc {
    std::stack<std::clock_t> _bnc_prof_tstart;

    void tic() {
	_bnc_prof_tstart.push(clock());
    }

    void toc() {
	std::clock_t tend = clock();
	std::cout << "toc: " << (static_cast<double>(tend-_bnc_prof_tstart.top())/CLOCKS_PER_SEC)
		  << " second(s)." << std::endl;
	_bnc_prof_tstart.pop();
    }
    
}  // namespace bnc

#endif /* PROFILING_H */
