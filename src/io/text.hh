#ifndef TEXT_H
#define TEXT_H

#include <fstream>
#include <string>
using namespace std;

namespace bnc {
    namespace io {
	namespace text {
	    template<class T>
	    void dump(const T &d,
		      typename enable_if<is_class<T>::value, string>::type fn)
	    {
		ofstream ofs(fn);
		ofs << d;
		ofs.flush();
		ofs.close();
	    }
	}  // namespace text
    } // namespace io
}  // namespace bnc

#endif /* TEXT_H */
