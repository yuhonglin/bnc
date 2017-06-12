#ifndef TEXT_H
#define TEXT_H

#include <fstream>
#include <string>
#include <vector>

#include <matrix/matrix.hh>

namespace bnc {
    namespace io {
	namespace text {
	    template<class T>
	    void dump(const T &d,
		      typename enable_if<is_class<T>::value, std::string>::type fn)
	    {
		std::ofstream ofs(fn);
		ofs << d;
		ofs.flush();
		ofs.close();
	    }

	    bnc::Vector load(const std::string& fn) {
		std::ifstream infile(fn);
		double x;
		std::vector<double> tmp;
		while(infile >> x) {
		    tmp.push_back(x);
		}
		Vector ret(tmp.size());
		for (int i=0; i<tmp.size(); i++) {
		    ret(i) = tmp[i];
		}
		return ret;
	    }
	}  // namespace text
    } // namespace io
}  // namespace bnc

#endif /* TEXT_H */
