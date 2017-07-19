#ifndef TEXT_H
#define TEXT_H

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iomanip>

#include <matrix/matrix.hh>

namespace bnc {
    namespace io {
	namespace text {
	    template<class T>
	    void dump(const T &d,
		      typename std::enable_if<std::is_class<T>::value, std::string>::type fn)
	    {
		std::ofstream ofs(fn);
		ofs << std::setprecision(16) << d;
		ofs.flush();
		ofs.close();
	    }

	    bnc::Matrix load(const std::string& fn) {
		std::ifstream infile;
		infile.open(fn);
		if(!infile.good()) {
		  LOG_WARNING(fn);
		  return Matrix(0,0);
		}
		
		int lineIdx = 0;
		int ncols = 0;
		int prevlen = 0;
		double n;
		std::vector<double> nums;
		std::string line;
		while (std::getline(infile, line)) {
		  std::stringstream lineStream(line);		  
		  while(lineStream >> n) {
		    nums.push_back(n);
		  }
		  if (ncols==0) {
		    ncols = nums.size();
		  } else if (ncols != nums.size()-prevlen) {
		    LOG_WARNING(std::string("not squared: ") + fn
				+ std::string(", return empty matrix"));
		    return Matrix(0,0);
		  }
		  prevlen = nums.size();
		}

		Matrix ret(nums.size()/ncols,ncols);
		for (int i=0; i<nums.size()/ncols; i++) {
		  for (int j=0; j<ncols; j++) {
		    ret(i,j) = nums[i*ncols+j];
		  }
		}
		
		return ret;
	    }
	}  // namespace text
    } // namespace io
}  // namespace bnc

#endif /* TEXT_H */
