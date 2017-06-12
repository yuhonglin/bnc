#ifndef OPTIM_H
#define OPTIM_H

#include <iostream>

#include <matrix/matrix.hh>

namespace bnc {

  namespace optim {

    enum INFO {
      SUCCESS = 0,
      ABNORMAL = 1,
      INVALIDINPUT = 2,
      MAXITERREACHED = 3
    };

    struct Result {
      Vector   x;
      double   f;
      int  nIter;
      int   info;
    };

    std::ostream& operator<< (std::ostream & out,  const Result& data) {
      out << "x     = " << data.x.transpose() << std::endl;
      out << "f     = " << data.f             << std::endl;
      out << "nIter = " << data.nIter         << std::endl;
      switch(data.info) {
      case SUCCESS:
	out << "info  = SUCCESS" << std::endl;
	break;
      case ABNORMAL:
	out << "info  = ABNORMAL" << std::endl;
	break;
      case INVALIDINPUT:
	out << "info  = INVALIDINPUT" << std::endl;
	break;
      case MAXITERREACHED:
	out << "info  = MAXITERREACHED" << std::endl;
	break;
	  
      }
      return out;
    }    
  }  // optim

}  // bnc

#endif /* OPTIM_H */
