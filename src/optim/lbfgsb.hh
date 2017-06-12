#ifndef LBFGSB_H
#define LBFGSB_H

#include <memory>

#include <cstring>

#include <util/logger.hh>
#include <util/fortran.hh>
#include <matrix/matrix.hh>
#include <optim/optim.hh>

extern "C" {
  void FCALL(setulb)(int* n, int* m, double *x, double *l, 
		    double *u, int *nbd, double* f, double *g, 
		    double* factr, double* pgtol, double *wa, 
		    int *iwa, char *task, int* iprint, 
		    char *csave, bool *lsave, int *isave, 
		    double *dsave);
}

namespace bnc {
  namespace optim {
    
    class Lbfgsb
    {
    public:
      enum BTYPE {
        UNBND = 0,
	LOWER = 1,
	BOTH  = 2,
        UPPER = 3 
      };

      static void copyCStrToCharArray (const char* source, char* dest, int ndest=60) {
	// Get the length of the source C string.
	int nsource = strlen(source);
	// Only perform the copy if the source can fit into the destination.
	if (nsource < ndest) {
	  // Copy the string.
	  strcpy(dest,source);
	  // Fill in the rest of the string with blanks.
	  for (int i = nsource; i < ndest; i++)
	    dest[i] = ' ';
	}
      }

      template <typename F, typename G>
      static Result
      min(F fFunc, G gFunc, Vector& x0, Vector& lb,
	  Vector& ub, iVector& btype,
	  double* wa=NULL, int* iwa=NULL,
	  int maxIter=1000, int m=5,
	  double factr=1e7, double pgtol=1e-5,
	  int iprint=-1)
      {
	char   task[60];
	copyCStrToCharArray("START", task);
	char   csave[60];
	bool   lsave[4];
	int    isave[44];
	double dsave[29];

	Vector g = Vector::Zero(x0.size());    // gradient
	Result ret;
	ret.x = x0;
	ret.nIter = 0;
	ret.f = 0.;
	ret.info = SUCCESS;

	int n = x0.size();
	
	// prepare working array
	//if (wa == NULL)
	//  wa = std::unique_ptr<double>
	//    (new double[(2*m + 4)*n + 11*m*m + 8*m]).get();
	wa = new double[(2*m + 4)*n + 11*m*m + 8*m];

	//if (iwa == NULL)
	//  iwa = std::unique_ptr<int>(new int[3*n]).get();
	iwa = new int[3*n];

	FCALL(setulb)(&n, &m, ret.x.data(), lb.data(),
		      ub.data(), btype.data(), &(ret.f), g.data(),
		      &factr, &pgtol, wa, iwa, task, &iprint,
		      csave, lsave, isave, dsave);

	for (int i=0; i<maxIter; i++) {
	  if (task[0] == 'F') {
	    // evaluate F and G functions
	    ret.f = fFunc(ret.x);
	    g = gFunc(ret.x);
	  } else if (task[0] == 'N') {
	    // new iteration
	    ret.nIter++;
	    if (ret.nIter == maxIter) {
	      copyCStrToCharArray("STOP", task);
	      ret.info = MAXITERREACHED;
	      break;
	    }
	  } else if (task[0] == 'C') {
	    // converged
	    break;
	  } else if (task[0] == 'A') {
	    // abnormal
	    ret.info = ABNORMAL;
	    break;
	  } else if (task[0] == 'E') {
	    ret.info = INVALIDINPUT;
	    break;
	  }
	  FCALL(setulb)(&n, &m, ret.x.data(), lb.data(),
			ub.data(), btype.data(), &(ret.f), g.data(),
			&factr, &pgtol, wa, iwa, task, &iprint,
			csave, lsave, isave, dsave);
	}

	return ret;
      } // function solve()
      
    };

  }  // optim
}  // bnc

#undef int

#endif /* LBFGSB_H */
