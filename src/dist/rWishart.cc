/*
 * This file is used to test rwishart.
 */
#include <cmath>

using namespace std;

#include "dist/dist.hh"
#include "dist/wishart.hh"

#include "rng/rng.hh"

#include "io/text.hh"

#include "dist/mvnorm.hh"

#include <Eigen/Eigenvalues>
#include <cstdlib>

#define F77_CALL(x) x

using namespace bnc;

extern "C" {
    void F77_CALL(dtrmm)(char* SIDE, char* UPLO, char* TRANSA,
			 char* DIAG, int* M, int* N, double *ALPHA,
			 double* A, int* LDA, double* B, int* LDB);
    void F77_CALL(dsyrk)(char* UPLO, char *TRANS, int* N, int* K,
			 double* ALPHA, double* A, int* LDA, double* BETA,
			 double* C, int* LDC);
    void F77_CALL(dpotrf)(char*, int*, double*, int*, int*);
}

template<class RNGType>
static double
*std_rWishart_factor(double nu, int p, int upper, double ans[], RNGType* rng)
{
    int pp1 = p + 1;

    memset(ans, 0, p * p * sizeof(double));
    for (int j = 0; j < p; j++) {	/* jth column */
	ans[j * pp1] = sqrt(rchisq(nu - j, rng));
	for (int i = 0; i < j; i++) {
	    int uind = i + j * p, /* upper triangle index */
		lind = j + i * p; /* lower triangle index */
	    ans[(upper ? uind : lind)] = rng->normal();
	    ans[(upper ? lind : uind)] = 0;
	}
    }
    return ans;
}

template<class RNGType>
Matrix rWishart(const int& nu, Matrix scale, RNGType* rng)
{
    double* scal = const_cast<double*>(scale.data());
    char R[2]="R";
    char U[2]="U";
    char N[2]="N";
    char T[2]="T";
    Matrix ret(scale.rows(),scale.rows());
    double *ans = const_cast<double*>(ret.data());
    int dims[2];
    dims[0] = scale.rows();
    dims[1] = scale.rows();
    double *scCp = new double[scale.rows()*scale.rows()];
    double *tmp = new double[scale.rows()*scale.rows()];
    double one = 1;
    double zero = 0;
    int info;
    int psqr = scale.rows()*scale.rows();
    for (int i=0; i<scale.rows()*scale.rows(); i++) {
	scCp[i] = scal[i];
    }

    F77_CALL(dpotrf)(U, &(dims[0]), scCp, &(dims[0]), &info);
    cout << info << endl;
    double *ansp = ans;

    for (int j = 0; j < 1; j++) {
	double *ansj = ansp + j * psqr;
	std_rWishart_factor(nu, dims[0], 1, tmp, rng);

	cout << tmp[0] << ' ' << tmp[3] << ' ' << tmp[6] << endl;
	cout << tmp[1] << ' ' << tmp[4] << ' ' << tmp[7] << endl;
	cout << tmp[2] << ' ' << tmp[5] << ' ' << tmp[8] << endl;	

	F77_CALL(dtrmm)(R, U, N, N, dims, dims,
			&one, scCp, dims, tmp, dims);
	cout << "--------------------------------------" << endl;
	cout << tmp[0] << ' ' << tmp[3] << ' ' << tmp[6] << endl;
	cout << tmp[1] << ' ' << tmp[4] << ' ' << tmp[7] << endl;
	cout << tmp[2] << ' ' << tmp[5] << ' ' << tmp[8] << endl;	
	
	F77_CALL(dsyrk)(U, T, &(dims[1]), &(dims[1]),
			&one, tmp, &(dims[1]),
			&zero, ansj, &(dims[1]));

//	for (int i = 1; i < dims[0]; i++)
//	    for (int k = 0; k < i; k++)
//		ansj[i + k * dims[0]] = ansj[k + i * dims[0]];
    }

    return ret;
}


int main(int argc, char *argv[])
{
    Matrix m(3,3);
    m << 10,1,0,1,5,0,0,0,1;
    MersenneTwisterRNG<AHRENS_DIETER> rng(10);

//    for(int i=0; i<10; i++)
//	cout << rng.normal() << endl;
    
    auto n = rWishart(10, m.inverse(), &rng);
    cout << n << endl;
    return 0;
}

