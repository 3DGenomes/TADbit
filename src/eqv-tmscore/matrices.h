#ifndef __matrices_h
#define __matrices_h 1

#include <math.h>
//#include "error.h"

//7 instead of 6 to have indexes 1<=i<=size

//extern void eigen ( double mat[7][7], double eigvec[7][7], double*  eigval, int ndim);
extern void eigen ( double (* mat)[7], double (* eigvec)[7], double*  eigval, int ndim);

//extern void jacobi(double a[7][7], int n, int np, double *d, double v[7][7], int& nrot, double tiny, int& idiag);
// Y extern void jacobi(double (* a)[7], int n, int np, double *d, double (*v)[7], int& nrot, double tiny, int& idiag);
extern void jacobi(double (* a)[7], int n, int np, double *d, double (*v)[7], int nrot, double tiny, int idiag);

// Sort eigenvalues in ascending order.
// Eigenvectors are ordered correspondingly.
//extern void eigsrt (double* d, double  v[7][7], int n,int np);
extern void eigsrt (double* d, double  (*v)[7], int n,int np);

#endif
