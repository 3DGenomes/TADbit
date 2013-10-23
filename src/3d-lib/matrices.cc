#include "matrices.h"

// eigen values
// call in rmsdist:
// eigen(omega,o,d,7);

//void eigen ( double mat[7][7], double eigvec[7][7], double*  eigval, int ndim)
void eigen ( double (* mat)[7], double (* eigvec)[7], double*  eigval, int ndim)
{
	int i, j, n;
	int nrotval = 0;
	int idiagval = 0;
//Y	int &nrot = nrotval; // passed by ref but not initialized in fortran
//	int &idiag = idiagval; // passed by ref but not initialized in fortran
	int nrot = nrotval; // passed by ref but not initialized in fortran
	int idiag = idiagval; // passed by ref but not initialized in fortran
	double tiny, amat[7][7];
	
	tiny = 1.0e-20;
	n = ndim;
	
	// copy matrix mat in amat
	for (i=1; i <= ndim; i++) {
		for (j=1; j<= ndim; j++) {
			amat[i][j] = mat[i][j];
//cout << "eigen i j mat--" << i << "--" << j << "--" << amat[i][j] << endl;
		}
	}
	
	jacobi(amat, ndim, n, eigval, eigvec, nrot, tiny, idiag);

/*
for (i=1;i<=6;i++) {	
	cout << "eigvec i j eigval i --" << i;
	for (j=1;j<=6;j++) {
		cout << "--" << eigvec[j][i];
	}
	cout << "=======" << eigval[i] << endl;
}
*/

	// sort
	eigsrt(eigval, eigvec, ndim, n);
//for (i=1;i<=6;i++) {	
//	cout << "eigvec i j eigval i --" << i;
//	for (j=1;j<=6;j++) {
//		cout << "--" << eigvec[j][i];
//	}
//	cout << "=======" << eigval[i] << endl;
//}
}

// call in eigen 
//void jacobi(          amat,  ndim,      n,    eigval,        eigvec,       nrot,        tiny,      idiag);
//void jacobi(double a[7][7], int n, int np, double *d, double v[7][7], int& nrot, double tiny, int& idiag)
//Y void jacobi(double (* a)[7], int n, int np, double *d, 
//			double (* v)[7], int& nrot, double tiny, int& idiag)
void jacobi(double (* a)[7], int n, int np, double *d, 
			double (* v)[7], int nrot, double tiny, int idiag)
{
	
	const int maxpw = 129;
	int i, j, ip, iq;
	double b[maxpw+1];
	double z[maxpw+1];
	double sm;
	double tresh;
	double g;
	double h;
	double t;
	double theta;
	double c;
	double s, tau;

	for (ip=1; ip<=n; ip++) {
		for (iq=1; iq<=n; iq++) {
			v[ip][iq] = 0.;
		}
		v[ip][ip] = 1.;
	}
	for (ip=1; ip<=n; ip++) {
		b[ip] = a[ip][ip];
		d[ip] = b[ip];
		z[ip] = 0.;
//cout << "jac1 ip b d z --" << ip << "--" << b[ip] << "--" << d[ip] << "--" << z[ip] << endl;
	}
	nrot = 0;
	for (i=1; i<=50; i++) {
		sm = 0.;
		for (ip=1; ip<= n-1; ip++) {
			for (iq=ip+1; iq<=n; iq++) {
				sm += fabs(a[ip][iq]);
			}
		}
		idiag = i;
		if (sm < tiny) {
			return;
		}
		if (i < 4) {
			tresh = 0.2 * sm /(n*n); 
		} else {
			tresh = 0;
		}
		for (ip=1; ip<=n-1; ip++) {
			for (iq=ip+1; iq<=n; iq++) {
				//g = 1.0d2 XX
				g = 1e2 * fabs(a[ip][iq]);
				if ( (i > 4) && ( fabs(d[ip]) + g == fabs(d[ip]) )
					&& ( fabs(d[iq]) + g == fabs(d[iq]) ) ) {
					a [ip][iq] = 0.;
				}
				else if ( fabs(a[ip][iq]) > tresh ) {
					h = d[iq] - d[ip];
					if (fabs(h) + g == fabs(h)) {
						t = a[ip][iq]/h;
					}
					else {
						theta = 0.5 * h / a[ip][iq];
						t = 1. / ( fabs(theta) + sqrt(theta * theta + 1.) );
						if (theta < 0.) t = -t;
					}
					c = 1. / sqrt(1 + t * t);
					s = t * c;
					tau = s / (1. + c);
					h = t * a[ip][iq];
					z [ip] = z[ip] - h;
					z [iq] = z[iq] + h; 	
					d [ip] = d[ip] - h;
					d [iq] = d[iq] + h; 	
					a [ip][iq] = 0.;
					for (j=1; j<=ip-1; j++) {
						g = a[j][ip];
						h = a[j][iq] ;
						a[j][ip] = g - s * (h + g * tau);
						a[j][iq] = h + s * ( g - h * tau);
					}
					for (j = ip+1; j <= iq-1; j++) {
						g = a[ip][j];
						h = a[j][iq] ;
						a[ip][j] = g - s * (h + g * tau);
						a[j][iq] = h + s * (g - h * tau);
					}
					for (j = iq+1; j <= n; j++) {
						g = a[ip][j];
						h = a[iq][j] ;
						a[ip][j] = g - s * (h + g * tau);
						a[iq][j] = h + s * (g - h * tau);
					}
					for (j=1; j <= n ; j++) {
						g = v[j][ip];
						h = v[j][iq] ;
						v[j][ip] = g - s * (h + g * tau);
						v[j][iq] = h + s * (g - h * tau);
//cout << "jacobi j ip iq v v" << j << "--" << ip << "--" << iq << "--" << v[j][ip] << "--" << v[j][iq]<< endl;
					}
					nrot = nrot + 1;
				}
			}
		}
		for (ip=1; ip <= n; ip++) {
			b[ip] = b[ip] + z[ip];
			d[ip] = b[ip];
			z[ip] = 0.;
		}
	}
//	error("jacobi : no convergence in 50 iterations. stopping...");
}

// 
// Sort eigenvalues in ascending order.
// Eigenvectors are ordered correspondingly.
//

// call in eigen:
// eigstr(eigval, eigvec, ndim, n);
//void eigsrt (double* d, double  v[7][7], int n,int np)
void eigsrt (double* d, double  (*v)[7], int n,int np)
{
	//implicit real*8 (a-h,o-z)
	//dimension d(np), v(np,np)
	int i, k, j;
	double p;
	for (i=1; i <= n-1; i++) {
		k = i;
		p = d[i];
		for (j=i+1; j <= n; j++) { 
			if (d[j] >= p) {
				k = j;
				p = d[j];
			}
		}
		if (k != i) {
			d[k] = d[i];
			d[i] = p;
			for (int j=1; j<=n; j++) {
				p = v[j][i];
				v[j][i] = v[j][k];
				v[j][k] = p;
			}
		}
	}
	return;
}

