#include "matrices.h"
#include <iostream>
using namespace std;
// Set origin to center of mass
void massCenter(float** xyz, int *zeros, int size) {
  float xm, ym, zm;
  xm = ym = zm = 0.;
  int i;

  for( i = 0; i < size; i++ ) {
    if (zeros[i])
      continue;
    xm += xyz[i][0];
    ym += xyz[i][1];
    zm += xyz[i][2];
  }
  xm /= size;
  ym /= size;
  zm /= size;
  for( i = 0; i < size; i++ ) {
    xyz[i][0] -= xm;
    xyz[i][1] -= ym;
    xyz[i][2] -= zm;
  }
}


void align(float** xyzA, float** xyzB, int *zeros, int size){
  massCenter(xyzA, zeros, size);
  massCenter(xyzB, zeros, size);
  // PStruct
  //  double dist;
  double u[4][4];
  double omega[7][7]; // 7 to use fortran like notation in loops
  double o[7][7]; // eigen vector SET BY EIGEN
  double ha[4][4];
  double ka[4][4];
  double r[4][4]; // rotation
  double op[4];
  double s;
  double det;

  double d [6+1]; // eigen val SET BY EIGEN and not used...

  // =====  Allocate rotated structure
  double (*xyzn)[3] = new double [size][3];

  // ===== Compute rotation matrix
  int i;
  for (i = 1; i <= 3; i++) {
    for (int j = 1;  j <= 3; j++) {
      u[i][j] = 0.0;
      //o[i][j] = 0.0;
      ha[i][j] = 0; //
      ka[i][j] = 0; //
      r[i][j] = 0; //
    }
    d[i] = 0.;
  }
  // additional init PF TEST
  for (i = 1; i <= 6; i++) {
    for (int j = 1;  j <= 6; j++) {
      o[i][j] = 0.0;
      // omega done in original code later
    }
    d[i] = 0.;
    op[i] = 0.;
  }
  int k;
  cout << "skipping " << flush;
  for (k = 0; k < size; k++) {
    if (zeros[i])
      cout << i << " " << flush;
    continue;
    u[1][1] = u[1][1] + xyzA[k][0] * xyzB[k][0];
    u[1][2] = u[1][2] + xyzA[k][0] * xyzB[k][1];
    u[1][3] = u[1][3] + xyzA[k][0] * xyzB[k][2];
    u[2][1] = u[2][1] + xyzA[k][1] * xyzB[k][0];
    u[2][2] = u[2][2] + xyzA[k][1] * xyzB[k][1];
    u[2][3] = u[2][3] + xyzA[k][1] * xyzB[k][2];
    u[3][1] = u[3][1] + xyzA[k][2] * xyzB[k][0];
    u[3][2] = u[3][2] + xyzA[k][2] * xyzB[k][1];
    u[3][3] = u[3][3] + xyzA[k][2] * xyzB[k][2];
  }
  cout << endl << flush;

  det = u[1][1] * u[2][2] * u[3][3] +
        u[1][2] * u[2][3] * u[3][1] +
        u[1][3] * u[2][1] * u[3][2] -
        u[1][3] * u[2][2] * u[3][1] -
        u[1][1] * u[2][3] * u[3][2] -
        u[1][2] * u[2][1] * u[3][3];

  //cout << "DET= " << det << endl << flush;
  for (i = 1; i <= 6; i++) {
    for (int j = 1; j <= 6; j++) {
      omega[i][j] = 0.0;
    }
  }
  for (i = 1; i <= 3; i ++) {
    for (int j = 4; j <= 6; j++) {
      omega[i][j] = u[i][j-3];
    }
  }
  for (i = 4; i <= 6; i++) {
    for (int j = 1; j <= 3; j++) {
      omega[i][j] = u[j][i-3];
    }
  }

  // eigen sets o (eigen vector) and d (eigen val)
  eigen(omega,o,d,6);

  // o values set by eigen used below
  // d not used
  for (k=1; k<=3; k++) {
    for (i = 1; i<=3; i++) {
      ha[i][k] = o[i][k];
      ka[i][k] = o[i + 3][k];
    }
  }

  for (k = 1; k<=3; k++) {
    s = 0.;
    for (i = 1; i<=3; i++) {
      s += ha[i][k] * ha[i][k] ;
    }
    s = sqrt(s);
    for (i = 1; i<=3; i++) {
      ha[i][k] = ha[i][k]/s;
    }
  }
  for (k = 1; k <= 3; k++) {
    s = 0.0;
    for (i = 1; i<=3; i++) {
      s += ka[i][k]* ka[i][k];
    }
    s = sqrt(s);
    for (i = 1; i <= 3; i++) {
      ka[i][k] = ka[i][k]/s;
    }
  }

  op[1] = ha[2][1] * ha[3][2]-ha[3][1] * ha[2][2];
  op[2] = ha[3][1] * ha[1][2]-ha[1][1] * ha[3][2];
  op[3] = ha[1][1] * ha[2][2]-ha[2][1] * ha[1][2];
  s = op[1] * ha[1][3] + op[2] * ha[2][3] + op[3] * ha[3][3];
  if(s < 0.) {
    for (k=1; k <= 3; k++) {
      ha[k][3] = -ha[k][3];
    }
  }
  op[1] = ka[2][1] * ka[3][2]-ka[3][1] * ka[2][2];
  op[2] = ka[3][1] * ka[1][2]-ka[1][1] * ka[3][2];
  op[3] = ka[1][1] * ka[2][2]-ka[2][1] * ka[1][2];
  s = op[1] * ka[1][3] + op[2] * ka[2][3] + op[3] * ka[3][3];


  if(det > 0.) {
    if(s < 0.) {
      for (k = 1; k<=3; k++) {
	ka[k][3] = -ka[k][3];
      }
    }
  }
  else {
    if(s > 0.0) {
      for (k = 1; k<=3; k++) {
	ka[k][3] = -ka[k][3];
      }
    }
  }

  k = 1;
  if(det < 0.) k = -1;
  // rotation matrix
  for (i = 1; i<=3; i++) {
    for (int j = 1; j<=3; j++) {
      r[i][j] = ka[i][1] * ha[j][1] + ka[i][2] * ha[j][2] + k * ka[i][3] * ha[j][3];
    }
  }

  for (i=0; i < size; i++) {
    xyzn[i][0]=r[1][1]*xyzA[i][0]+r[1][2]*xyzA[i][1]+r[1][3]*xyzA[i][2];
    xyzn[i][1]=r[2][1]*xyzA[i][0]+r[2][2]*xyzA[i][1]+r[2][3]*xyzA[i][2];
    xyzn[i][2]=r[3][1]*xyzA[i][0]+r[3][2]*xyzA[i][1]+r[3][3]*xyzA[i][2];
    xyzA[i][0] = xyzn[i][0];
    xyzA[i][1] = xyzn[i][1];
    xyzA[i][2] = xyzn[i][2];
  }

  if (xyzn) delete [] xyzn;
}
