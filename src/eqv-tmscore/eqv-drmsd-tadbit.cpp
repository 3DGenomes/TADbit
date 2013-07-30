#include "matrices.cc"
// #include <iostream>
// using namespace std;


float distance(float* p1, float* p2) {
  float x, y, z;

  x = p1[0] - p2[0];
  y = p1[1] - p2[1];
  z = p1[2] - p2[2];
  
  return sqrt(x*x + y*y + z*z);
}

// Set origin to center of mass
void massCenter(float** xyz, int size) {
  float xm, ym, zm;
  xm = ym = zm = 0.;
  int i;

  for( i = 0; i < size; i++ ) {
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

void rmsdRMSD(float** xyzA, float** xyzB, int size, float thres, 
	      int &eqv, float &rms, float &drms, int * &cons_list, int consistency) {
  float dist;
  int last;
  dist = .0;
  rms = 0.;
  drms = .0;
  eqv = 0;
  last = size - 1;


  massCenter(xyzA, size);
  massCenter(xyzB, size);

  // PStruct
  //  double dist;
  double u[4][4];
  double omega[7][7]; // 7 to use fortran like notation in loops
  double o[7][7]; // eigen vector SET BY EIGEN
  double ha[4][4];
  double ka[4][4];
  double oa[7][7];
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
      oa[i][j] = 0.0;
      o[i][j] = 0.0;
      // omega done in original code later
    }
    d[i] = 0.;
    op[i] = 0.;
  }
  int k;
  for (k = 0; k < size; k++) {
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

  det = u[1][1] * u[2][2] * u[3][3] +
    u[1][2] * u[2][3] * u[3][1] +
    u[1][3] * u[2][1] * u[3][2]-
    u[1][3] * u[2][2] * u[3][1]-
    u[1][1] * u[2][3] * u[3][2]-
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
      s = s + ha[i][k] * ha[i][k] ;
    }
    for (i = 1; i<=3; i++) {
      ha[i][k] = ha[i][k]/sqrt(s);
    }
  }
  for (k = 1; k <= 3; k++) {
    s = 0.0;
    for (i = 1; i<=3; i++) {
      s = s + ka[i][k]* ka[i][k];
    }
    for (i = 1; i <= 3; i++) {
      ka[i][k] = ka[i][k]/sqrt(s);
    }
  }

  for (k = 1; k <= 3; k++) {
    for (i = 1; i <= 3; i++) {
      oa[i][k] = ha[i][k];
      oa[i + 3][k] = ka[i][k];
      oa[i][k + 3] = ha[i][k];
      oa[i + 3][k + 3] = -ka[i][k];
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
  }

  // RMS dist
  for (i= 0; i < size; i++) {
    for (int j=0; j<3; j++) {
      xyzA[i][j] = xyzn[i][j];
    }
  }

  if (xyzn) delete [] xyzn;
  // PStruct

  // rmsd last particle in the model since loop1 stops at the second last
  dist = distance(xyzA[last],xyzB[last]) * distance(xyzA[last],xyzB[last]);

  if (consistency==1){
    if (sqrt(dist) < thres) cons_list[last]++;
    // eqv remaining particles
    // cout <<"start"<<endl;
    for (int i=0; i < size; i++) {
      dist = distance(xyzA[i],xyzB[i]) * distance(xyzA[i],xyzB[i]);
      if (sqrt(dist) < thres) {
	cons_list[i]=1;
      }else{cons_list[i]=0;}
    }
    // cout <<"end"<<endl;
  }else{
    //dist = distance(xyzn[last],xyzB[last]) * distance(xyzn[last],xyzB[last]);
    if (sqrt(dist) < thres) eqv++;
    rms += dist;
    // rmsd remaining particles
    for (int i=0; i < size-1; i++) {
      dist = distance(xyzA[i], xyzB[i]) * distance(xyzA[i],xyzB[i]);
      //dist = distance(xyzn[i],xyzB[i]) * distance(xyzn[i],xyzB[i]);
      if (sqrt(dist) < thres) eqv++;
      rms += dist;
      // drmsd
      for (int j=i+1; j < size; j++) {
	dist =(distance(xyzA[i],xyzA[j]) - distance(xyzB[i],xyzB[j]))
	  * (distance(xyzA[i],xyzA[j]) - distance(xyzB[i],xyzB[j]));
	//dist =(distance(xyzn[i],xyzn[j]) - distance(xyzB[i],xyzB[j]))
	//* (distance(xyzn[i],xyzn[j]) - distance(xyzB[i],xyzB[j]));
	drms += dist;
      }
    }
    drms = sqrt(drms / (size*(size-1)/2));
    rms = sqrt(rms / size);
  }
}

