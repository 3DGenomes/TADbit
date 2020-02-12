#include <sstream>
#include <map>
#include <set>
#include <math.h>
#include "align.h"
#include <cstring>

// #include <iostream>
using namespace std;

float distance(float* p1, float* p2) {
  float x, y, z;

  x = p1[0] - p2[0];
  y = p1[1] - p2[1];
  z = p1[2] - p2[2];
  
  return sqrt(x*x + y*y + z*z);
}

float ddistance(float* p1, float* p2) {
  float x, y, z;

  x = p1[0] - p2[0];
  y = p1[1] - p2[1];
  z = p1[2] - p2[2];
  
  return x*x + y*y + z*z;
}


void avgCoord(map<string, float**>::iterator it1, map<string, float**>::iterator it2,
	      int *zeros, int size, set<string> modelList, bool add_first, float** avg) {
  bool is_in; 

  align(it2->second, it1->second, zeros, size);
  // PStruct
  if(add_first) {
    // cout << it2->first << endl;
    for (int i=0; i < size; i++) {
      avg[i][0] += it1->second[i][0];
      avg[i][1] += it1->second[i][1];
      avg[i][2] += it1->second[i][2];
    }
  }
  is_in = modelList.find(it2->first) != modelList.end();
  if(!is_in) {
    // cout << it2->first << endl;
    for (int i=0; i < size; i++) {
      avg[i][0] += it2->second[i][0];
      avg[i][1] += it2->second[i][1];
      avg[i][2] += it2->second[i][2];
    }
   modelList.insert(it2->first);
  }
}


void rmsdRMSD(float** xyzA, float** xyzB, int *zeros, int size, float thres, 
	      int &eqv, float &rms, float &drms) {
  float dist;
  int last;
  dist = .0;
  rms = 0.;
  drms = .0;
  eqv = 0;
  last = size - 1;
  thres *= thres;

  align(xyzA, xyzB, zeros, size);
  // PStruct

  // rmsd last particle in the model since loop1 stops at the second last
  dist = ddistance(xyzA[last],xyzB[last]);

  // cout <<"end"<<endl;
  //dist = distance(xyzn[last],xyzB[last]) * distance(xyzn[last],xyzB[last]);
  if (dist < thres) eqv++;
  rms += dist;
  // rmsd remaining particles
  for (int i=0; i < size-1; i++) {
    dist = ddistance(xyzA[i], xyzB[i]);
    //dist = distance(xyzn[i],xyzB[i]) * distance(xyzn[i],xyzB[i]);
    if (dist < thres) eqv++;
    rms += dist;
    // drmsd
    for (int j=i+1; j < size; j++) {
	dist = distance(xyzA[i],xyzA[j]) - distance(xyzB[i],xyzB[j]);
	drms += dist*dist;
    }
  }
  drms = sqrt(drms / (size*(size-1)/2));
  rms  = sqrt(rms / size);
}


void consistency(float** xyzA, float** xyzB, int *zeros, int size, float thres, 
		 int * &cons_list) {
  float dist;
  int last;
  dist = .0;
  last = size - 1;
  thres *= thres;
  align(xyzA, xyzB, zeros, size);
  // PStruct

  // rmsd last particle in the model since loop1 stops at the second last
  dist = ddistance(xyzA[last],xyzB[last]);

  if (dist < thres) cons_list[last]++;
  // eqv remaining particles
  // cout <<"start"<<endl;
  for (int i=0; i < size; i++) {
    dist = ddistance(xyzA[i],xyzB[i]);
    if (dist < thres) {
	cons_list[i]=1;
    }else{cons_list[i]=0;}
  }
}


// void findCenrtroid (map<string, float**>::iterator it1, float** avg, int size, map<string, float> dist2Centroid ) {
float findCenrtroid (map<string, float**>::iterator it1, float** avg, int size) {
  float dist, rms;
  dist = .0;
  rms = .0;  

  // rmsd
  for (int i=0; i < size; i++) {
    dist = ddistance(avg[i],it1->second[i]);
    rms += dist;
  }
  rms = sqrt(rms / size);

  return rms;
}


float** populateMap(int size, float** xyz) {
  float **tmpxyz;

  tmpxyz = new float*[size];
  for(int i=0; i<size; i++) {
    tmpxyz[i] = new float[3];
    memset(tmpxyz[i], 0, 3*sizeof(float));
  }
  for(int i=0; i<size; i++) {
    for (int j=0; j<3; j++) {
      tmpxyz[i][j] = xyz[i][j];
    }
  }

  return tmpxyz;

//  for (int i=0; i<size; i++) {
//    delete[] tmpxyz[i];
//  }
//  delete[] tmpxyz;

}

