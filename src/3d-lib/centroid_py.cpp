#include "Python.h"
#include "centroid.cpp"
#include <iostream>
#include <string>
using namespace std;
// cout << "START" << endl << flush;

/* The function doc string */
PyDoc_STRVAR(centroid_wrapper__doc__,
"From a list of lists of xyz positions, and a given threshold (nm), return the \n\
number of equivalent positions, the RMSD and the dRMSD.\n\
   :param xyzs: list of lists of tuples of (x, y, z ) coordinates.\n\
   :param dcutoff: distance cutoff to consider 2 particles as equivalent \n\
      in position (nm)\n\
   :param nmodels: number of models or list of lists passed as first argument\n\
\n\
   :returns: the index of the model that is found to be the centroid\n\
");


static PyObject* centroid_wrapper(PyObject* self, PyObject* args)
{
  PyObject *py_xs;
  PyObject *py_ys;
  PyObject *py_zs;
  int size;
  int nmodels;
  int verbose;

  if (!PyArg_ParseTuple(args, "OOOiii", &py_xs, &py_ys, &py_zs, &size, &nmodels, &verbose))
    return NULL;
 
  float **xyz;
  int i;
  int j;
  int numP;
  float dist2Avg;
  map<string, float**>::iterator it1;
  map<string, float**>::iterator it2;
  map<float, string>::iterator it3;
  float **avg;
  bool add_first;
  map<float, string> dist2Centroid;
  set<string> modelList;
  string modelId;
  ostringstream tmpStr;


  //map<string, float**> xyzlist;
  map<string, float**> xyzlist;
  
  xyz = new float*[size];
  for(int i=0; i<size; i++) {
    xyz[i] = new float[3];
    memset(xyz[i], 0, 3*sizeof(float));
  }

  avg = new float*[size];
  for(int i=0; i<size; i++) {
    avg[i] = new float[3];
    memset(avg[i], 0, 3*sizeof(float));
  }

  for (j=0; j<nmodels; j++){

    for (i=0; i<size; i++){
      xyz[i][0] = PyFloat_AS_DOUBLE(PyList_GET_ITEM(PyList_GET_ITEM(py_xs, j), i));
      xyz[i][1] = PyFloat_AS_DOUBLE(PyList_GET_ITEM(PyList_GET_ITEM(py_ys, j), i));
      xyz[i][2] = PyFloat_AS_DOUBLE(PyList_GET_ITEM(PyList_GET_ITEM(py_zs, j), i));
    }
    tmpStr.str("");
    tmpStr.clear();
    tmpStr << j;
    modelId = tmpStr.str();
    // cout << modelId << " model id "  << xyz[1][0]<<" "<<xyz[1][1]<<" "<<xyz[1][2]<<endl;
    xyzlist.insert(make_pair(modelId, populateMap(size, xyz)));
  }

  numP = 1; 
  add_first = 1;
  it1=xyzlist.begin();
  for ((it2=it1)++; it2!=xyzlist.end(); it2++) {
    //cout << it1->first << " " << it2->first << endl;
    // avgCoord(it1->second, it2->second, size, avg);
    avgCoord(it1, it2, size, modelList, add_first, avg);
    add_first = 0;
    numP++;
    // cout << "\n";
  }

  for (int i = 0; i < size; ++i) {
    for (int j = 0; j < 3; ++j) {
      avg[i][j] /= numP;
    }
  }

  for (it1=xyzlist.begin(); it1!=xyzlist.end(); it1++) {
    dist2Avg = findCenrtroid(it1, avg, size);
    dist2Centroid.insert(make_pair(dist2Avg, it1->first));
  }

  if (verbose){
    for (it3=dist2Centroid.begin(); it3!=dist2Centroid.end(); it3++) {
      cout << it3->second << " rmsd2avg " << it3->first << endl;
    }
  }

  for (int i=0; i<size; i++) {
    delete[] xyz[i];
  }
  delete[] xyz;

  for (int i=0; i<size; i++) {
    delete[] avg[i];
  }
  delete[] avg;

  // give it to me
  // PyObject * py_result = NULL;
  // py_result = PyDict_New();
  // i=0;
  // for (it3=dist2Centroid.begin(); it3!=dist2Centroid.end(); it3++) {
  //   PyDict_SetItem(py_result, PyInt_FromLong(atoi(it3->second.c_str())), PyFloat_FromDouble(it3->first));
  //   i++;
  // }
  // return py_result;
  return PyInt_FromLong(atoi(dist2Centroid.begin()->second.c_str()));
}


 
static PyMethodDef centroidMethods[] =
  {
    {"centroid_wrapper", centroid_wrapper, METH_VARARGS, 
    centroid_wrapper__doc__},
    {NULL, NULL, 0, NULL}
  };

PyMODINIT_FUNC
 
initcentroid(void)
{
  (void) Py_InitModule3("centroid", centroidMethods, 
			"Functions to get the centroid of a given group of models.");
}
