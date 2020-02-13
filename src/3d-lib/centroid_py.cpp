#include "Python.h"
#include "3dStats.h"
#include <iostream>
// #include <string>
// using namespace std;
// cout << "START" << endl << flush;

#if PY_MAJOR_VERSION >= 3
    #define PyInt_FromLong PyLong_FromLong
#endif

/* The function doc string */
PyDoc_STRVAR(centroid_wrapper__doc__,
"From a list of lists of xyz positions, and a given threshold (nm), return the \n\
number of equivalent positions, the RMSD and the dRMSD.\n\
   :param xyzs: list of lists of tuples of (x, y, z ) coordinates.\n\
   :param dcutoff: distance cutoff to consider 2 particles as equivalent \n\
      in position (nm)\n\
   :param nmodels: number of models or list of lists passed as first argument\n\
   :param verbose: prints the distance of each model to average model (in stderr)\n\
   :param getavg: return a list for each x, y, z coordinates, representing the average model\n\
\n\
   :returns: the index of the model that is found to be the centroid\n\
");


static PyObject* centroid_wrapper(PyObject* self, PyObject* args)
{
  PyObject *py_xs;
  PyObject *py_ys;
  PyObject *py_zs;
  PyObject *py_zeros;
  int size;
  int nmodels;
  int verbose;
  int getavg;

  if (!PyArg_ParseTuple(args, "OOOOiiii", &py_xs, &py_ys, &py_zs, &py_zeros, &size, 
			&nmodels, &verbose, &getavg))
    return NULL;
 
  float **xyz;
  int zeros[size];
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


  for (i=0; i<size; i++)
    zeros[i] = PyObject_IsTrue(PyTuple_GET_ITEM(py_zeros, i));

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
    avgCoord(it1, it2, zeros, size, modelList, add_first, avg);
    add_first = 0;
    numP++;
    // cout << "\n";
  }

  for (int i = 0; i < size; ++i) {
    for (int j = 0; j < 3; ++j) {
      avg[i][j] /= numP;
    }
    // cout << "   p" << i << " " << i << " " << avg[i][0] << " " << avg[i][1] << " " << avg[i][2]  << endl;
  }

  for (it1=xyzlist.begin(); it1!=xyzlist.end(); it1++) {
    dist2Avg = findCenrtroid(it1, avg, size);
    dist2Centroid.insert(make_pair(dist2Avg, it1->first));
  }

  if (verbose){
    for (it3=dist2Centroid.begin(); it3!=dist2Centroid.end(); it3++) {
      cerr << it3->second << " rmsd2avg " << it3->first << endl;
    }
  }

  for (int i=0; i<size; i++) {
    delete[] xyz[i];
  }
  delete[] xyz;

  // give it to me

  if (getavg){
    PyObject * py_result = NULL;
    PyObject * py_subresult = NULL;
    py_result = PyList_New(3);
    i=0;
    for (int j = 0; j < 3; ++j) {
      py_subresult = PyList_New(size);
      for (int i = 0; i < size; ++i) {
	PyList_SetItem(py_subresult, i, PyFloat_FromDouble(avg[i][j]));
      }
      PyList_SetItem(py_result, j, py_subresult);
    }
    for (int i=0; i<size; i++) {
      delete[] avg[i];
    }
    delete[] avg;
    
    return py_result;
  }

  return PyInt_FromLong(atoi(dist2Centroid.begin()->second.c_str()));
}


 
static PyMethodDef centroidMethods[] =
  {
    {"centroid_wrapper", centroid_wrapper, METH_VARARGS, 
    centroid_wrapper__doc__},
    {NULL, NULL, 0, NULL}
  };

#if PY_MAJOR_VERSION >= 3
  #define MOD_ERROR_VAL NULL
  #define MOD_SUCCESS_VAL(val) val
  #define MOD_INIT(name) PyMODINIT_FUNC PyInit_##name(void)
  #define MOD_DEF(ob, name, doc, methods) \
          static struct PyModuleDef moduledef = { \
            PyModuleDef_HEAD_INIT, name, doc, -1, methods, }; \
          ob = PyModule_Create(&moduledef);
#else
  #define MOD_ERROR_VAL
  #define MOD_SUCCESS_VAL(val)
  #define MOD_INIT(name) PyMODINIT_FUNC init##name(void)
  #define MOD_DEF(ob, name, doc, methods) \
          ob = Py_InitModule3(name, methods, doc);
#endif

MOD_INIT(centroid) {

	PyObject *m;

	MOD_DEF(m, "centroid", "Functions to get the centroid of a given group of models.",
			centroidMethods)
	if (m == NULL)
		return MOD_ERROR_VAL;

	return MOD_SUCCESS_VAL(m);

}
