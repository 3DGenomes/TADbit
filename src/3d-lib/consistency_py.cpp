#include "Python.h"
#include "3dStats.h"
// #include <iostream>
// using namespace std;

 //cout << "START" << endl << flush;

/* The function doc string */
PyDoc_STRVAR(consistency_wrapper__doc__,
"From 2 lists of xyz positions, and a given threshold (nm), return the \n\
number of equivalent positions, the RMSD and the dRMSD.\n\
   :param xyz_a: first list of tuples of (x, y, z ) coordinates.\n\
   :param xyz_b: second list of tuples of (x, y, z ) coordinates.\n\
   :param dcutoff: distance cutoff to consider 2 particles as equivalent \n\
      in position (nm)\n\
   :param nmodels: number of models passed\n\
\n\
   :returns: a list with the number of equivalent positions, the RMSD and \n\
      the dRMSD, of the alignment. If consistency is True, returns list of \n\
      0 or 1 if a given particle is in equivalent position in both strands.\n\
");

float maximumValue(float *vals, int size)
{
  float max = vals[0];
  for(int i = 1; i<size; i++)
    if(vals[i] > max)
      max = vals[i];
  return max;
}

static PyObject* consistency_wrapper(PyObject* self, PyObject* args)
{
  PyObject **py_xs;
  PyObject **py_ys;
  PyObject **py_zs;
  PyObject **py_models;
  int size;
  int nmodels;
  float thres;
  //cout << "START" << endl << flush;
 
  if (!PyArg_ParseTuple(args, "OOOifOi", &py_xs, &py_ys, &py_zs, &size, &thres, &py_models, &nmodels))
    return NULL;
 
  float ***xyzn;
  int   *cons_list;
  int **scores;
  int i;
  int j;
  int jj;
  int k;
  int msize;
  //cout << "START1" << endl << flush;

  msize = nmodels*(nmodels-1)/2;
  xyzn = new float **[nmodels];
  for (i=0;i<nmodels;i++)
    xyzn[i] = new float *[size];
  //cout << "START2" << endl << flush;


  for (j=0; j<nmodels; j++){
    for (i=0; i<size; i++){
      xyzn[j][i] = new float[3];
      memset(xyzn[j][i], 0, 3*sizeof(float));
      xyzn[j][i][0] = PyFloat_AS_DOUBLE(PyList_GET_ITEM(PyList_GET_ITEM(py_xs, j ), i));
      xyzn[j][i][1] = PyFloat_AS_DOUBLE(PyList_GET_ITEM(PyList_GET_ITEM(py_ys, j ), i));
      xyzn[j][i][2] = PyFloat_AS_DOUBLE(PyList_GET_ITEM(PyList_GET_ITEM(py_zs, j ), i));
    }
  }
  //cout << "START3" << endl << flush;
  scores = new int*[msize];

  k = 0;
  for (j=0; j<nmodels-1; j++){
    for (jj=j+1; jj<nmodels; jj++){
      cons_list = new int[size];
      scores[k] = new int[size];
      consistency(xyzn[j], xyzn[jj], size, thres, cons_list);
      scores[k] = cons_list;
      k++;
      // scores[j+jj-1] = cons_list;
    }
  }

  //cout << "START4" << endl << flush;
  PyObject * py_result = NULL;
  PyObject * py_subresult = NULL;
  py_result = PyList_New(msize);
  k = 0;
  for (j=0; j<nmodels; j++){
    for (jj=j+1; jj<nmodels; jj++){
	py_subresult = PyList_New(size);
	for (i=0; i<size; i++)
	  PyList_SetItem(py_subresult, i, PyInt_FromLong(scores[k][i]));
	PyList_SetItem(py_result,k , py_subresult);
	k++;
    }
  }

  // free
  delete[] cons_list;
  //cout << "START5" << endl << flush;
  for (int j=0; j<nmodels; j++){
    for (int i=0; i<size; i++)
      delete[] xyzn[j][i];
    delete[] xyzn[j];
  }
  delete[] xyzn;
  
  for (int i=0; i<msize-1; i++){
    //cout << i << " "<<msize<<endl << flush;
    delete[] scores[i];
  }
  //cout << "START6" << endl << flush;
  delete[] scores;

  //cout << "START7" << endl << flush;
  //cout << "START5" << endl << flush;
  Py_DECREF(py_xs);
  //cout << "START5" << endl << flush;
  Py_DECREF(py_ys);
  Py_DECREF(py_zs);
  Py_DECREF(py_models);
  Py_CLEAR(py_xs);
  Py_CLEAR(py_ys);
  Py_CLEAR(py_zs);
  Py_CLEAR(py_models);
  
  // give it to me
  return py_result;
}
 
static PyMethodDef ConsistencyMethods[] =
  {
    {"consistency_wrapper", consistency_wrapper, METH_VARARGS, 
    consistency_wrapper__doc__},
    {NULL, NULL, 0, NULL}
  };

PyMODINIT_FUNC
 
initconsistency(void)
{
  (void) Py_InitModule3("consistency", ConsistencyMethods, 
			"Functions to compaire two Chromatin strands.");
}
