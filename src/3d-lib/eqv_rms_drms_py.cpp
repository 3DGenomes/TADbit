#include "Python.h"
#include "3dStats.h"
// #include <iostream>
// using namespace std;

// cout << "START" << endl << flush;

/* The function doc string */
PyDoc_STRVAR(rmsdRMSD_wrapper__doc__,
"From 2 lists of xyz positions, and a given threshold (nm), return the \n\
number of equivalent positions, the RMSD and the dRMSD.\n\
   :param xyz_a: first list of tuples of (x, y, z ) coordinates.\n\
   :param xyz_b: second list of tuples of (x, y, z ) coordinates.\n\
   :param dcutoff: distance cutoff to consider 2 particles as equivalent \n\
      in position (nm)\n\
   :param nmodels: number of models passed\n\
\n\
   :returns: a list with the number of equivalent positions, the RMSD and \n\
      the dRMSD, of the alignment.\n\
");

float maximumValue(float *vals, int size)
{
  float max = vals[0];
  for(int i = 1; i<size; i++)
    if(vals[i] > max)
      max = vals[i];
  return max;
}

static PyObject* rmsdRMSD_wrapper(PyObject* self, PyObject* args)
{
  PyObject **py_xs;
  PyObject **py_ys;
  PyObject **py_zs;
  PyObject **py_models;
  int size;
  int one;
  int nmodels;
  float thres;
  char *what;
  int normed;
  // cout << "START" << endl << flush;
 
  if (!PyArg_ParseTuple(args, "OOOifOiisi", &py_xs, &py_ys, &py_zs, &size, &thres, &py_models, &nmodels, &one, &what, &normed))
    return NULL;
 
  float ***xyzn;
  float *nrmsds;
  float *drmsds;
  float *scores;
  float max_normed;
  int i;
  int j;
  int jj;
  int k;
  int msize;
  float rms = 0;
  float drms = 0;
  int   eqv = 0;
  // cout << "START" << endl << flush;

  msize = nmodels*(nmodels-1)/2;
  xyzn = new float **[nmodels];
  for (i=0;i<nmodels;i++)
    xyzn[i] = new float *[size];
  nrmsds = new float[msize];
  drmsds = new float[msize];
  scores = new float[msize];
  // cout << "START" << endl << flush;


  PyObject * py_result = NULL;
  PyObject * py_subresult = NULL;
  py_result = PyDict_New();
  for (j=0; j<nmodels; j++){
    for (i=0; i<size; i++){
      xyzn[j][i] = new float[3];
      memset(xyzn[j][i], 0, 3*sizeof(float));
      xyzn[j][i][0] = PyFloat_AS_DOUBLE(PyList_GET_ITEM(PyList_GET_ITEM(py_xs, j ), i));
      xyzn[j][i][1] = PyFloat_AS_DOUBLE(PyList_GET_ITEM(PyList_GET_ITEM(py_ys, j ), i));
      xyzn[j][i][2] = PyFloat_AS_DOUBLE(PyList_GET_ITEM(PyList_GET_ITEM(py_zs, j ), i));
    }
  }
  // cout << "START2" << endl << flush;

  k = 0;
  for (j=0; j<nmodels; j++){
    for (jj=j+1; jj<nmodels; jj++){
      rms = 0;
	drms = 0;
      eqv = 0;
      rmsdRMSD(xyzn[j], xyzn[jj], size, thres, eqv, rms, drms);
      nrmsds[k] = rms;
      drmsds[k] = drms;
      scores[k] = eqv * drms / rms;
	// cout << j<<" " <<jj<<" " <<eqv<<" " <<drms<<" " <<rms<<" "<<scores[k]<< endl << flush;
	k++;
    }
  }
  // cout << "START5" << endl << flush;
  if (one){
    // free
    for (int j=0; j<nmodels; j++){
      for (int i=0; i<size; i++){
	delete[] xyzn[j][i];
      }
      delete[] xyzn[j];
    }
    delete[] xyzn;
    
    // give it to me
    return PyFloat_FromDouble(drmsds[0]);
  }

  if (strcmp(what,"rmsd")==0){
    k=0;
    for (j=0; j<nmodels; j++){
      for (jj=j+1; jj<nmodels; jj++){
	if (normed)
	  py_subresult = PyFloat_FromDouble(1-nrmsds[k]/maximumValue(nrmsds, msize));
	else
	  py_subresult = PyFloat_FromDouble(nrmsds[k]);
	PyDict_SetItem(py_result, PyTuple_Pack(2, PyList_GET_ITEM(py_models, j ), PyList_GET_ITEM(py_models, jj)), py_subresult);
	PyDict_SetItem(py_result, PyTuple_Pack(2, PyList_GET_ITEM(py_models, jj), PyList_GET_ITEM(py_models, j )), py_subresult);
	k++;
      }
    }
  }else if (strcmp(what,"drmsd")==0){
    k=0;
    for (j=0; j<nmodels; j++){
      for (jj=j+1; jj<nmodels; jj++){
	if (normed)
	  py_subresult = PyFloat_FromDouble(1-drmsds[k]/maximumValue(drmsds, msize));
	else
	  py_subresult = PyFloat_FromDouble(drmsds[k]);
	PyDict_SetItem(py_result, PyTuple_Pack(2, PyList_GET_ITEM(py_models, j ), PyList_GET_ITEM(py_models, jj)), py_subresult);
	PyDict_SetItem(py_result, PyTuple_Pack(2, PyList_GET_ITEM(py_models, jj), PyList_GET_ITEM(py_models, j )), py_subresult);
	k++;
      }
    }
  }else if (strcmp(what,"eqv")==0){
    k=0;
    for (j=0; j<nmodels; j++){
      for (jj=j+1; jj<nmodels; jj++){
	py_subresult = PyFloat_FromDouble(scores[k] * nrmsds[k] / drmsds[k]);
	PyDict_SetItem(py_result, PyTuple_Pack(2, PyList_GET_ITEM(py_models, j ), PyList_GET_ITEM(py_models, jj)), py_subresult);
	PyDict_SetItem(py_result, PyTuple_Pack(2, PyList_GET_ITEM(py_models, jj), PyList_GET_ITEM(py_models, j )), py_subresult);
	k++;
      }
    }
  }else if (strcmp(what,"score")==0){
    // cout << "START6" << endl << flush;
    max_normed = maximumValue(nrmsds, msize) / maximumValue(drmsds, msize);
    // cout << "START7" << endl << flush;
    k=0;
    for (j=0; j<nmodels; j++){
      for (jj=j+1; jj<nmodels; jj++){
	py_subresult = PyFloat_FromDouble(scores[k] * max_normed);
	// py_subresult = PyFloat_FromDouble(scores[k]);
	// cout << " " << j << " "<<jj<<" "<<scores[k] << " " << scores[k] * max_normed << " " << max_normed<<endl << flush;
	PyDict_SetItem(py_result, PyTuple_Pack(2, PyList_GET_ITEM(py_models, j ), PyList_GET_ITEM(py_models, jj)), py_subresult);
	PyDict_SetItem(py_result, PyTuple_Pack(2, PyList_GET_ITEM(py_models, jj), PyList_GET_ITEM(py_models, j )), py_subresult);
	k++;
      }
    }
  }
  // cout << "START8" << endl << flush;

  // free
  // cout << "START5" << endl << flush;
  delete[] drmsds;
  // cout << "START5" << endl << flush;
  delete[] nrmsds;
// cout << "START5" << endl << flush;
  delete[] scores;
// cout << "START5" << endl << flush;
  for (int j=0; j<nmodels; j++){
    for (int i=0; i<size; i++){
      delete[] xyzn[j][i];
    }
    delete[] xyzn[j];
  }
// cout << "START5" << endl << flush;
  delete[] xyzn;
// cout << "START5" << endl << flush;
  Py_DECREF(py_xs);
// cout << "START5" << endl << flush;
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
 
static PyMethodDef Eqv_rms_drmsMethods[] =
  {
    {"rmsdRMSD_wrapper", rmsdRMSD_wrapper, METH_VARARGS, 
    rmsdRMSD_wrapper__doc__},
    {NULL, NULL, 0, NULL}
  };

PyMODINIT_FUNC
 
initeqv_rms_drms(void)
{
  (void) Py_InitModule3("eqv_rms_drms", Eqv_rms_drmsMethods, 
			"Functions to compaire two Chromatin strands.");
}
