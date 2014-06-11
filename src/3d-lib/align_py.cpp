#include "Python.h"
#include "align.h"


/* The function doc string */
PyDoc_STRVAR(aligner3d__doc__,
"From a list of lists of xyz positions,  the \n\
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


static PyObject* aligner3d_wrapper(PyObject* self, PyObject* args)
{
  PyObject *py_xs1;
  PyObject *py_ys1;
  PyObject *py_zs1;
  PyObject *py_xs2;
  PyObject *py_ys2;
  PyObject *py_zs2;
  int size;

  if (!PyArg_ParseTuple(args, "OOOOOOi", &py_xs1, &py_ys1, &py_zs1, 
			&py_xs2, &py_ys2, &py_zs2, &size))
    return NULL;
 
  float **xyz1;
  float **xyz2;
  int i;

 
  xyz1 = new float*[size];
  xyz2 = new float*[size];
  for(int i=0; i<size; i++) {
    xyz1[i] = new float[3];
    memset(xyz1[i], 0, 3*sizeof(float));
    xyz2[i] = new float[3];
    memset(xyz2[i], 0, 3*sizeof(float));
  }


  for (i=0; i<size; i++){
    xyz1[i][0] = PyFloat_AS_DOUBLE(PyList_GET_ITEM(py_xs1, i));
    xyz1[i][1] = PyFloat_AS_DOUBLE(PyList_GET_ITEM(py_ys1, i));
    xyz1[i][2] = PyFloat_AS_DOUBLE(PyList_GET_ITEM(py_zs1, i));
    xyz2[i][0] = PyFloat_AS_DOUBLE(PyList_GET_ITEM(py_xs2, i));
    xyz2[i][1] = PyFloat_AS_DOUBLE(PyList_GET_ITEM(py_ys2, i));
    xyz2[i][2] = PyFloat_AS_DOUBLE(PyList_GET_ITEM(py_zs2, i));
  }

  align(xyz2, xyz1, size);

  // give it to me
  PyObject * py_result = NULL;
  PyObject * py_subresult = NULL;
  py_result = PyList_New(3);
  for (int j = 0; j < 3; ++j) {
    py_subresult = PyList_New(size);
    for (int i = 0; i < size; ++i) {
      PyList_SetItem(py_subresult, i, PyFloat_FromDouble(xyz2[i][j]));
    }
    PyList_SetItem(py_result, j, py_subresult);
  }
  for (int i=0; i<size; i++) {
    delete[] xyz1[i];
    delete[] xyz2[i];
  }
  delete[] xyz1;
  delete[] xyz2;
  
  return py_result;
}


 
static PyMethodDef aligner3dMethods[] =
  {
    {"aligner3d_wrapper", aligner3d_wrapper, METH_VARARGS, 
    aligner3d__doc__},
    {NULL, NULL, 0, NULL}
  };

PyMODINIT_FUNC
 
initaligner3d(void)
{
  (void) Py_InitModule3("aligner3d", aligner3dMethods, 
			"Functions to align of a given group of models.");
}
