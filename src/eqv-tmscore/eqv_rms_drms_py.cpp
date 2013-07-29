#include "Python.h"
#include "eqv-drmsd-tadbit.cpp"

/* The function doc string */
PyDoc_STRVAR(rmsdRMSD_wrapper__doc__,
"From 2 lists of xyz positions, and a given threshold (nm), return the \n\
number of equivalent positions, the RMSD and the dRMSD.\n\
   :param xyz_a: first list of tuples of (x, y, z ) coordinates.\n\
   :param xyz_b: second list of tuples of (x, y, z ) coordinates.\n\
   :param dcutoff: distance cutoff to consider 2 particles as equivalent \n\
      in position (nm)\n\
   :param consistency: whether to get consistency per particle\n\
\n\
   :returns: a list with the number of equivalent positions, the RMSD and \n\
      the dRMSD, of the alignment. If consistency is True, returns list of \n\
      0 or 1 if a given particle is in equivalent position in both strands.\n\
");


static PyObject* rmsdRMSD_wrapper(PyObject* self, PyObject* args)
{
  PyObject **py_xyzA;
  PyObject **py_xyzB;
  int size;
  float thres;
  int consistency;
 
  if (!PyArg_ParseTuple(args, "OOifi", &py_xyzA, &py_xyzB, &size, &thres, &consistency))
    return NULL;
 
  float **list1 = NULL;
  float **list2 = NULL;
  int i, j;
  list1 = (float**)malloc(size * sizeof(float*));
  for (i = 0 ; i < size ; i++ )
    list1[i] = (float*)malloc(3 * sizeof(float));
  list2 = (float**)malloc(size * sizeof(float*));
  for (i = 0 ; i < size ; i++ )
    list2[i] = (float*)malloc(3 * sizeof(float));

  for (i=0; i<size; i++){
    for (j=0; j<=2; j++){
      list1[i][j] = PyFloat_AS_DOUBLE(PyTuple_GET_ITEM(PyList_GET_ITEM(py_xyzA, i), j));
      list2[i][j] = PyFloat_AS_DOUBLE(PyTuple_GET_ITEM(PyList_GET_ITEM(py_xyzB, i), j));
    }
  }
  float rms, drms;

  int *cons_list;
  cons_list = new int[size];

  int eqv;
  rmsdRMSD(list1, list2, size, thres, eqv, rms, drms, cons_list, consistency);

  PyObject * py_result = NULL;
  if (consistency){
    py_result = PyList_New(size);
    for (i=0; i<size; i++)
      PyList_SetItem(py_result, i, PyInt_FromLong(cons_list[i]));
  }else{
    py_result = PyList_New(3);
    PyList_SetItem(py_result, 0, PyInt_FromLong(eqv));
    PyList_SetItem(py_result, 1, PyFloat_FromDouble(rms));
    PyList_SetItem(py_result, 2, PyFloat_FromDouble(drms));
  }

  free(cons_list);
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
			"Functions to compaire two Chromatine strands.");
}
