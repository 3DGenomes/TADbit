
// testing:
// gcc -shared tadbit_py.c -I/usr/include/python2.7 -lm -lpthread -std=gnu99 -fPIC -g -O3 -Wall -o tadbit_py.so

#include "Python.h"
#include "tadbit.c"

/* The module doc string */
PyDoc_STRVAR(tadbit_py__doc__,
"here is a wrapper to tadbit function");

/* The function doc string */
PyDoc_STRVAR(_tadbit_wrapper__doc__,
"Run tadbit function in tadbit.c.\n\
    :argument obs: a python list of lists of floats, representing a list of linearized matrices.\n\
    :argument 0 n: number of rows or columns in the matrix\n\
    :argument 0 m: number of matrices\n\
    :argument 0 n_threads: number of threads to use\n\
    :argument 0 verbose: wether to display more/less information about process\n\
    :argument 0 speed: can be 0, 2, 3 or 4. Divide the calculated area by 2**(speed-1)\n\
    :argument 0 heuristic: whether to use or not some heuristics\n\
    :returns: a python list with each\n");


/* The wrapper to the underlying C function */
static PyObject *_tadbit_wrapper (PyObject *self, PyObject *args){
  PyObject **obs;
  int n=0;
  int m=0;
  int n_threads=0;
  const int verbose=0;
  const int speed=0;
  const int heuristic=0;
  /* output */
  tadbit_output *seg = (tadbit_output *) malloc(sizeof(tadbit_output));

  if (!PyArg_ParseTuple(args, "Oiiiiii:tadbit", &obs, &n, &m, &n_threads, &verbose, &speed, &heuristic))
    return NULL;

  // convert list of lists to pointer o pointers
  // if something goes wrong, it is probably from there :S
  int i;
  int j;
  double ** list;
  list = malloc(m * sizeof(double **));
  for (i = 0 ; i < m ; i++ )
    list[i] = malloc(n*n * sizeof(double*));
  for (i = 0 ; i < m ; i++){
    printf("\ni\n");
    PyObject * tmplist = PyList_GET_ITEM(obs, i);
    for (j = 0 ; j < n*n ; j++){
      printf(" j\n");
      list[i][j] =  PyFloat_AsDouble(PyTuple_GET_ITEM(tmplist, j));
      printf(" %f", list[i][j]);
    }
  }

  // run tadbit
  tadbit(list, n, m, n_threads, verbose, speed, heuristic, seg);

  // store each tadbit output
  int mbreaks      = seg->maxbreaks;
  int nbreaks_opt  = seg->nbreaks_opt;
  double * llikmat = seg->llikmat;
  double * mllik   = seg->mllik;
  int    * bkpts   = seg->bkpts;

  // declare python objects to store lists
  PyObject * py_bkpts;
  PyObject * py_llikmat;
  PyObject * py_mllik;
  PyObject * py_result;

  // get bkpts
  int dim = nbreaks_opt*n;
  py_bkpts = PyList_New(dim+n);
  for(i = 0 ; i < dim+n; i++){
    PyList_SetItem(py_bkpts, i, PyInt_FromLong(bkpts[i]));
  }
  /* This is to return directly the list of breaks found
  j = 0;
  py_bkpts = PyList_New(nbreaks_opt);
  for(i = 0 ; i < n; i++){
    if (bkpts[i+dim]==1){
      PyList_SetItem(py_bkpts, j, PyInt_FromLong(i));
      j++;
    }
  }
  */

  // get llikmat
  py_llikmat = PyList_New(n*n+n);
  for(i = 0 ; i < n*n+n; i++){
    PyList_SetItem(py_llikmat, i, PyFloat_FromDouble(llikmat[i]));
  }

  // get mllik
  py_mllik = PyList_New(mbreaks);
  for(i = 0 ; i < mbreaks ; i++){
    PyList_SetItem(py_mllik, i, PyFloat_FromDouble(mllik[i]));
  }

  // group results into a python list
  py_result = PyList_New(5);

  PyList_SetItem(py_result, 0, PyInt_FromLong(mbreaks));
  PyList_SetItem(py_result, 1, PyInt_FromLong(nbreaks_opt));
  PyList_SetItem(py_result, 2, py_llikmat);
  PyList_SetItem(py_result, 3, py_mllik);
  PyList_SetItem(py_result, 4, py_bkpts);

  return py_result;
}

/* A list of all the methods defined by this module. */
/* The {NULL, NULL} entry indicates the end of the method definitions */
static PyMethodDef tadbit_py_methods[] = {
	{"_tadbit_wrapper",  _tadbit_wrapper, METH_VARARGS, _tadbit_wrapper__doc__},
	{NULL, NULL}      /* sentinel */
};

/* When Python imports a C module named 'X' it loads the module */
/* then looks for a method named "init"+X and calls it.  Hence */
/* for the module "mandelbrot" the initialization function is */
/* across operating systems and between C and C++ compilers */
PyMODINIT_FUNC
inittadbit_py(void)
{
	/* There have been several InitModule functions over time */
	Py_InitModule3("tadbit_py", tadbit_py_methods,
                   tadbit_py__doc__);
}
