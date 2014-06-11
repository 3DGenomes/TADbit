
// testing:
// gcc -shared tadbit_py.c -I/usr/include/python2.7 -lm -lpthread -std=gnu99 -fPIC -g -O3 -Wall -o tadbit_py.so

#include "Python.h"
#include "tadbit.c"
#include "norm-lib/visibility.h"

/* The module doc string */
PyDoc_STRVAR(tadbit_py__doc__,
"here is a wrapper to tadbit function");

/* The function doc string */
PyDoc_STRVAR(_tadbit_wrapper__doc__,
"Run tadbit function in tadbit.c.\n\
    :argument obs: a python list of lists of int, representing a list of linearized matrices.\n\
    :argument weights: a python list of lists of floats, representing a list of linearized matrices.\n\
    :argument remove: a python list of lists of booleans mapping positively columns to remove.\n\
    :argument 0 n: number of rows or columns in the matrix\n\
    :argument 0 m: number of matrices\n\
    :argument 0 n_threads: number of threads to use\n\
    :argument 0 verbose: whether to display more/less information about process\n\
    :argument 0 max_tad_size: an integer defining maximum size of TAD. Default defines it to the number of rows/columns.\n\
    :argument 1 do_not_use_heuristic: whether to use or not some heuristics\n\
    :returns: a python list with each\n");


/* The wrapper to the underlying C function */
static PyObject *_tadbit_wrapper (PyObject *self, PyObject *args){
  PyObject **py_obs;
  PyObject *py_remove;
  int normalization;
  int n;
  int m;
  int n_threads;
  const int verbose;
  const int max_tad_size;
  const int nbks;
  const int do_not_use_heuristic;
  /* output */
  tadbit_output *seg = (tadbit_output *) malloc(sizeof(tadbit_output));

  if (!PyArg_ParseTuple(args, "OOiiiiiiii:tadbit", &py_obs, &py_remove, 
			&normalization, &n, &m, &n_threads, &verbose, 
			&max_tad_size, &nbks, &do_not_use_heuristic))
    return NULL;
  // convert list of lists to pointer o pointers
  // if something goes wrong, it is probably from there :S
  int i, j;
  int **obs;
  obs = malloc(m * sizeof(int*));
  for (i = 0 ; i < m ; i++ )
    obs[i] = malloc(n*n * sizeof(int));
  for (i = 0 ; i < m ; i++)
    for (j = 0 ; j < n*n ; j++)
      obs[i][j] = PyInt_AS_LONG(PyTuple_GET_ITEM(PyList_GET_ITEM(py_obs, i), j));

  double **weights = visibility(obs, m, n);

  char *remove = (char *) malloc (n * sizeof(char));
  for (j = 0 ; j < n ; j++){
    remove[j] = PyInt_AS_LONG(PyTuple_GET_ITEM(py_remove, j)); // automatic casting into char
  }

  // run tadbit
  tadbit(obs, weights, remove, n, m, n_threads, verbose, max_tad_size, nbks, do_not_use_heuristic, seg);

  // store each tadbit output

  // declare python objects to store lists
  PyObject * py_bkpts;
  PyObject * py_llikmat;
  PyObject * py_mllik;
  PyObject * py_result;
  PyObject * py_passages;

  // get bkpts
  const int MAXBREAKS = nbks ? nbks : n/5;
  int dim = MAXBREAKS * n;
  py_bkpts = PyList_New(dim);
  for(i = 0 ; i < dim; i++)
    PyList_SetItem(py_bkpts, i, PyInt_FromLong(seg->bkpts[i]));

  // get passages
  py_passages = PyList_New(n);
  for(i = 0 ; i < n; i++)
    PyList_SetItem(py_passages, i, PyFloat_FromDouble(seg->passages[i]));

  // get llikmat
  py_llikmat = PyList_New(n*n);
  for(i = 0 ; i < n*n; i++)
    PyList_SetItem(py_llikmat, i, PyFloat_FromDouble(seg->llikmat[i]));

  // get mllik
  py_mllik = PyList_New(seg->maxbreaks);
  for(i = 0 ; i < seg->maxbreaks ; i++)
    PyList_SetItem(py_mllik, i, PyFloat_FromDouble(seg->mllik[i]));

  // group results into a python list
  py_result = PyList_New(6);

  PyList_SetItem(py_result, 0, PyInt_FromLong(seg->maxbreaks));
  PyList_SetItem(py_result, 1, PyInt_FromLong(seg->nbreaks_opt));
  PyList_SetItem(py_result, 2, py_passages);
  PyList_SetItem(py_result, 3, py_llikmat);
  PyList_SetItem(py_result, 4, py_mllik);
  PyList_SetItem(py_result, 5, py_bkpts);

  // free many things... no leaks here!!
  for (i = 0 ; i < m ; i++){
    free(obs[i]);
    free(weights[i]);
  }
  free(obs);
  free(weights);

  destroy_tadbit_output(seg);

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
