
// gcc tadbit_py.c -I/usr/include/python2.7 -lm -lpthread
// testing:
// gcc -shared tadbit_py.c -I/usr/include/python2.7 -lm -lpthread -std=gnu99 -fPIC -g -O3 -Wall -o tadbit_py.so

#include "Python.h"
#include "tadbit.c"

/* The module doc string */
PyDoc_STRVAR(tadbit_py__doc__,
"tadbitpyc doc");

/* The function doc string */
PyDoc_STRVAR(py_tadbit__doc__,
"run tadbit in tadbit.c");


/* The wrapper to the underlying C function */
static PyObject *py_tadbit (PyObject *self, PyObject *args){
  PyObject **obs;
  int n;
  const int m;
  int n_threads;
  const int verbose;
  const int speed;
  const int heuristic;
  /* output */
  tadbit_output *seg = (tadbit_output *) malloc(sizeof(tadbit_output));

  if (!PyArg_ParseTuple(args, "Oiiiiii:tadbit", &obs, &n, &m, &n_threads, &verbose, &speed, &heuristic))
    return NULL;

  // convert list of lists to pointer o pointers
  int i;
  int j;
  PyObject * list2;
  PyObject * pyfloat;
  double **list = (double **) malloc(m * sizeof(double **));
  for (i = 0 ; i < m ; i++) {
    list2 = PyList_GetItem(list, i);
    for (j = 0 ; j < 9 ; j++){
      printf("%d\n", j);
      pyfloat = PyList_GetItem(list2, j);
      list[i][j] = PyFloat_AsDouble(pyfloat);
    }
  }
  printf("->%d matrices, %d rows\n", m,n);
  tadbit(list, n, m, n_threads, verbose, speed, heuristic, seg);
  return 1;
}

/* A list of all the methods defined by this module. */
/* "iterate_point" is the name seen inside of Python */
/* "py_iterate_point" is the name of the C function handling the Python call */
/* "METH_VARGS" tells Python how to call the handler */
/* The {NULL, NULL} entry indicates the end of the method definitions */
static PyMethodDef tadbit_py_methods[] = {
	{"py_tadbit",  py_tadbit, METH_VARARGS, py_tadbit__doc__},
	{NULL, NULL}      /* sentinel */
};

/* When Python imports a C module named 'X' it loads the module */
/* then looks for a method named "init"+X and calls it.  Hence */
/* for the module "mandelbrot" the initialization function is */
/* "initmandelbrot".  The PyMODINIT_FUNC helps with portability */
/* across operating systems and between C and C++ compilers */
PyMODINIT_FUNC
inittadbit_py(void)
{
	/* There have been several InitModule functions over time */
	Py_InitModule3("tadbit_py", tadbit_py_methods,
                   tadbit_py__doc__);
}
