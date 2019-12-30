#include "Python.h"
#include "squared_distance_matrix_calculation.c"

/* The function doc string */
PyDoc_STRVAR(squared_distance_matrix_calculation__doc__,
"From lists of lists of xyz positions (the cartesian positions of all the particles) \n\
in the model, we get the distance matrix between particles pairs \n\
   :param model: list of tuples of (x, y, z) coordinates.\n\
   :param nparticles: number of particles to analyse\n\
\n\
   :returns: list of lists of floats of the distances between all\n\
the pairs of particles in the model \n\
");


static PyObject* squared_distance_matrix_calculation_wrapper(PyObject* self, PyObject* args)
{
  /* 
     These are definitions of Python objects needed to import the coordinates from TADbit generated models.
     In TADbit each model has the coordinates stored in a dictionary composed of 3 lists ['x'], ['y'], ['z']. 
     py_x : List of the x coordinates;
     py_y : List of the y coordinates;
     py_z : List of the z coordinates;
     nparticles : Number of particles in a single model;
  */
  PyObject *py_x;
  PyObject *py_y;
  PyObject *py_z;
  int nparticles;
  int particle, particle1, particle2;
  
  /* 
     Here we import the input parameters into the Python objects
     created above 
  */
  if (!PyArg_ParseTuple(args, "OOOi", &py_x, &py_y, &py_z, 
			&nparticles))
    return NULL;

  
  /* These are definition of C variables to perform the calculation */
  double **model;
  double **squared_distance_matrix;

  /* These are commands to allocate the memory for the C arrays */
  model = matrix2_d(nparticles, 3);  
  squared_distance_matrix = matrix2_d(nparticles, nparticles);  
  
  /* These are the commands to store the information (in this case the
     coordinates) from Python objects to C arrays */
  for(particle=0; particle < nparticles; particle++)
    {
      model[particle][0]=PyFloat_AS_DOUBLE(PyList_GET_ITEM(py_x, particle));
      model[particle][1]=PyFloat_AS_DOUBLE(PyList_GET_ITEM(py_y, particle));
      model[particle][2]=PyFloat_AS_DOUBLE(PyList_GET_ITEM(py_z, particle));
      //fprintf(stderr, "%d %f %f %f\n", particle, model[particle][0], model[particle][1], model[particle][2]);
    }
  compute_squared_distance_matrix(model, nparticles, squared_distance_matrix); 
      
  /* Store the results in a Python object to be used int TADbit */
  PyObject * py_squared_distance_matrix = NULL;
  PyObject * py_tmp_list = NULL;
  py_squared_distance_matrix = PyList_New(nparticles);
  for (particle1 = 0; particle1 < nparticles; ++particle1)
    {
      py_tmp_list = PyList_New(nparticles);
      for (particle2 = 0; particle2 < nparticles; ++particle2)
	{	      
	  PyList_SetItem(py_tmp_list, particle2, PyFloat_FromDouble(squared_distance_matrix[particle1][particle2]));
	}
      PyList_SetItem(py_squared_distance_matrix, particle1, py_tmp_list);
    }

  free2_d(model, nparticles);
  free2_d(squared_distance_matrix, nparticles);
  
  return py_squared_distance_matrix; 
}


 
static PyMethodDef squared_distance_matrix_Methods[] =
  {
    {"squared_distance_matrix_calculation_wrapper", squared_distance_matrix_calculation_wrapper, METH_VARARGS, 
     squared_distance_matrix_calculation__doc__},
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

MOD_INIT(squared_distance_matrix) {

	PyObject *m;

	MOD_DEF(m, "squared_distance_matrix", "Functions to compute the distance map of a given model.",
			squared_distance_matrix_Methods)
	if (m == NULL)
		return MOD_ERROR_VAL;

	return MOD_SUCCESS_VAL(m);

}
