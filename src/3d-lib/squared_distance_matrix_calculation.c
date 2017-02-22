#include<stdio.h>
#include<stdlib.h>

void compute_squared_distance_matrix(double **model, int nparticles, double **squared_distance_matrix);
double squared_distance_between_particles(double *P1, double *P2);
double *vector1_d(int dim);
double **matrix2_d(int dim1, int dim2);
void free1_d(double *v);
void free2_d(double **m, int dim1);
void free3_d(double ***t, int dim1, int dim2);

/*************************/

void compute_squared_distance_matrix(double **model, int nparticles, double **squared_distance_matrix)
{
  int particle1, particle2;    /* particles, and components indeces */

  for(particle1 = 0; particle1 < nparticles; particle1++)
    {
      for(particle2 = particle1+1; particle2 < nparticles; particle2++)
	{
	  squared_distance_matrix[particle1][particle2] = squared_distance_between_particles(model[particle1], model[particle2]);
	  squared_distance_matrix[particle2][particle1] = squared_distance_matrix[particle1][particle2];
	}
    }
}

/*************************/

double squared_distance_between_particles(double *P1, double *P2)
{
  int   component;
  double squared_distance;

  squared_distance = 0.0;
  for(component = 0; component < 3; component++)
    {
      squared_distance += (P1[component]-P2[component])*(P1[component]-P2[component]);
    }
  return(squared_distance);
}

/*************************/

/* Memory allocation and initialization to zero for an double vector of size dim */
double *vector1_d(int dim)
{
  /* MDS 2011 */
  int i;
  double *v;

  if((v = (double *)malloc((size_t) dim * sizeof(double))) == NULL)
    {
      fprintf(stderr, "Memory allocation problems");
      exit(1);
    }

  for(i = 0; i < dim; ++i)
    {
      v[i] = 0.0;
    }
  return(v);
}

/*************************/

/* Memory allocation for an double matrix of size dim1 * dim2 */
double **matrix2_d(int dim1, int dim2)
{
  /* MDS 2011 */
  int i;
  double **m;

  m = (double **)malloc(dim1 * sizeof(double *));
  if(m == NULL)
    {
      fprintf(stderr, "Memory allocation problems");
      exit(1);
    }

  for(i = 0; i < dim1; ++i)
    {
      m[i] = vector1_d(dim2);
    }

  return(m);
}

/*******************************/

/* Memory allocation for an double tensor of size dim1 * dim2 * dim3 */
double ***tensor3_d(int dim1, int dim2, int dim3)
{
  /* MDS 2011 */
  int i;
  double ***t;

  t = (double ***)malloc(dim1 * sizeof(double **));
  if(t == NULL)
    {
      fprintf(stderr, "Memory allocation problems");
      exit(1);
    }

  for(i = 0; i < dim1; ++i)
    {
      t[i] = matrix2_d(dim2, dim3);
    }

  return(t);
}

/*******************************/

void free1_d(double *v)
{
  free(v);
}

/*******************************/

void free2_d(double **m, int dim1)
{
  int i;

  for(i = 0; i < dim1; ++i)
    {
      free1_d(m[i]);
    }

  free(m);
}

/*******************************/

void free3_d(double ***t, int dim1, int dim2)
{
  int i;

  for(i = 0; i < dim1; ++i)
    {
      free2_d(t[i], dim2);
    }
  free(t);
}
