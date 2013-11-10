#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <pthread.h>
#include <ctype.h>
#include <unistd.h>
#include <float.h>
#include <assert.h>
#include <pthread.h>

#ifndef _TADBIT_LOADED
#define _TADBIT_LOADED

#define TOLERANCE 1e-6
#define MAXITER 10000

typedef struct {
   const int n;
   const int m;
   const int **k;
   const double *d;
   const double **w;
   const double **lg;
   const char *skip;
   double *llikmat;
   const int verbose;
} llworker_arg;

typedef struct {
   const int n;
   const double *llikmat;
   double *old_llik;
   double *new_llik;
   int nbreaks;
   int *new_bkpt_list;
   const int *old_bkpt_list;
} dpworker_arg;



// 'tadbit' output struct.
typedef struct {
   int m;
   int maxbreaks;
   int nbreaks_opt;
   int *passages;
   double *llikmat;
   double *mllik;
   int *bkpts;
   double **weights;
} tadbit_output;



void
tadbit(
  /* input */
  int **obs,
  int n,
  const int m,
  int n_threads,
  const int verbose,
  //const int speed,
  const int max_tad_size,
  const int do_not_use_heuristic,
  /* output */
  tadbit_output *seg
);


void
destroy_tadbit_output(
   tadbit_output *seg
);
#endif
