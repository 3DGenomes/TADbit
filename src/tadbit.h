#include <pthread.h>

#ifndef TADBIT_LOADED
#define TADBIT_LOADED

#define TOLERANCE 1e-6
#define MAXITER 10000


// Struct to hold slice data.
typedef struct {
   int size;
   double *lgamma;
   double *counts;
   double *dist;
   double *weights;
} ml_block;


typedef struct {
   int n;
   int m;
   double **obs;
   double *dist;
   double **log_gamma;
   int *skip;
   double *llikmat;
   int verbose;
   int speed;
} thread_arg;


// 'tadbit' output struct.
typedef struct {
   int maxbreaks;
   int nbreaks_opt;
   double *llikmat;
   double *mllik;
   int *bkpts;
} tadbit_output;


void
tadbit(
  /* input */
  double **obs,
  int n,
  const int m,
  int n_threads,
  const int verbose,
  const int speed,
  const int heuristic,
  /* output */
  tadbit_output *seg
);
#endif
