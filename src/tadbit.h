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
   double **k;
   double *d;
   double **w;
   double **lg;
   int *skip;
   double *llikmat;
   int verbose;
} thread_arg;


// 'tadbit' output struct.
typedef struct {
   int maxbreaks;
   int nbreaks_opt;
   int *passages;
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
  //const int speed,
  const int max_tad_size,
  const int do_not_use_heuristic,
  /* output */
  tadbit_output *seg
);
#endif
