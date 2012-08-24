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
   ml_block *blocks[3];
} ml_slice;


void
tadbit(
  double **obs,
  const int n,
  const int m,
  double max_tad_size,
  int n_threads,
  const int verbose,
  int *breaks,
  double *llik,
  int heuristic
);
#endif
