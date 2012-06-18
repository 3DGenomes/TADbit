#ifndef TADBIT
#define TADBIT
#define TOLERANCE 1e-6
#define MAXITER 10000


// Struct to hold slice data.
typedef struct {
   int size[3];
   double *k[3];
   double *d[3];
   double *w[3];
} ml_slice;


void
tadbit(
  double **obs,
  const int n,
  const int m,
  double max_tad_size,
  int n_threads,
  const int verbose,
  int *return_val
);
#endif
