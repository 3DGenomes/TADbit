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

typedef struct {
   int n;
   int m;
   double **obs;
   double *dist;
   double **log_gamma;
   double *llikmat;
   int verbose;
} thread_arg;

// 'tadbit' output struct.
typedef struct {
   int nbreaks_opt;
   double *llikmat;
   double *mllik;
   int *bkpts;
} tadbit_output;


void
free_tadbit_ouput(
  tadbit_output *seg
);

void
tadbit(
  /* input */
  double **obs,
  int n,
  const int m,
  double max_tad_size,
  int n_threads,
  const int verbose,
  const int heuristic,
  /* output */
  tadbit_output *seg
);
#endif
