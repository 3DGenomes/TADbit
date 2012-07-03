#ifndef TADBIT_LOADED
#define TADBIT_LOADED
#define TOLERANCE 1e-6
#define MAXITER 10000


// Struct to hold slice data.
typedef struct {
   int size[3];
   double *k[3];
   double *d[3];
   double *w[3];
} ml_slice;

// Struct to hold block data.
typedef struct {
   // Diagonal block data.
   int n_diag_blk;
   int *sizes_diag_blk;
   double **k_diag;
   double **d_diag;
   // Off-diagonal block data.
   int n_off_blk;
   int *sizes_off_blk;
   double **k_off;
   double **d_off;
} ml_block;

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
