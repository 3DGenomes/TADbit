#ifndef TADBIT_LOADED
#define TADBIT_LOADED
#define TOLERANCE 1e-6

// A struct to hold block data.
typedef struct {
   int size[3];
   double *k[3];
   double *d[3];
} ml_blocks;


// A struct for thread argument.
typedef struct {
   int m;
   int n;
   double *llik;
   double *dis;
   double **obs;
   int *bkpts;
   int *done;
} thread_arg;


int *tadbit(double **, int, int, int, int);
#endif
