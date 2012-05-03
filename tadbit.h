#ifndef TADBIT_LOADED
#define TADBIT_LOADED
#define TOLERANCE 1e-6

// A struct to hold block data.
typedef struct {
   int size[3];
   double *k[3];
   double *d[3];
} ml_blocks;

int *tadbit(double **obs, int n, int m, int fast);
#endif
