#ifndef TADBIT_LOADED
#define TADBIT_LOADED
#define TOLERANCE 1e-6

// A struct to hold block data.
typedef struct {
   int size[3];
   double *k[3];
   double *d[3];
} ml_blocks;

// Useful functions.
int n_proc(void);
double ml_ab(double *, double *, double *, int);
int *tadbit(double **, int, int, int, int);
#endif
