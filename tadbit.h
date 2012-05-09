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
double ml_ab(const double *, const double *, double *, const int);
int *tadbit(const double **, const int, const int, const int, double,
      const int, const int);
#endif
