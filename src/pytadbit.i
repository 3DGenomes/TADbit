%module pytadbit

%{
#define SWIG_FILE_WITH_INIT
#include "tadbit.h"
%}

void tadbit(double **obs, const int n, const int m,
  double max_tad_size, int n_threads, const int verbose, int *return_val);
