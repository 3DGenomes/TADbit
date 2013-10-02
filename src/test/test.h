#define _GNU_SOURCE
#include <stdlib.h>
#include <time.h>
#include <glib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include "tadbit.h"

int _max_cache_index;
int *_cache_index;

double
ll(
  const int    n,
  const int    i_,
  const int    _i,
  const int    j_,
  const int    _j,
  const int    diag,
  const int    *k,
  const double *d,
  const double *w,
  const double *lg,
        double *c
);

