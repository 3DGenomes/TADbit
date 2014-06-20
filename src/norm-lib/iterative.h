/* @(#)iterative.h
 */
/*
 *  Created by Ivan Junier
 *  
 */

#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include <boost/regex.hpp>
#include <boost/regex/v4/regex.hpp>


#ifndef _ITERATIVE_H
#define _ITERATIVE_H 1

#ifdef __cplusplus
    extern "C" {
#endif
      double **iterative(int **obs, int m, int n);
#ifdef __cplusplus
    }
#endif

#endif  /* _ITERATIVE_H */

