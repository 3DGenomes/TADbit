/* @(#)iterative.h
 */
/*
 *  Created by Ivan Junier
 *  
 */

#include <sys/types.h>
#include <unistd.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <pstream.h>
#include <cfloat>
#include <climits>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>
#include <boost/algorithm/string.hpp>
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

