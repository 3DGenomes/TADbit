#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include "tadbit.h"

// Declare and register R/C interface.
SEXP tadbit_R_call(SEXP list, SEXP fast_yn, SEXP threads);
R_CallMethodDef callMethods[] = {
   {"tadbit_R_call", (DL_FUNC) &tadbit_R_call, 2},
   {NULL, NULL, 0}
};

void R_init_tadbit(DllInfo *info) {
   R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}


SEXP tadbit_R_call(SEXP list, SEXP fast_yn, SEXP threads) {
/*
   * This is a tadbit wrapper for R. The matrices have to passed
   * in a list (in R). Checks that the input consists of numeric
   * square matrices, with identical dimensions. The list is
   * is converted to pointer of pointers to doubles and passed
   * to 'tadbit'.
   * Assume that NAs can be passed from R and are ignored in the
   * computation.
*/

   R_len_t i, m = length(list);
   int first = 1, n, *dim;
   int fast = INTEGER(fast_yn)[0], n_threads = INTEGER(threads)[0];

   // Convert 'obs_list' to pointer of pointer to double.
   double **obs = (double **) malloc(m * sizeof(double **));
   for (i = 0 ; i < m ; i++) {
      // This fails is list element is not numeric.
      obs[i] = REAL(coerceVector(VECTOR_ELT(list, i), REALSXP));
      // Check that input is a matrix.
      if (!isMatrix(VECTOR_ELT(list, i))) {
         error("input must be square matrix");
      }
      // Check the dimension.
      dim = INTEGER(getAttrib(VECTOR_ELT(list, i), R_DimSymbol));
      if (dim[0] != dim[1]) {
         error("input must be square matrix");
      }
      if (first) {
         n = dim[0];
         first = 0;
      }
      else {
         if (n != dim[0]) {
            error("all matrices must have same dimensions");;
         }
      }
   }

   // If 'n_threads' is 0, allocate thread number automatically.
   if (!n_threads) {
      // Keep one core free and allocate one thread per core.
      n_threads = n_proc() ? n_proc() - 1 : 1;
   }

   // Call 'tadbit'.
   int *bkpts = tadbit((const double **) obs, n, m, fast, n_threads);

   // Wrap it up.
   SEXP return_val_sexp;
   PROTECT(return_val_sexp = allocVector(INTSXP, n));
   int *return_val = INTEGER(return_val_sexp);
   // Copy output from 'tadbit'.:60

   for (i = 0 ; i < n ; i++) {
      return_val[i] = bkpts[i];
   }
   free (bkpts);
   UNPROTECT(1);

   return return_val_sexp;

}

