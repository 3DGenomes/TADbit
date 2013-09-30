#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "tadbit.h"

// Declare and register R/C interface.
SEXP
tadbit_R_call(
  SEXP list,
  SEXP n_threads,
  SEXP verbose,
  SEXP max_tad_size,
  SEXP heuristic
);

R_CallMethodDef callMethods[] = {
   {"tadbit_R_call", (DL_FUNC) &tadbit_R_call, 5},
   {NULL, NULL, 0}
};

void R_init_tadbit(DllInfo *info) {
   R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}


SEXP
tadbit_R_call(
  SEXP list,
  SEXP n_threads,
  SEXP verbose,
  SEXP max_tad_size,
  SEXP do_not_use_heuristic
){

/*
   * This is a tadbit wrapper for R. The matrices have to be passed
   * in a list (in R). Checks that the input consists of numeric
   * square matrices, with identical dimensions. The list is
   * is converted to pointer of pointers to doubles and passed
   * to 'tadbit'.
   * Assume that NAs can be passed from R and are ignored in the
   * computation.
*/

   R_len_t i, m = length(list);
   int first = 1, N, *dim_int;

   SEXP dim;
   PROTECT(dim = allocVector(INTSXP, 2));

   // Convert 'obs_list' to pointer of pointer to double.
   int **obs = (int **) malloc(m * sizeof(int **));
   for (i = 0 ; i < m ; i++) {
      // This fails if list element is not numeric.
      obs[i] = INTEGER(coerceVector(VECTOR_ELT(list, i), INTSXP));
      // Check that input is a matrix.
      if (!isMatrix(VECTOR_ELT(list, i))) {
         error("input must be square matrix");
      }
      // Check the dimension.
      dim = getAttrib(VECTOR_ELT(list, i), R_DimSymbol);
      dim_int = INTEGER(dim);
      if (dim_int[0] != dim_int[1]) {
         error("input must be square matrix");
      }
      if (first) {
         N = dim_int[0];
         first = 0;
      }
      else {
         if (N != dim_int[0]) {
            error("all matrices must have same dimensions");
         }
      }
   }

   UNPROTECT(1);

   tadbit_output *seg = (tadbit_output *) malloc(sizeof(tadbit_output));
   
   // Call 'tadbit'.
   tadbit(obs, N, m, INTEGER(n_threads)[0], INTEGER(verbose)[0],
         INTEGER(max_tad_size)[0], INTEGER(do_not_use_heuristic)[0], seg);

   int maxbreaks = seg->maxbreaks;

   // Check that tadbit exited successfully.
   if (maxbreaks == -1) {
      free(obs);
      free(seg);
      SEXP fail_SEXP;
      PROTECT(fail_SEXP = allocVector(INTSXP, 1));
      int *fail = INTEGER(fail_SEXP);
      fail[0] = -1;
      SEXP list_SEXP;
      PROTECT(list_SEXP = allocVector(VECSXP, 1));
      SET_VECTOR_ELT(list_SEXP, 0, fail_SEXP);
      UNPROTECT(2);
      return list_SEXP;
   }

   // Copy output to R-readable variables.
   SEXP nbreaks_SEXP;
   SEXP passages_SEXP;
   SEXP llikmat_SEXP;
   SEXP mllik_SEXP;
   SEXP bkpts_SEXP;

   PROTECT(nbreaks_SEXP = allocVector(INTSXP, 1));
   PROTECT(passages_SEXP = allocVector(INTSXP, N));
   PROTECT(llikmat_SEXP = allocVector(REALSXP, N*N));
   PROTECT(mllik_SEXP = allocVector(REALSXP, maxbreaks));
   PROTECT(bkpts_SEXP = allocVector(INTSXP, N*(maxbreaks-1)));

   int *nbreaks_opt = INTEGER(nbreaks_SEXP); 
   int *passages = INTEGER(passages_SEXP); 
   double *llikmat = REAL(llikmat_SEXP);
   double *mllik = REAL(mllik_SEXP);
   int *bkpts = INTEGER(bkpts_SEXP);

   nbreaks_opt[0] = seg->nbreaks_opt;
   for (i = 0 ; i < N ; i++) passages[i] = seg->passages[i];
   for (i = 0 ; i < N*N ; i++) llikmat[i] = seg->llikmat[i];
   for (i = 0 ; i < maxbreaks ; i++) mllik[i] = seg->mllik[i];
   // Remove first column associated with 0 breaks. Itcontains only
   // 0s and shifts the index in R (vectors start at position 1).
   for (i = N ; i < N*(maxbreaks-1) ; i++) bkpts[i-N] = seg->bkpts[i];


   // Set 'dim' attributes.
   SEXP dim_llikmat;
   PROTECT(dim_llikmat = allocVector(INTSXP, 2));
   INTEGER(dim_llikmat)[0] = N;
   INTEGER(dim_llikmat)[1] = N;
   setAttrib(llikmat_SEXP, R_DimSymbol, dim_llikmat);

   SEXP dim_breaks;
   PROTECT(dim_breaks = allocVector(INTSXP, 2));
   INTEGER(dim_breaks)[0] = N;
   INTEGER(dim_breaks)[1] = maxbreaks-1;
   setAttrib(bkpts_SEXP, R_DimSymbol, dim_breaks);

   free(obs);
   free(seg->passages);
   free(seg->llikmat);
   free(seg->mllik);
   free(seg->bkpts);
   free(seg);

   SEXP list_SEXP;
   PROTECT(list_SEXP = allocVector(VECSXP, 5));
   SET_VECTOR_ELT(list_SEXP, 0, nbreaks_SEXP);
   SET_VECTOR_ELT(list_SEXP, 1, llikmat_SEXP);
   SET_VECTOR_ELT(list_SEXP, 2, mllik_SEXP);
   SET_VECTOR_ELT(list_SEXP, 3, bkpts_SEXP);
   SET_VECTOR_ELT(list_SEXP, 4, passages_SEXP);
   UNPROTECT(8);

   return list_SEXP;

}
