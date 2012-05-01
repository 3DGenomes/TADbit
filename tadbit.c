#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#define TOLERANCE 1e-6

// Declare and register R/C interface.
SEXP tadbit_R_call(SEXP list, SEXP fast_yn);
R_CallMethodDef callMethods[] = {
   {"tadbit_R_call", (DL_FUNC) &tadbit_R_call, 2},
   {NULL, NULL, 0}
};

void R_init_tadbit(DllInfo *info) {
   R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}


double ml_ab(double *k_orig, double *d_orig, double *ab, int n) {
/*
   * The 2-array 'ab' is upated in place and the log-likelihood
   * is returned.
   * The fitted model (by maximum likelihood) is Poisson with lambda
   * paramter such that lambda = exp(a + b*d). So the full log-likelihood
   * of the model is Sigma -exp(a + b*d_i) + k_i(a + b*d_i).
*/

   if (n < 1) {
      return 0.0;
   }

   int i, j = 0;
   double a = ab[0], b = ab[1], llik;
   double da = 0.0, db = 0.0, oldgrad;
   double f, g, dfda = 0.0, dfdb = 0.0, dgda = 0.0, dgdb = 0.0;
   double denom, tmp; // 'tmp' is used as computation intermediate.

   double k[n], d[n];

   // Filter out NAs while copying.
   for (i = 0 ; i < n ; i++) {
      if (isnan(d_orig[i]) || isnan(k_orig[i])) {
         continue;
      }
      else {
         d[j] = d_orig[i];
         k[j] = k_orig[i];
         j++;
      }
   }

   if (j < 3) {
      return NAN;
   }
   else {
      n = j;
   }


   // Comodity function.
   void recompute_fg() {
      f = 0.0; g = 0.0;
      for (i = 0 ; i < n ; i++) {
         tmp  =  exp(a+da+(b+db)*d[i])-k[i];
         f   +=  tmp;
         g   +=  tmp * d[i];
      }
   }

   recompute_fg();

   // Newton-Raphson until gradient function < TOLERANCE.
   while ((oldgrad = f*f + g*g) > TOLERANCE) {

      // Compute the derivatives.
      dfda = dfdb = dgda = dgdb = 0.0;
      for (i = 0 ; i < n ; i++) {
         tmp   =   exp(a+b*d[i]);
         dfda +=   tmp;
         dgda +=   tmp * d[i];
         dgdb +=   tmp * d[i]*d[i];
      }
      dfdb = dgda;

      denom = dfdb*dgda - dfda*dgdb;
      da = (f*dgdb - g*dfdb) / denom;
      db = (g*dfda - f*dgda) / denom;

      recompute_fg();
      // Gradient test. Traceback if not going down the gradient.
      while (f*f + g*g > oldgrad) {
         da /= 2;
         db /= 2;
         recompute_fg();
      }

      // Update 'a' and 'b'.
      a += da;
      b += db;


   }

   // Compute log-likelihood (using 'dfda').
   llik = 0.0;
   for (i = 0 ; i < n ; i++) {
      llik += exp(a+b*d[i]) + k[i] * (a + b*d[i]);
   }

   // Update 'ab' in place (to make the estimates available).
   ab[0] = a; ab[1] = b;

   return llik;

}

double **break_in_blocks(double *mat, int n, int i, int j, double **blocks) {
/*
   *  Break up 'mat' in three blocks delimited by 'i' and 'j'.
   *  The upper block is (0,i-1)x(i,j), the triangular block is
   *  the upper triangular block without diagonal (i,j)x(i,j)
   *  and the bottom block is (j+1,n)x(i,j).
*/

   int row, col;
   int top_counter = 0, tri_counter = 0, bot_counter = 0;

   double *top = blocks[0];
   double *tri = blocks[1];
   double *bot = blocks[2];


   // Fill vertically.
   for (col = i ; col < j+1 ; col++) {
      // Skip if 'i' is 0.
      for (row = 0 ; row < i ; row++) {
         top[top_counter++] = mat[row+col*n];
      }

      // Skip if 'col' is i.
      for (row = i ; row < col ; row++) {
         tri[tri_counter++] = mat[row+col*n];
      }

      // Skip if 'j' is n-1.
      for (row = j+1 ; row < n ; row++) {
         bot[bot_counter++] = mat[row+col*n];
      }

   }

   return blocks;

}


void remove_non_local_maxima (double **obs, double *dis, int n, int m,
    double **k_blk, double **d_blk, int *bkpts) {
/*
   * Segment the data with one breakpoint, compute the
   * likelihood and return the local maxima of that function.
   * Final breakpoints are likely to be one of those local
   * maxima.
*/

   int i, j, k;
   double llik[n];
   double ab[3][2] = {{0.0,0.0}, {0.0,0.0}, {0.0,0.0}};

/*
   * The code is not very DRY, those lines are almost duplicated
   * in the function 'tadbit'.
*/

   // Compute the log-lik of the first segment forward...
   for (j = 2 ; j < n-3 ; j++) {
      // Cut the (i,j)-blocks.
      d_blk = break_in_blocks(dis, n, 0, j, d_blk);

      llik[j] = 0.0;
      for (k = 0 ; k < m ; k++) {
         k_blk = break_in_blocks(obs[k], n, 0, j, k_blk);

         // Compute the likelihood and sum.
         llik[j] +=
             ml_ab(k_blk[1], d_blk[1], ab[1], j*(j+1)/2)          +
             ml_ab(k_blk[2], d_blk[2], ab[2], (n-j-1)*(j+1)) / 2;
      }
   }

   // START DEBUG
   printf("done forward\n");
   // END DEBUG

   // ... and the second segment backward.
   for (j = 3 ; j < n-2 ; j++) {
      // Cut the (i,j)-blocks.
      d_blk = break_in_blocks(dis, n, j, n-1, d_blk);

      for (k = 0 ; k < m ; k++) {
         k_blk = break_in_blocks(obs[k], n, j, n-1, k_blk);
         // Compute the likelihood and sum.
         llik[j-1] += 
             ml_ab(k_blk[0], d_blk[0], ab[0], (j-1)*(n-j))  / 2  +
             ml_ab(k_blk[1], d_blk[1], ab[1], (n-j-1)*(n-j)/2);
      }
   }

   bkpts[n-1] = 1;
   for (i = 0 ; i < n-1 ; i++) {
      bkpts[i] = 0;
   }
   for (i = 3 ; i < n-1 ; i++) {
      if (llik[i] > llik[i-1] && llik[i] > llik[i+1]) {
         bkpts[i] = 1;
      }
   }

}

int *get_breakpoints(double *llik, int n, int *all_breakpoints) {
/*
   * Find the maximum likelihood estimate by a dynamic algorithm.
*/

   int i,j, nbreaks = 0;
   int new_breakpoint = 0;

   double tmp;

   double new_llik[n];
   double old_llik[n];
   // Initialize to first line of 'llik'.
   for (j = 0 ; j < n ; j++) {
      old_llik[j] = llik[j*n];
   }

   // Breakpoint lists. The first index is the end of the segment,
   // the second is 1 if this position is an end (breakpoint).
   
   // int new_bkpt_list[n][n];
   // int old_bkpt_list[n][n];
   int *new_bkpt_list = (int *) malloc(n*n * sizeof(int));
   int *old_bkpt_list = (int *) malloc(n*n * sizeof(int));

   // Initialize to 0.
   for (i = 0 ; i < n*n ; i++) {
      new_bkpt_list[i] = old_bkpt_list[i] = 0;
   }

   double new_full_llik = old_llik[n-1];
   double old_full_llik = -INFINITY;

   while (old_full_llik < new_full_llik) {

      // Update breakpoints.
      nbreaks++;
      for (i = 0 ; i < n ; i++) {
         for (j = 0 ; j < n ; j++) {
            old_bkpt_list[i+j*n] = new_bkpt_list[i+j*n];
         }
      }

      // Cycle over end point 'j'.
      for (j = 3 * nbreaks + 2 ; j < n ; j++) {
         new_llik[j] = -INFINITY;

         // Cycle over start point 'i'.
         for (i = 3 * nbreaks ; i < j - 1 ; i++) {

            // NAN if not a breakpoint, so next line evaluates to false.
            tmp = old_llik[i-1] + llik[i+j*n];
            if (tmp > new_llik[j]) {
               new_llik[j] = tmp;
               new_breakpoint = i-1;
            }
         }

         // Update breakpoint list.
         if (new_llik[j] > -INFINITY) {
            for (i = 0 ; i < n ; i++) {
               new_bkpt_list[j+i*n] = old_bkpt_list[new_breakpoint+i*n];
            }
            new_bkpt_list[j+new_breakpoint*n] = 1;
         }

      }

      // Update full log-likelihoods.
      old_full_llik = new_full_llik;
      new_full_llik = new_llik[n-1];
      for (i = 0 ; i < n ; i++) {
         old_llik[i] = new_llik[i];
      }

   }
   

   for (i = 0 ; i < n ; i++) {
      all_breakpoints[i] = old_bkpt_list[n-1+i*n];
   }

   free(new_bkpt_list);
   free(old_bkpt_list);

   return all_breakpoints;

}


int *tadbit(double **obs, int n, int m, int fast) {

   int i, j, k;

/*
   * Allocate memory and initialize variables. The distance
   * matrix 'dis' is the distance to the main diagonal. Every
   * element of coordinate (i,j) is on a diagonal; the distance
   * is the shift to the main diagonal |i-j|.
   * 'd_blk' and 'k_blk' will hold the distance data ('d_blk')
   * and observation data ('k_blk') when the matrices are
   * segmented. Each segmentation defines 3 regions, which is
   * why there are 3 such matrices. They are allocated the maximum
   * size they can have upon segmentation for simplicity.
   * 'ab' contains parameters 'a' and 'b' for the maximum likelihood
   * model. Because each segmentation defines 3 regions we need
   * 3 pairs of parameters.
*/

   double *dis = (double *) malloc(n*n * sizeof(double));
   double *llik = (double *) malloc(n*n * sizeof(double));

   k = 0;
   for (i = 0; i < n ; i++) {
      for (j = 0; j < n ; j++) {
         llik[k] = NAN;
         dis[k] = abs(i-j);
         k++;
      }
   }

   // Allocate max possible size to blocks matrices.
   double **d_blk = (double **) malloc(3 * sizeof(double *));
   d_blk[0] = (double *) malloc((n+1)*(n+1)/4 * sizeof(double));
   d_blk[1] = (double *) malloc(n*(n+1)/2 * sizeof(double));
   d_blk[2] = (double *) malloc((n+1)*(n+1)/4 * sizeof(double));

   double **k_blk = (double **) malloc(3 * sizeof(double *));
   k_blk[0] = (double *) malloc((n+1)*(n+1)/4 * sizeof(double));
   k_blk[1] = (double *) malloc(n*(n+1)/2 * sizeof(double));
   k_blk[2] = (double *) malloc((n+1)*(n+1)/4 * sizeof(double));


   // Initialize 'a' and 'b' to 0.
   double ab[3][2] = {{0.0,0.0}, {0.0,0.0}, {0.0,0.0}};

/*
   * If 'fast' is true, a heuristic is used to speed up the
   * algorithm. The log-likelihood is computed by inserting a
   * single break and local maxima are used as only candidate
   * breakpoints. Because tadbit is O(n^2) the gain is of the
   * same order.
*/
   
   int bkpts[n];
   // By default, all breakpoints are candidates.
   for (i = 0 ; i < n ; i++) {
      bkpts[i] = 1;
   }

   // If 'fast', only local maxima are candidates.
   if (fast) {
      remove_non_local_maxima(obs, dis, n, m, k_blk, d_blk, bkpts);
   }


/*
   * Compute the log-likelihood of the segments. the element
   * (i,j) of the matrix-like array 'llik' will contain the
   * log-likelihood of the segment starting at i and ending
   * at j. the matrix is initialized with nan because not all
   * elements will be computed. the lower triangular part is
   * left out and possibily most of the elements if fast is
   * true.
*/

   for (i = 0 ; i < n-2 ; i++) {
      // Skip if not a potential breakpoint.
      if ((i > 0) && (!bkpts[i-1])) {
         continue;
      }

      for (j = i+2 ; j < n ; j++) {
         // Skip if not a potential breakpoint.
         if (!bkpts[j]) {
            continue;
         }
         
         // Segment the (i,j)-blocks.
         d_blk = break_in_blocks(dis, n, i, j, d_blk);

         llik[i+j*n] = 0.0;
         for (k = 0 ; k < m ; k++) {
            k_blk = break_in_blocks(obs[k], n, i, j, k_blk);
            // Get the likelihood per block and sum.
            llik[i+j*n] +=
                ml_ab(k_blk[0], d_blk[0], ab[0], i*(j-i+1))       / 2  +
                ml_ab(k_blk[1], d_blk[1], ab[1], (j-i)*(j-i+1)/2)      +
                ml_ab(k_blk[2], d_blk[2], ab[2], (n-j-1)*(j-i+1)) / 2;
         }
      }
   }


/*
   * The matrix 'llik' contains the log-likelihood of the
   * segments. The breakpoints are found by the dynamic
   * programming routine 'get_breakpoints'.
*/

   int *all_breakpoints = (int *) malloc(n * sizeof(n));
   all_breakpoints = get_breakpoints(llik, n, all_breakpoints);

   // Free allocated memory (not sure this is needed).
   for (i = 0 ; i < 3 ; i++) {
      free(d_blk[i]);
      free(k_blk[i]);
   }
   free(d_blk);
   free(k_blk);
   free(dis);
   free(llik);

   // Done!!
   return all_breakpoints;

}

// R/C interface.
SEXP tadbit_R_call(SEXP list, SEXP fast_yn) {
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
   int fast = INTEGER(fast_yn)[0];

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

   // Call 'tadbit'.
   int *bkpts = tadbit(obs, n, m, fast);

   // Wrap it up.
   SEXP return_val_sexp;
   PROTECT(return_val_sexp = allocVector(INTSXP, n));
   int *return_val = INTEGER(return_val_sexp);
   // Copy output from 'tadbit'.
   for (i = 0 ; i < n ; i++) {
      return_val[i] = bkpts[i];
   }
   free (bkpts);
   UNPROTECT(1);

   return return_val_sexp;

}

int *read_int_matrix (const char *fname, int *n_ptr) {
/*
   * Read integers from a file and return a line-matrix (a pointer
   * of pointers to int). Assume that the first line is a header
   * containing as many space characters as columns, that the number
   * of rows is the same as the number of columns, and that each
   * row starts with the row name.
*/

   int i, j = 0, check_row = 0, check_col, n = 0;
   char line[65536], row_name[256], *a, *z;

   int *array;
   long u;

   // Garbage collection in case of early return.
   int *failure_exit() {
      free(array);
      return NULL;
   }

   FILE *f = fopen(fname, "r");
   if (f == NULL) {
      return NULL;
   }

   // Use first line to get 'n' by counting tabs.
   fgets(line, sizeof(line), f);
   for (i = 0 ; i < strlen(line)-1 ; i++) {
      if (isspace(line[i])) {
         n++;
      }
   }

   array = (int *) malloc (n*n * sizeof(int));

   while (fgets(line, sizeof(line), f)) {
      if (++check_row > n) {
         printf("check_row = %d\n", check_row);
         return failure_exit();
      }

      a = line;
      z = NULL;

      // Skip row name.
      if (sscanf(line, "%s", row_name)) {
         a += strlen(row_name);
      }

      check_col = 0;
      u = strtol(a, &z, 10);
      while (a != z) {
         check_col++;
         array[j++] = (int) u;
         a = z;
         u = strtol(a, &z, 10);
      }

      if (check_col != n) {
         printf("n = %d, check_col = %d\n", n, check_col);
         return failure_exit();
      }

   }

   if (check_row != n) {
      return failure_exit();
   }

   fclose(f);
   *n_ptr = n;
   return array;

}


int main (int argc, const char* argv[]) {
   int i, n, *some_ints = read_int_matrix(argv[1], &n);
   if (some_ints != NULL) {
      double *counts = (double *) malloc(n*n * sizeof(double));
      for (i = 0 ; i < n*n ; i++) {
         counts[i] = (double) some_ints[i];
      }
      double **obs = &counts;
      int *breakpoints = tadbit(obs, n, 1, 1);
   }
}
