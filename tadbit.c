#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <pthread.h>
#include <ctype.h>
#include <unistd.h>

#include "tadbit.h"


double
poiss_reg (
  const double *k,
  const double *d,
  double *ab,
  const int n
){

/*
   * The 2-array 'ab' is upated in place and the log-likelihood
   * is returned.
   * The fitted model (by maximum likelihood) is Poisson with lambda
   * paramter such that lambda = exp(a + b*d). So the full log-likelihood
   * of the model is Sigma -exp(a + b*d_i) + k_i(a + b*d_i).
*/

   if (n < 3) {
      return 0.0;
   }

   int i;
   double a = ab[0], b = ab[1], llik;
   double da = 0.0, db = 0.0, oldgrad;
   double f, g, dfda = 0.0, dfdb = 0.0, dgda = 0.0, dgdb = 0.0;
   double denom;

   // 'tmp' is a computation intermediate that will be the return
   // value of 'exp'. This can call '__slowexp' which on 64-bit machines
   // can return a long double (causing segmentation fault if 'tmp' is
   // declared as long).
   long double tmp;

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

void
slice(
  const double *k,
  const double *d,
  const int n,
  const int i,
  const int j,
  ml_slice *blocks
){

/*
   *  Break up 'mat' in three blocks delimited by 'i' and 'j'.
   *  The upper block is (0,i-1)x(i,j), the triangular block is
   *  the upper triangular block without diagonal (i,j)x(i,j)
   *  and the bottom block is (j+1,n)x(i,j).
*/

   int l, row, col;

   for (l = 0 ; l < 3 ; l++) {
      blocks->size[l] = 0;
   }

   // Fill vertically.
   for (col = i ; col < j+1 ; col++) {
      // Skip if 'i' is 0.
      for (row = 0 ; row < i ; row++) {
         if (!isnan(k[row+col*n])) {
            blocks->k[0][blocks->size[0]] = k[row+col*n];
            blocks->d[0][blocks->size[0]] = d[row+col*n];
            blocks->size[0]++;
         }
      }

      // Skip if 'col' is i.
      for (row = i ; row < col ; row++) {
         if (!isnan(k[row+col*n])) {
            blocks->k[1][blocks->size[1]] = k[row+col*n];
            blocks->d[1][blocks->size[1]] = d[row+col*n];
            blocks->size[1]++;
         }
      }

      // Skip if 'j' is n-1.
      for (row = j+1 ; row < n ; row++) {
         if (!isnan(k[row+col*n])) {
            blocks->k[2][blocks->size[2]] = k[row+col*n];
            blocks->d[2][blocks->size[2]] = d[row+col*n];
            blocks->size[2]++;
         }
      }
   }
}

int
get_breakpoints(
  double *llik,
  int n,
  int m,
  int *all_breakpoints
){

/*
   * Dynamic programming algorithm. Compute the most likely position
   * of breakpoints for up to n/4 breakpoints.
   * Return the minimum number of breakpoint that has a defined
   * maximum log-likelihood.
*/

   int i,j, nbreaks = 0;
   int new_bkpt;

   double tmp;

   double new_llik[n];
   double old_llik[n];
   // Initialize to first line of 'llik' (the end may be NAN).
   for (j = 0 ; j < n ; j++) {
      old_llik[j] = llik[j*n];
      new_llik[j] = -INFINITY;
   }

   // Breakpoint lists. The first index is the end of the slice,
   // the second is 1 if this position is an end (breakpoint).
   // These must be allocated from the heap because 'n' can be large.
   int *new_bkpt_list = (int *) malloc(n*n * sizeof(int));
   int *old_bkpt_list = (int *) malloc(n*n * sizeof(int));

   // Initialize to 0.
   for (i = 0 ; i < n ; i++) {
   for (j = 0 ; j < n ; j++) {
      all_breakpoints[i+j*n] = new_bkpt_list[i+j*n] = \
            old_bkpt_list[i+j*n] = 0;
   }
   }

   double new_full_llik = old_llik[n-1];
   double old_full_llik = -INFINITY;

   int n_params;
   const int N = n*(n-1)/2;
   double BIC = -INFINITY;

   for (nbreaks = 1 ; nbreaks < n/4 ; nbreaks++) {
      // Update breakpoint lists.
      for (i = 0 ; i < n ; i++) {
      for (j = 0 ; j < n ; j++) {
         old_bkpt_list[i+j*n] = new_bkpt_list[i+j*n];
      }
      }

      // Cycle over end point 'j'.
      for (j = 3 * nbreaks + 2 ; j < n ; j++) {
         new_llik[j] = -INFINITY;
         new_bkpt = -1;

         // Cycle over start point 'i'.
         for (i = 3 * nbreaks ; i < j-3 ; i++) {

            // NAN if not a potential breakpoint (so evaluates to false).
            tmp = old_llik[i-1] + llik[i+j*n];
            if (tmp > new_llik[j]) {
               new_llik[j] = tmp;
               new_bkpt = i-1;
            }
         }

         // Update breakpoint list (skip if log-lik is undefined).
         if (new_llik[j] > -INFINITY) {
            for (i = 0 ; i < n ; i++) {
               new_bkpt_list[j+i*n] = old_bkpt_list[new_bkpt+i*n];
            }
            new_bkpt_list[j+new_bkpt*n] = 1;
         }

      }

      // Update full log-likelihoods.
      old_full_llik = new_full_llik;
      new_full_llik = new_llik[n-1];

      for (i = 0 ; i < n ; i++) {
         old_llik[i] = new_llik[i];
         all_breakpoints[nbreaks+i*n] = new_bkpt_list[n-1+i*n];
      }

      // Return when max BIC is reached.
      n_params = m * (8 + 6*(nbreaks-1));
      if (2*new_full_llik - n_params*log(N) < BIC) {
         free(new_bkpt_list);
         free(old_bkpt_list);
         return nbreaks-1;
      }
      else {
         BIC = 2*new_full_llik - n_params*log(N);
      }

   }

   // Did not maximize BIC.
   return -1;

}


void
tadbit(
  const double **obs,
  const int n,
  const int m,
  double max_tad_size,
  int n_threads,
  const int verbose,
  int *return_val
){

   int i, j, k;

   // Get thread number if set to 0(automatic).
   if (n_threads < 1) {
      #ifdef _SC_NPROCESSORS_ONLN
         n_threads = (int) sysconf(_SC_NPROCESSORS_ONLN);
      #else
         n_threads = 1;
      #endif
   }

   // Get absolute max tad size.
   max_tad_size = max_tad_size > 0 ? 
      max_tad_size < 1 ? max_tad_size *n : max_tad_size :
      INFINITY;

   // Create an array of thread ids.
   pthread_t tid[n_threads];
   int *assignment = (int *) malloc(n*n * sizeof(int));
   int processed = 0, to_process = 0;
   
   // The pointer 'assignment' specified which thread has to
   // compute which part of 'llik'.
   // NB: this is not optimal as this may cause some threads
   // to lag behind if their task takes more time. However this
   // is much simpler to implement and prevents clobbering and
   // memory errors.
   for (i = 0 ; i < n-3 ; i++) {
   for (j = i+3 ; j < n ; j++) {
      if (j-i < max_tad_size) {
         assignment[i+j*n] = to_process % n_threads;
         to_process++;
      }
      else {
         // Do not process if slice too large.
         assignment[i+j*n] = -1;
      }
   }
   }

   // Allocate memory and initialize variables. The distance
   // matrix 'dis' is the distance to the main diagonal. Every
   // element of coordinate (i,j) is on a diagonal; the distance
   // is the log-shift to the main diagonal |i-j|.
   // 'ab' contains parameters 'a' and 'b' for the maximum likelihood
   // model. Because each segmentation defines 3 regions we need
   // 3 pairs of parameters.


   double *dis = (double *) malloc(n*n * sizeof(double));
   double *llik = (double *) malloc(n*n * sizeof(double));

   for (k = i = 0; i < n ; i++) {
   for (j = 0; j < n ; j++) {
      llik[k] = NAN;
      dis[k] = log(1+abs(i-j));
      k++;
   }
   }

   // Inline function for threads.
   void *fill_matrix(void *arg) {
 
   /*
      * Compute the log-likelihood of the slices. the element
      * (i,j) of the matrix-like pointer 'llik' will contain the
      * log-likelihood of the slice starting at i and ending
      * at j. the matrix is initialized with nan because not all
      * elements will be computed. the lower triangular part is
      * left out.
   */
  
      int i, j, k;
      pthread_t myid = pthread_self();
   
      // Allocate max possible size to blocks.
      ml_slice *slc = (ml_slice *) malloc(sizeof(ml_slice));
      // Readability variable.
      int nmax = (n+1)*(n+1)/4;
      
      slc->k[0] = (double *) malloc(nmax     * sizeof(double));
      slc->k[1] = (double *) malloc(nmax * 2 * sizeof(double));
      slc->k[2] = (double *) malloc(nmax     * sizeof(double));
      slc->d[0] = (double *) malloc(nmax     * sizeof(double));
      slc->d[1] = (double *) malloc(nmax * 2 * sizeof(double));
      slc->d[2] = (double *) malloc(nmax     * sizeof(double));
       
      // Initialize 'a' and 'b' to 0.
      double ab[3][2] = {{0.0,0.0}, {0.0,0.0}, {0.0,0.0}};
      
      // Readability variables.
      int do_not_process, assigned_to_me;
      
      for (i = 0 ; i < n-3 ; i++) {
      for (j = i+3 ; j < n ; j++) {
         do_not_process = assignment[i+j*n] < 0;
         assigned_to_me = pthread_equal(myid, tid[assignment[i+j*n]]);
         if (do_not_process || !assigned_to_me) {
            continue;
         }
         // Distinct parts of the array, not lock needed.
         llik[i+j*n] = 0.0;
         for (k = 0 ; k < m ; k++) {
            slice(obs[k], dis, n, i, j, slc);
            // Get the likelihood per slice and sum.
            llik[i+j*n] +=
               poiss_reg(slc->k[0], slc->d[0], ab[0], slc->size[0]) / 2 +
               poiss_reg(slc->k[1], slc->d[1], ab[1], slc->size[1])     +
               poiss_reg(slc->k[2], slc->d[2], ab[2], slc->size[2]) / 2;
         }
         // No synchronization needed (no mutex).
         processed++;
         if (verbose) {
            fprintf(stderr, "Computing likelihood (%0.f%% done)\r",
               99 * processed / (float) to_process);
         }
      }
      } // End of the (i,j) for loop.

      // Free allocated memory.
      for (i = 0 ; i < 3 ; i++) {
         free(slc->k[i]);
         free(slc->d[i]);
      }
      free(slc);
   
      // Thread exit.
      return NULL;
   
   }

  
   // Start the threads.
   for (i = 0 ; i < n_threads ; i++) {
      if (pthread_create(&(tid[i]), NULL, &fill_matrix, NULL) != 0) {
         return;
      }
   }

   // Wait for threads to return.
   for (i = 0 ; i < n_threads ; i++) {
      pthread_join(tid[i], NULL);
   }
   if (verbose) {
      fprintf(stderr, "Computing likelihood (100%% done)\n");
   }

   free(assignment);

   // The matrix 'llik' now contains the log-likelihood of the
   // segments. The breakpoints are found by dynamic
   // programming.

   int *all_breakpoints = (int *) malloc(n*n * sizeof(int));
   int nbreaks = get_breakpoints(llik, n, m, all_breakpoints);

   for (i = 0 ; i < n ; i++) {
      return_val[i] = all_breakpoints[nbreaks+i*n];
   }

   free(all_breakpoints);
   free(dis);
   free(llik);

   // Done!!

}
