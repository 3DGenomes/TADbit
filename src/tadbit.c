#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <pthread.h>
#include <ctype.h>
#include <unistd.h>

#include "tadbit.h"

#define FIT_SLICE(slc,ab) \
   poiss_reg(slc->k[0], slc->d[0], slc->w[0], ab[0], slc->size[0]) / 2 + \
   poiss_reg(slc->k[1], slc->d[1], slc->w[1], ab[1], slc->size[1])     + \
   poiss_reg(slc->k[2], slc->d[2], slc->w[2], ab[2], slc->size[2]) / 2

double
poiss_reg (
  const double *k,
  const double *d,
  const double *w,
  double *ab,
  const int n
){

/*
   * Arguments:
   *   k:  Counts (size 'n').
   *   d:  Distances (size 'n').
   *   w:  Weights (size 'n').
   *   ab: Starting values of 'a' and 'b' (size 2).
   *   n:  Length of the arrays 'k', 'd' and 'w'.
   *
   * Return:
   *   Log-likelihood.
   *
   * Side-effects:
   *   ab:  updated in place with the estimates.
   *
   * The 2-array 'ab' is upated in place and the log-likelihood
   * is returned.
   * The fitted model (by maximum likelihood) is Poisson with lambda
   * paramter such that lambda = w*exp(a + b*d). So the full
   * log-likelihood of the model is the sum of terms
   *
   *         - w_i exp(a + b*d_i) + k_i(log(w_i) + a + b*d_i)
   *
*/

   if (n < 1) {
      return 0.0;
   }
   if (n < 3) {
      return NAN;
   }

   int i, iter = 0;
   double a = ab[0], b = ab[1], llik;
   double da = 0.0, db = 0.0, oldgrad;
   double f, g, dfda = 0.0, dfdb = 0.0, dgda = 0.0, dgdb = 0.0;
   double denom;

   // 'tmp' is a computation intermediate that will be the return
   // value of 'exp'. This can call '__slowexp' which on 64-bit machines
   // can return a long double (causing segmentation fault if 'tmp' is
   // declared as long).
   long double tmp;

   // Comodity function. The function f is -dl/da and g is -dl/db.
   void recompute_fg(void) {
      f = 0.0; g = 0.0;
      for (i = 0 ; i < n ; i++) {
         tmp  =  w[i] * exp(a+da+(b+db)*d[i]) - k[i];
         f   +=  tmp;
         g   +=  tmp * d[i];
      }
   }

   recompute_fg();

   // Newton-Raphson until gradient function < TOLERANCE.
   while ((oldgrad = f*f + g*g) > TOLERANCE && iter++ < MAXITER) {

      // Compute the derivatives.
      dfda = dfdb = dgda = dgdb = 0.0;
      for (i = 0 ; i < n ; i++) {
         tmp   =   w[i] * exp(a+b*d[i]);
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

   if (iter < MAXITER) {
      // Compute log-likelihood (using 'dfda').
      llik = 0.0;
      for (i = 0 ; i < n ; i++) {
         llik += exp(a+b*d[i]) + k[i] * (a + b*d[i]);
      }
      // Update 'ab' in place (to make the estimates available).
      ab[0] = a; ab[1] = b;
   }
   else {
      // Something probably went wrong. Reset 'ab' and return NAN.
      llik = NAN;
      ab[0] = ab[1] = 0.0;
   }

   return llik;

}

void
slice(
  const double *obs,
  const double *dis,
  const int n,
  const int start,
  const int end,
  ml_slice *blocks
){

/*
   * Arguments:
   *   obs:    Observation array (size 'n'^2).
   *   dis:    Distance array (size 'n'^2).
   *   n:      Size of th 'obs' and 'dis'.
   *   start:  Start index of the slice (included).
   *   end:    Stop index of the slice (included).
   *   blocks: Where to write the 3 blocks of the slice.
   * Side-effects:
   *   Update 'blocks' in place.
   *    
   * Observations are weighted by the geometric mean of the counts on
   * the diagonals.
*/

   int l, row, col;

   // Comodity function for dryness.
   void add_to_block(int l) {
      if (!isnan(obs[row+col*n])) {
         blocks->k[l][blocks->size[l]] = obs[row+col*n];
         blocks->d[l][blocks->size[l]] = dis[row+col*n];
         blocks->w[l][blocks->size[l]] = \
             sqrt(obs[row+row*n]*obs[col+col*n]);
         blocks->size[l]++;
      }
   }

   for (l = 0 ; l < 3 ; l++) {
      blocks->size[l] = 0;
   }

   // Fill vertically.
   for (col = start ; col < end+1 ; col++) {
      // Skipped if 'start' is 0.
      for (row = 0 ; row < start ; row++) {
         add_to_block(0);
      }
      // Skipped if 'col' is 'start'.
      for (row = start ; row < col ; row++) {
         add_to_block(1);
      }
      // Skipped if 'end' is 'n-1'.
      for (row = end+1 ; row < n ; row++) {
         add_to_block(2);
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

   // Breakpoint lists. The first index (line) is the end of the slice,
   // the second is 1 if this position is an end (breakpoint).
   // These must be allocated from the heap because 'n' can be large.
   int *new_bkpt_list = (int *) malloc(n*n * sizeof(int));
   int *old_bkpt_list = (int *) malloc(n*n * sizeof(int));

   // Initialize to 0.
   for (i = 0 ; i < n*n ; i++) {
      all_breakpoints[i] = new_bkpt_list[i] = old_bkpt_list[i] = 0;
   }

   double new_full_llik = old_llik[n-1];

   int n_params;
   const int N = n*(n-1)/2;
   double BIC = -INFINITY;

   for (nbreaks = 1 ; nbreaks < n/4 ; nbreaks++) {
      // Update breakpoint lists.
      for (i = 0 ; i < n*n ; i++) {
         old_bkpt_list[i] = new_bkpt_list[i];
      }

      // Cycle over end point 'j'.
      for (j = 3 * nbreaks + 2 ; j < n ; j++) {
         new_llik[j] = -INFINITY;
         new_bkpt = -1;

         // Cycle over start point 'i'.
         for (i = 3 * nbreaks ; i < j-3 ; i++) {

            // NAN if not a potential breakpoint, so the following
            // lines evaluates to false.
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
      new_full_llik = new_llik[n-1];

      for (i = 0 ; i < n ; i++) {
         old_llik[i] = new_llik[i];
         all_breakpoints[nbreaks+i*n] = new_bkpt_list[n-1+i*n];
      }

      // Return when max BIC is reached.
      n_params = nbreaks + m * (8 + 6*(nbreaks-1));
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
  double **obs,
  int n,
  const int m,
  double max_tad_size,
  int n_threads,
  const int verbose,
  int *return_val
){


   // Get absolute max tad size.
   max_tad_size = max_tad_size > 0 ? 
      max_tad_size < 1 ? max_tad_size *n : max_tad_size :
      INFINITY;

   int i, j, k, l;
   const int init_n = n;

   // Allocate memory and initialize variables. The distance
   // matrix 'dis' is the distance to the main diagonal. Every
   // element of coordinate (i,j) is on a diagonal; the distance
   // is the log-shift to the main diagonal 'i-j'.
   // 'ab' contains parameters 'a' and 'b' for the maximum likelihood
   // model. Because each segmentation defines 3 regions we need
   // 3 pairs of parameters.

   double *init_dis = (double *) malloc(init_n*init_n * sizeof(double));

   for (l = 0, i = 0; i < init_n ; i++) {
   for (j = 0; j < init_n ; j++) {
      init_dis[l] = log(1+abs(i-j));
      l++;
   }
   }

   // Simplify input. Remove line and column if 0 on the diagonal.
   int remove[init_n];
   for (i = 0 ; i < init_n ; i++) {
      remove[i] = 0;
      for (k = 0 ; k < m ; k++) {
         if (obs[k][i+i*init_n] < 1.0) {
            remove[i] = 1;
         }
      }
   }

   // Update the dimension.
   for (i = 0 ; i < init_n ; i++) {
      n -= remove[i];
   }

   // Allocate and copy.
   double **new_obs = (double **) malloc(m * sizeof(double *));
   double *dis = (double *) malloc(n*n * sizeof(double));
   for (k = 0 ; k < m ; k++) {
      l = 0;
      new_obs[k] = (double *) malloc(n*n * sizeof(double));
      for (j = 0 ; j < init_n ; j++) {
      for (i = 0 ; i < init_n ; i++) {
         if (remove[i] || remove[j]) {
            continue;
         }
         new_obs[k][l] = obs[k][i+j*init_n];
         dis[l] = init_dis[i+j*init_n];
         l++;
      }
      }
   }

   // We will not need the initial observations any more.
   free(init_dis);
   obs = new_obs;


   // Get thread number if set to 0 (automatic).
   if (n_threads < 1) {
      #ifdef _SC_NPROCESSORS_ONLN
         n_threads = (int) sysconf(_SC_NPROCESSORS_ONLN);
      #else
         n_threads = 1;
      #endif
   }


   double *llik = (double *) malloc(n*n * sizeof(double));
   for (l = 0 ; l < n*n ; l++) {
      llik[l] = NAN;
   }

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
  
      int i = 0, j, k;
      pthread_t myid = pthread_self();
   
      // Allocate max possible size to blocks.
      ml_slice *slc = (ml_slice *) malloc(sizeof(ml_slice));
      // Readability variable.
      int nmax = (n+1)*(n+1)/4;

      // TODO: Reduce max size to minimum actually used.
      // My last attempts caused segmentation faults.
      slc->k[0] = (double *) malloc(nmax     * sizeof(double));
      slc->k[1] = (double *) malloc(nmax * 2 * sizeof(double));
      slc->k[2] = (double *) malloc(nmax     * sizeof(double));
      slc->d[0] = (double *) malloc(nmax     * sizeof(double));
      slc->d[1] = (double *) malloc(nmax * 2 * sizeof(double));
      slc->d[2] = (double *) malloc(nmax     * sizeof(double));
      slc->w[0] = (double *) malloc(nmax     * sizeof(double));
      slc->w[1] = (double *) malloc(nmax * 2 * sizeof(double));
      slc->w[2] = (double *) malloc(nmax     * sizeof(double));
       
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

         // Distinct parts of the array, no lock needed.
         llik[i+j*n] = 0.0;
         for (k = 0 ; k < m ; k++) {
            // Get the (i,j) slice (stored in 'slc').
            slice(obs[k], dis, n, i, j, slc);
            // Get the likelihood and sum (see macro definition).
            llik[i+j*n] += FIT_SLICE(slc, ab);
         }
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

   for (l = 0, i = 0 ; i < init_n ; i++) {
      if (remove[i]) {
         return_val[i] = 0;
      }
      else {
         return_val[i] = all_breakpoints[nbreaks+l*n];
         l++;
      }
   }

   for (k = 0 ; k < m ; k++) {
      free(new_obs[k]);
   }
   free(new_obs);
   free(all_breakpoints);
   free(dis);
   free(llik);

   // Done!! The results is in 'return_val'.

}
