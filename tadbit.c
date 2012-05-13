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

void
get_breakpoints(
  double *llik,
  int n,
  int *all_breakpoints,
  int *minmax
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

   minmax[0] = minmax[1] = 0;

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

      if ((new_full_llik > -INFINITY) && (minmax[0] == 0)) {
         minmax[0] = nbreaks;
      }
      if (!(new_full_llik > -INFINITY) && (minmax[0] > 0)) {
         minmax[1] = nbreaks-1;
         break;
      }

   }

   if (minmax[1] == 0) {
      minmax[1] = n/4;
   }
   
   free(new_bkpt_list);
   free(old_bkpt_list);

}


void
free_ml_block(
  ml_block *blk
){
   int i;
   // Start releasing diagonal blocks.
   for (i = 0 ; i < blk->n_diag_blk; i++) {
      free(blk->k_diag[i]);
      free(blk->d_diag[i]);
   }
   free(blk->sizes_diag_blk);
   free(blk->k_diag);
   free(blk->d_diag);

   // Then release off-diagonal blocks.
   for (i = 0 ; i < blk->n_off_blk; i++) {
      free(blk->k_off[i]);
      free(blk->d_off[i]);
   }
   free(blk->sizes_off_blk);
   free(blk->k_off);
   free(blk->d_off);

   // Done.

}

void
malloc_and_dice(
  const double *obs,
  const double *dis,
  int bk_pts[],
  const int nbreaks,
  ml_block *blk,
  const int n
){

/*
   * Fill ml_block data according to breakpoints.
*/

   int a, b, c, i, j, k;
   int size;

   blk->n_diag_blk = nbreaks+1;
   blk->n_off_blk = nbreaks*(nbreaks+1)/2;

   blk->sizes_diag_blk = (int *) malloc(blk->n_diag_blk * sizeof(int));
   blk->sizes_off_blk = (int *) malloc(blk->n_off_blk * sizeof(int));
   blk->k_diag = (double **) malloc(blk->n_diag_blk * sizeof(double *));
   blk->d_diag = (double **) malloc(blk->n_diag_blk * sizeof(double *));
   blk->k_off = (double **) malloc(blk->n_off_blk * sizeof(double *));
   blk->d_off = (double **) malloc(blk->n_off_blk * sizeof(double *));

   // Start with diagonal blocks.
   for (c = 0 ; c < nbreaks+1 ; c++) {
      // Get the size of the block
      size = (bk_pts[c+1]-bk_pts[c]) * (bk_pts[c+1]-bk_pts[c]-1) / 2;

      // Allocate the memory.
      blk->k_diag[c] = (double *) malloc(size * sizeof(double));
      blk->d_diag[c] = (double *) malloc(size * sizeof(double));

      // Fill in the values.
      blk->sizes_diag_blk[c] = 0;
      for (k = 0, j = bk_pts[c] + 1 ; j < bk_pts[c+1] + 1 ; j++) {
      for (i = bk_pts[c] + 1 ; i < j ; i++) {
         if (isnan(obs[i+j*n])) {
            continue;
         }
         blk->k_diag[c][k] = obs[i+j*n];
         blk->d_diag[c][k] = dis[i+j*n];
         blk->sizes_diag_blk[c]++;
         k++;
      }
      }
   }

   // Then do off-diagonal blocks.
   for (c = 0, b = 0 ; b < nbreaks ; b++) {
   for (a = b+1 ; a < nbreaks+1 ; a++) {
      // Get the size of the block
      size = (bk_pts[a+1]-bk_pts[a]) * (bk_pts[b+1]-bk_pts[b]);

      // Allocate the memory.
      blk->k_off[c] = (double *) malloc(size * sizeof(double));
      blk->d_off[c] = (double *) malloc(size * sizeof(double));

      // Fill in the values.
      blk->sizes_off_blk[c] = 0;
      for (k = 0, j = bk_pts[a]+1 ; j < bk_pts[a+1] + 1 ; j++) {
      for (i = bk_pts[b]+1 ; i < bk_pts[b+1] + 1 ; i++) {
         if (isnan(obs[i+j*n])) {
            continue;
         }
         blk->k_off[c][k] = obs[i+j*n];
         blk->d_off[c][k] = dis[i+j*n];
         blk->sizes_off_blk[c]++;
         k++;
      }
      }
      c++;
   }
   }
   
   // Done.
   
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

   int *all_breakpoints = return_val;
   int minmax[2];
   get_breakpoints(llik, n, all_breakpoints, minmax);

   free(llik);
   llik = (double *) malloc(n/4 * sizeof(double));

   // Assign tasks to threads.
   assignment = (int *) malloc(n/4 * sizeof(int));
   processed = 0;
   to_process = minmax[1] - minmax[0] + 1;
   for (i = minmax[0] ; i < minmax[1] ; i++) {
      assignment[i] = processed % n_threads;
      processed++;
   }

   // 'all_breakpoints' consists of a matrix with breakpoint
   // position for each number of breakpoint.
   void *get_block_ml(void *arg) {

      int j, k, nbreaks;
      int bk_pts[n/4];
      double ab[2] = {0.0, 0.0};

      pthread_t myid = pthread_self();

      ml_block *blk = (ml_block *) malloc(sizeof(ml_block));

      // Readability variable.
      int assigned_to_me;
      for (nbreaks = minmax[0] ; nbreaks  < minmax[1] ; nbreaks++) {
         assigned_to_me = pthread_equal(myid, tid[assignment[nbreaks]]);
         if (!assigned_to_me) {
            continue;
         }

         // Recode breakpoint positions in an (n+2)-array.
         bk_pts[0] = -1;
         bk_pts[nbreaks+1] = n-1;
         for (k = 1, j = 0 ; j < n ; j++) {
            if (all_breakpoints[nbreaks+j*n] == 1) {
               bk_pts[k] = j;
               k++;
            }
         }

         llik[nbreaks] = 0.0;
         for (k = 0 ; k < m ; k++) {
            // Allocate memory and dice matrices.
            malloc_and_dice(obs[k], dis, bk_pts, nbreaks, blk, n);

            // Start with diagonal blocks.
	    for (j = 0 ; j < blk->n_diag_blk ; j++) {
               ab[0] = ab[1] = 0.0;
               llik[nbreaks] += poiss_reg(blk->k_diag[j],
                     blk->d_diag[j], ab, blk->sizes_diag_blk[j]);
	    }
            // Then off-diagonal blocks.
	    for (j = 0 ; j < blk->n_off_blk ; j++) {
               ab[0] = ab[1] = 0.0;
               llik[nbreaks] += poiss_reg(blk->k_off[j],
                     blk->d_off[j], ab, blk->sizes_off_blk[j]);
	    }

            // Free the memory (not efficent).
            free_ml_block(blk);
         }

         // Wrapper for complex memory relase.

      }

      return NULL;

   }

   // Start the threads once more...
   for (i = 0 ; i < n_threads ; i++) {
      if (pthread_create(&(tid[i]), NULL, &get_block_ml, NULL) != 0) {
         return;
      }
   }

   // ... and wait for them.
   for (i = 0 ; i < n_threads ; i++) {
      pthread_join(tid[i], NULL);
   }

   // START DEBUG.
   for (i = 0 ; i < minmax[1] ; i++) {
      printf("%f, ", llik[i]);
   }
   printf("\n");
   // END DEBUG.

   free(dis);

   // Done!!

}
