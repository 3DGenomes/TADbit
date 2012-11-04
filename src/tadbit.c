#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include <ctype.h>
#include <unistd.h>

#include "tadbit.h"


// Globals variables. //

double **globs;               // Global pointer to observations.
double **globw;               // Global pointer to weights.
int globm;                    // Number of replicates.
int n_processed;              // Number of slices processed so far.
int n_to_process;             // Total number of slices to process.
int taskQ_i;                  // Index used for task queue.
pthread_mutex_t tadbit_lock;  // Mutex to access task queue.


/*
int
obscmp (
  const void *a,
  const void *b
){
   int k;
   double x = 0.0;
   double y = 0.0; 

   for (k = 0 ; k < globm ; k++) {
      x += globs[k][*(int *)a] / globw[k][*(int *)a];
      y += globs[k][*(int *)b] / globw[k][*(int *)b];
   }

   // Sort weighted observations in descending order.
   if (x == y) return 0;
   return (y > x) ? 1 : -1;
}
*/


void
recompute_fg(
  /* input */
  const ml_block *blk,
  const double a,
  const double b,
  const double da,
  const double db,
  /* output */
  double *f,
  double *g
){
// SYNOPSIS:                                                            
//   Subfroutine of 'poiss_reg' that computes 'f' and 'g' in Newton-    
//   Raphson cycles.                                                    
//                                                                      
// ARGUMENTS:                                                           
//   '*blk': the block to fit by Poisson regression.                    
//   'a': parameter 'a' of the Poisson regression (see 'poiss_reg').    
//   'b': parameter 'b' of the Poisson regression (see 'poiss_reg').    
//   'da': computed differential of 'a' (see 'poiss_reg').              
//   'db': computed differential of 'b' (see 'poiss_reg').              
//        -- output arguments --                                        
//   'f': first function to zero, recomputed by the routine.            
//   'g': second function to zero, recomputed by the routine.           
//                                                                      
// RETURN:                                                              
//   'void'                                                             
//                                                                      
// SIDE-EFFECTS:                                                        
//   Update 'f' and 'g' in place.                                       
//                                                                      

   const int n = blk->size;
   const double *k = blk->counts;
   const double *d = blk->dist;
   const double *w = blk->weights;

   // 'tmp' is a computation intermediate that will be the return
   // value of 'exp'. This can call '__slowexp' which on 64-bit machines
   // can return a long double (causing segmentation fault if 'tmp' is
   // declared as long).
   long double tmp;
   int i;

   *f = 0.0; *g = 0.0;
   for (i = 0 ; i < n ; i++) {
      tmp  =  w[i] * exp(a+da+(b+db)*d[i]) - k[i];
      *f  +=  tmp;
      *g  +=  tmp * d[i];
   }

   return;

}


double
poiss_reg (
  const ml_block *blk
){
// SYNOPSIS:                                                            
//   The fitted model (by maximum likelihood) is Poisson with lambda    
//   paramter such that lambda = w * exp(a + b*d). So the full          
//   log-likelihood of the model is the sum of terms                    
//                                                                      
//         - w_i exp(a + b*d_i) + k_i(log(w_i) + a + b*d_i)             
//                                                                      
// ARGUMENTS:                                                           
//   '*blk': the block to fit by Poisson regression.                    
//                                                                      
// RETURN:                                                              
//   The maximum log-likelihood of a block of hiC data.                 
//                                                                      

   const int n = blk->size;
   const double *log_gamma = blk->lgamma;
   const double *k = blk->counts;
   const double *d = blk->dist;
   const double *w = blk->weights;

   if (n < 1) {
      return 0.0;
   }
   if (n < 3) {
      return NAN;
   }

   int i;
   int iter = 0;
   double denom;
   double oldgrad;
   double f;
   double g;
   double a = 0.0;
   double b = 0.0;
   double da = 0.0;
   double db = 0.0;
   double dfda = 0.0;
   double dfdb = 0.0;
   double dgda = 0.0;
   double dgdb = 0.0;
   // See the comment about 'tmp' in 'recompute_fg'.
   long double tmp; 

   recompute_fg(blk, a, b, da, db, &f, &g);

   // Newton-Raphson until gradient function is less than TOLERANCE.
   // The gradient function is the square norm 'f*f + g*g'.
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

      recompute_fg(blk, a, b, da, db, &f, &g);

      // Traceback if we are not going down the gradient. Cut the
      // length of the steps in half until this step goes down
      // the gradient.
      while (f*f + g*g > oldgrad) {
         da /= 2;
         db /= 2;
         recompute_fg(blk, a, b, da, db, &f, &g);
      }

      // Update 'a' and 'b'.
      a += da;
      b += db;

   }

   if (iter >= MAXITER) {
      // Something probably went wrong. Return NAN.
      return NAN;
   }

   // Compute log-likelihood (using 'dfda').
   double llik = 0.0;
   for (i = 0 ; i < n ; i++) {
      llik += exp(a+b*d[i]) + k[i] * (a + b*d[i]) - log_gamma[i];
   }

   return llik;

}

double
fit_slice(
  ml_block *blocks[3]
){
// SYNOPSIS:                                                            
//   Wrapper for 'poiss_reg'. Fits the three regions of a slice         
//   by Poisson regression and return the log-likelihood.               
//                                                                      
// PARAMETERS:                                                          
//   '*blocks' : the slice to be fitted by Poisson regression.          
//                                                                      
// RETURN:                                                              
//   The total maximum log-likelihood of the slice.                     
//                                                                      

   double top = poiss_reg(blocks[0]);
   double mid = poiss_reg(blocks[1]);
   double bot = poiss_reg(blocks[2]);

   // The likelihood of 'top' and 'bot' blocks are divided by 2 because
   // the hiC map is assumed to be symmetric. Those data points are
   // used two times in a segmentation, while the data points of the
   // 'mid' block are used only one time.

   return top/2 + mid + bot/2;

}


void
slice(
  /* input */
  const double *log_gamma,
  const double *obs,
  const double *dist,
  const int n,
  const int start,
  const int end,
  /* output */
  ml_block *blocks[3]
){
// SYNOPSIS:                                                            
//   Extract from a single replicate of the hiC data the three blocks   
//   of a slice delimited by start and end positions.                   
//                                                                      
// PARAMETERS:                                                          
//   '*log_gamma' : pre-computed log-gammas of the counts.              
//   '*obs' : raw hiC count data (single replicate).                    
//   '*dist' : linear separations corresponding to counts.              
//   'n' : row/col number of the full data matrix.                      
//   'start' : start index of the slice.                                
//   'end' : end index of the slice.                                    
//        -- output arguments --                                        
//   '*blocks' : variable to store the slice.                           
//                                                                      
// RETURN:                                                              
//   'void'                                                             
//                                                                      
// SIDE-EFFECTS:                                                        
//   Update 'blocks' in place.                                          
//                                                                      

   int l;
   int pos;
   int row;
   int col;

   // Initialize block sizes to 0.
   for (l = 0 ; l < 3 ; l++) {
     blocks[l]->size = 0;
   }

   // Iterate over data vertically.
   for (col = start ; col < end+1 ; col++) {
   for (row = 0 ; row < n ; row++) {

      // Skip NAs in the data.
      if (isnan(obs[row+col*n]))
         continue;

      // Find which block to update ('l').
      // 0: top, 1: middle, 2: bottom, or none.                    
      if        (row < start)  l = 0;
      else if   (row < col)    l = 1;
      else if   (row > end)    l = 2;
      else                     continue;

      pos = blocks[l]->size;

      blocks[l]->counts[pos] = log_gamma[row+col*n];
      blocks[l]->counts[pos] = obs[row+col*n];
      blocks[l]->dist[pos] = dist[row+col*n];
      // The weight is the square root of the product of the
      // diagonal terms (the product of potencies).
      blocks[l]->weights[pos] = sqrt(obs[row+row*n]*obs[col+col*n]);

      blocks[l]->size++;

   }
   }
   // End of the double for loop.

   return;

}

void
mlwalk(
  /* input */
  const double *llik_mat,
  const int n,
  const int m,
  const int maxbreaks,
  /* output */
  double *mllik,
  int *breakpoints
){
// SYNOPSIS:                                                            
//   Dynamic programming algorithm to compute the most likely position  
//   of breakpoints given a matrix of slice maximum log-likelihood.     
//                                                                      
// PARAMETERS:                                                          
//   '*llik_mat' : matrix of maximum log-likelihood values.             
//   'n' : row/col number of 'llik_mat'.                                
//   'm' : number of replicates.                                        
//        -- output arguments --                                        
//   '*mllik' : maximum log-likelihood of the segmentations.            
//   '*breakpoints' : optimal breakpoints per number of breaks.         
//                                                                      
// RETURN:                                                              
//   The optimal number of breakpoints.                                 
//                                                                      
// SIDE-EFFECTS:                                                        
//   Update 'breakpoints' in place.                                     
//                                                                      

   int i;
   int j;
   int nbreaks;

   double new_llik[n];
   double old_llik[n];

   // Breakpoint lists. The first index (row) is the end of the slice,
   // the second is 1 if this position is an end (breakpoint).
   // These must be allocated from the heap because 'n' can be large.
   int *new_bkpt_list = (int *) malloc(n*n * sizeof(int));
   int *old_bkpt_list = (int *) malloc(n*n * sizeof(int));

   // Initializations.
   // 'breakpoints' is a 'n' x 'maxbreaks' array. The first index (row)
   // is 1 if there is a breakpoint at that location, the second index
   // (column) is the number of breakpoints.
   memset(breakpoints, 0, n*maxbreaks * sizeof(int));

   for (i = 0 ; i < n*n ; i++) {
      new_bkpt_list[i] = 0;
      old_bkpt_list[i] = 0;
   }

   for (i = 0 ; i < maxbreaks ; i++) {
      mllik[i] = NAN;
   }

   // Initialize 'old_llik' to the first line of 'llik_mat' containing
   // the log-likelihood of segments starting at index 0.
   for (i = 0 ; i < n ; i++) {
      old_llik[i] = llik_mat[i*n];
      new_llik[i] = -INFINITY;
   }


   // Dynamic programming.
   for (nbreaks = 1 ; nbreaks < maxbreaks ; nbreaks++) {
      // Update breakpoint lists.
      for (i = 0 ; i < n*n ; i++) {
         old_bkpt_list[i] = new_bkpt_list[i];
      }

      // Cycle over end point 'j'.
      for (j = 3 * nbreaks + 2 ; j < n ; j++) {
         new_llik[j] = -INFINITY;
         int new_bkpt = -1;

         // Cycle over start point 'i'.
         for (i = 3 * nbreaks ; i < j-3 ; i++) {

            // If NAN the following condition evaluates to false.
            double tmp = old_llik[i-1] + llik_mat[i+j*n];
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
      mllik[nbreaks] = new_llik[n-1];

      // Record breakpoints.
      for (i = 0 ; i < n ; i++) {
         old_llik[i] = new_llik[i];
         breakpoints[i+nbreaks*n] = new_bkpt_list[n-1+i*n];
      }

   }

   free(new_bkpt_list);
   free(old_bkpt_list);

   return;

}

void *
fill_llikmat(
   void *arg
){
// SYNOPSIS:                                                            
//   Compute the log-likelihood of the slices. The element (i,j) of     
//   the matrix 'llikmat' will contain the log-likelihood  of the       
//   slice starting at i and ending at j. the matrix is initialized     
//   with nan because not all elements will be computed. The lower      
//   triangular part is left out.                                       
//                                                                      
// PARAMETERS:                                                          
//   '*arg' : matrix of maximum log-likelihood values.                  
//                                                                      
// RETURN:                                                              
//   void                                                               
//                                                                      
// SIDE-EFFECTS:                                                        
//   Update 'llikmat' in place.                                         
//                                                                      

   thread_arg *myargs = (thread_arg *) arg;
   const int n = myargs->n;
   const int m = myargs->m;
   const double **obs = (const double **) myargs->obs;
   const double *dist = (const double*) myargs->dist;
   const double **log_gamma = (const double **) myargs->log_gamma;
   const int *skip = (const int *) myargs->skip;
   double *llikmat = myargs->llikmat;
   const int verbose = myargs->verbose;

   int i;
   int j;
   int k;

   int job_index;

   // Allocate max possible size to blocks.
   ml_block *blocks[3];

   // The middle block can be bigger than top and bottom.
   int nmax[3] = { (n+1)*(n+1)/4, (n+1)*(n+1)/2, (n+1)*(n+1)/4 };

   // Allocate and initialize 3 blocks of a slice.
   for (i = 0 ; i < 3 ; i++) {
      blocks[i] = (ml_block *) malloc(sizeof(ml_block));

      blocks[i]->lgamma  = (double *) malloc(nmax[i] * sizeof(double));
      blocks[i]->counts  = (double *) malloc(nmax[i] * sizeof(double));
      blocks[i]->dist    = (double *) malloc(nmax[i] * sizeof(double));
      blocks[i]->weights = (double *) malloc(nmax[i] * sizeof(double));

      memset(blocks[i]->lgamma,  0.0, nmax[i] * sizeof(double));
      memset(blocks[i]->counts,  0.0, nmax[i] * sizeof(double));
      memset(blocks[i]->dist,    0.0, nmax[i] * sizeof(double));
      memset(blocks[i]->weights, 0.0, nmax[i] * sizeof(double));
   }
   
   // Break out of the loop when task queue is empty.
   while (1) {

      pthread_mutex_lock(&tadbit_lock);
      while ((taskQ_i < n*n) && (skip[taskQ_i] > 0)) {
         taskQ_i++;
      }
      if (taskQ_i >= n*n) {
         // Task queue is empty. Exit while loop and wrap up.
         pthread_mutex_unlock(&tadbit_lock);
         break;
      }
      job_index = taskQ_i;
      taskQ_i++;
      pthread_mutex_unlock(&tadbit_lock);

      i = job_index % n;
      j = job_index / n;

      // Distinct parts of the array, no lock needed.
      llikmat[i+j*n] = 0.0;
      for (k = 0 ; k < m ; k++) {
         // Get the (i,j) slice (stored in 'blocks').
         slice(log_gamma[k], obs[k], dist, n, i, j, blocks);
         // Get the likelihood and sum (see macro definition).
         llikmat[i+j*n] += fit_slice(blocks);
      }
      n_processed++;
      if (verbose) {
         fprintf(stderr, "computing likelihood (%0.f%% done)\r",
            99 * n_processed / (float) n_to_process);
      }
   }

   // Free allocated memory.
   for (i = 0 ; i < 3 ; i++) {
      free(blocks[i]->lgamma);
      free(blocks[i]->counts);
      free(blocks[i]->dist);
      free(blocks[i]->weights);
      free(blocks[i]);
   }

   return NULL;

}



void
tadbit(
  /* input */
  double **obs,
  int n,
  const int m,
  int n_threads,
  const int verbose,
  //const int max_tad_size,
  int max_tad_size,
  const int heuristic,
  /* output */
  tadbit_output *seg
){

   const int N = n; // Original size.
   int errno;       // Used for error checking.

   int i;
   int j;
   int k;
   int l;


   // Allocate memory and initialize variables. The distance
   // matrix 'dist' is the distance to the main diagonal. Every
   // element of coordinate (i,j) is on a diagonal; the distance
   // is the log-shift to the main diagonal 'i-j'.

   double *init_dist = (double *) malloc(N*N * sizeof(double));

   for (l = 0, i = 0; i < N ; i++) {
   for (j = 0; j < N ; j++) {
      init_dist[l] = log(abs(i-j));
      l++;
   }
   }

   // Simplify input. Remove line and column if 0 on the diagonal.
   int remove[N];
   for (i = 0 ; i < N ; i++) {
      remove[i] = 0;
      for (k = 0 ; k < m ; k++) {
         if (obs[k][i+i*N] < 1.0) {
            remove[i] = 1;
         }
      }
   }

   // Update the dimension. 'N' is the original row/column number,
   // 'n' is the row/column number after removing rows and columns
   // with 0 on the diagonal.
   for (i = 0 ; i < N ; i++) {
      n -= remove[i];
   }

   // Allocate and copy.
   double **log_gamma = (double **) malloc(m * sizeof(double *));
   double **new_obs = (double **) malloc(m * sizeof(double *));
   double *dist = (double *) malloc(n*n * sizeof(double));
   for (k = 0 ; k < m ; k++) {
      l = 0;
      log_gamma[k] = (double *) malloc(n*n * sizeof(double));
      new_obs[k] = (double *) malloc(n*n * sizeof(double));
      for (j = 0 ; j < N ; j++) {
      for (i = 0 ; i < N ; i++) {
         if (remove[i] || remove[j]) {
            continue;
         }
         log_gamma[k][l] = lgamma(obs[k][i+j*N]);
         new_obs[k][l] = obs[k][i+j*N];
         dist[l] = init_dist[i+j*N];
         l++;
      }
      }
   }

   // We will not need the initial observations any more.
   free(init_dist);
   obs = new_obs;


   // Compute row/column sums (identical by symmetry).
   double **rowsums = (double **) malloc(m * sizeof(double *));
   for (k = 0 ; k < m ; k++) {
      rowsums[k] = (double *) malloc(n * sizeof(double));
      memset(rowsums[k], 0.0, n * sizeof(double));
   }
   for (i = 0 ; i < n ; i++) {
   for (k = 0 ; k < m ; k++) {
   for (l = 0 ; l < n ; l++) {
      rowsums[k][i] += obs[k][i+l*n];
   }
   }
   }

   // Compute the weights.
   double **weights = (double **) malloc(m * sizeof(double *));
   for (k = 0 ; k < m ; k++) {
      weights[k] = (double *) malloc(n*n * sizeof(double));
      memset(weights[k], 0.0, n*n * sizeof(double));
      // Compute scalar product.
      for (i = 0 ; i < n ; i++) {
      for (j = 0 ; j < n ; j++) {
         weights[k][i+j*n] = sqrt(rowsums[k][i]*rowsums[k][j]);
      }
      }
   }

   // We don't need the row/column sums any more.
   for (k = 0 ; k < m ; k++) {
      free(rowsums[k]);
   }
   free(rowsums);

   // Get thread number if set to 0 (automatic).
   if (n_threads < 1) {
      #ifdef _SC_NPROCESSORS_ONLN
         n_threads = (int) sysconf(_SC_NPROCESSORS_ONLN);
      #else
         n_threads = 1;
      #endif
   }

   double *llikmat = (double *) malloc(n*n * sizeof(double));
   for (l = 0 ; l < n*n ; l++) {
      llikmat[l] = NAN;
   }

   int *skip = (int *) malloc(n*n *sizeof(int));
   memset(skip, 1, n*n*sizeof(int));
   n_processed = 0;
   n_to_process = 0;
   for (i = 0 ; i < n-3 ; i++) {
   for (j = i+3 ; j < n ; j++) {
      skip[i+j*n] = 0;
      n_to_process++;
   }
   }

   if (heuristic) {
      // Determine TAD size.
      double total = 0.0;
      for (l = 0 ; l < m ; l++)
      for (j = 0 ; j < n ; j++)
      for (i = j+1 ; i < n ; i++)
         total += obs[l][i+j*n];

      double Q80 = 0.0;
      for (i = 1 ; i < n ; i++) {
         for (l = 0 ; l < m ; l++)
         for (j = 0 ; j < n-i ; j++)
               Q80 += obs[l][j+i+j*n];
         if (Q80 / total > .80) break;
      }

      max_tad_size = i;
      fprintf(stderr, "set 'max_tad_size' to %d\n", i);


      /*
      int *indices = (int *) malloc(n*(n-1)/2 * sizeof(int));
      for (i = 0, l = 0 ; i < n ; i++) {
      for (j = i+1 ; j < n ; j++) {
         indices[l++] = i+j*n;
      }
      }

      globs = obs;
      globw = weights;
      globm = m;
      qsort(indices, n*(n-1)/2, sizeof(int), obscmp);

      double *cumsums = (double *) malloc(n*(n-1)/2 * sizeof(double));
      memset(cumsums, 0.0, n*(n-1)/2 * sizeof(double));
      for (l = 0 ; l < n*(n-1)/2; l++) {
         if (l > 0) cumsums[l] = cumsums[l-1];
         for (k = 0 ; k < m ; k++) {
            cumsums[l] += obs[k][indices[l]] / weights[k][indices[l]];
         }
      }
      double fullsum = cumsums[n*(n-1)/2-1];

      //memset(skip, 1, n*n * sizeof(int));
      //n_to_process = 0;
      max_tad_size = 0;
      for (l = 0 ; l < n*(n-1)/2 ; l++) {
         int i0 = indices[l] % n;
         int j0 = indices[l] / n;
         if ((j0 - i0) > max_tad_size) max_tad_size = j0 - i0;

         for (i = i0 ; i < j0 ; i++) {
         for (j = i+1 ; j < j0 ; j++) {
            if (skip[i+j*n]) n_to_process++;
            skip[i+j*n] = 0;
         }
         }

         if (cumsums[l] / fullsum > 0.9) break;
      }
      printf("max_tad_size set to %d\n", max_tad_size);

      free(cumsums);
      free(indices);
      */

   } // End of heuristic pre-screen.

   for (i = 0 ; i < n-3 ; i++) {
   for (j = i+3 ; j < n ; j++) {
      if ((j-i) > max_tad_size) {
         skip[i+j*n] = 1;
         n_to_process--;
         continue;
      }
   }
   }


   // Allocate 'tid'.
   pthread_t *tid = (pthread_t *) malloc((1+n_threads)*sizeof(pthread_t));

   thread_arg arg = {
      .n = n,
      .m = m,
      .obs = obs,
      .dist = dist,
      .log_gamma = log_gamma,
      .skip = skip,
      .llikmat = llikmat,
      .verbose = verbose,
   };

   errno = pthread_mutex_init(&tadbit_lock, NULL);
   if (errno) {
      fprintf(stderr, "error initializing mutex (%d)\n", errno);
      return;
   }

   taskQ_i = 0;

   memset(tid, 0, (1 + n_threads) * sizeof(pthread_t *));
   for (i = 1 ; i < 1 + n_threads ; i++) {
      errno = pthread_create(&(tid[i]), NULL, &fill_llikmat, &arg);
      if (errno) {
         fprintf(stderr, "error creating thread (%d)\n", errno);
         return;
      }
   }

   // Wait for threads to return.
   for (i = 1 ; i < 1 + n_threads ; i++) {
      pthread_join(tid[i], NULL);
   }
   if (verbose) {
      fprintf(stderr, "computing likelihood (100%% done)\n");
   }

   pthread_mutex_destroy(&tadbit_lock);
   free(skip);
   free(tid);

   // The matrix 'llikmat' now contains the log-likelihood of the
   // segments. The breakpoints are found by dynamic programming.

   const int maxbreaks = n/4;
   double *mllik = (double *) malloc(maxbreaks * sizeof(double));
   int *bkpts = (int *) malloc(n*maxbreaks * sizeof(int));

   mlwalk(llikmat, n, m, maxbreaks, mllik, bkpts);

   // Get optimal number of breaks by AIC.
   double AIC = -INFINITY;
   int n_params;
   int nbreaks_opt;
   for (nbreaks_opt = 1 ; nbreaks_opt  < maxbreaks ; nbreaks_opt++) {
      n_params = nbreaks_opt + m*(8 + nbreaks_opt*6);
      if (AIC > mllik[nbreaks_opt] - n_params) {
         break;
      }
      else {
         AIC = mllik[nbreaks_opt] - n_params;
      }
   }
   nbreaks_opt -= 1;

   // XXX Implementation of the quality score.
   double *llikmatcpy = (double *) malloc (n*n * sizeof(double));
   double *mllikcpy = (double *) malloc(maxbreaks * sizeof(double));
   int *bkptscpy = (int *) malloc(n*maxbreaks * sizeof(int));
   int *passages = (int *) malloc(n * sizeof(int));
   memcpy(llikmatcpy, llikmat, n*n * sizeof(double));
   memcpy(bkptscpy, bkpts, n*maxbreaks * sizeof(int));
   memset(passages, 0, n * sizeof(int));

   for (l = 0 ; l < 10 ; l++) {
      i = 0;
      for (j = 0 ; j < n ; j++) {
         if (bkptscpy[j+nbreaks_opt*n]) {
            llikmatcpy[i+j*n] -= m*6;
            i = j+1;
            passages[j] += bkpts[j+nbreaks_opt*n];
         }
      }
      mlwalk(llikmatcpy, n, m, maxbreaks, mllikcpy, bkptscpy);
   }

   free(llikmatcpy);
   free(mllikcpy);
   free(bkptscpy);

   // Resize output to match original.
   int *resized_bkpts = (int *) malloc(N*maxbreaks * sizeof(int));
   int *resized_passages = (int *) malloc(N * sizeof(int));
   memset(resized_bkpts, 0, N*maxbreaks * sizeof(int));
   memset(resized_passages, 0, N * sizeof(int));

   for (l = 0, i = 0 ; i < N ; i++) {
      if (!remove[i]) {
         resized_passages[i] = passages[l];
         for (j = 0 ; j < maxbreaks ; j++) {
            resized_bkpts[i+j*N] = bkpts[l+j*n];
         }
         l++;
      }
   }

   free(passages);
   free(bkpts);

   double *resized_llikmat = (double *) malloc(N*N * sizeof(double));
   for (i = 0 ; i < N*N ; i++) {
      resized_llikmat[i] = NAN;
   }

   for (l = 0, i = 0 ; i < N ; i++) {
      if (remove[i]) continue;
      for (k = 0, j = 0 ; j < N ; j++) {
         if (remove[j]) continue;
         resized_llikmat[i+j*N] = llikmat[l+k*n];
         k++;
      }
      l++;
   }

   free(llikmat);

   for (k = 0 ; k < m ; k++) {
      free(new_obs[k]);
      free(log_gamma[k]);
      free(weights[k]);
   }
   free(weights);
   free(new_obs);
   free(log_gamma);
   free(dist);

   

   // Update output struct.
   seg->maxbreaks = maxbreaks;
   seg->nbreaks_opt = nbreaks_opt;
   seg->passages = resized_passages;
   seg->llikmat = resized_llikmat;
   seg->mllik = mllik;
   seg->bkpts = resized_bkpts;

   return;

}
