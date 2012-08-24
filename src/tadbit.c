#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include <ctype.h>
#include <unistd.h>

#include "tadbit.h"


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
  const ml_slice *slc
){
// SYNOPSIS:                                                            
//   Wrapper for 'poiss_reg'. Fits the three regions of a slice         
//   by Poisson regression and return the log-likelihood.               
//                                                                      
// PARAMETERS:                                                          
//   '*slc' : the slice to be fitted by Poisson regression.             
//                                                                      
// RETURN:                                                              
//   The total maximum log-likelihood of the slice.                     
//                                                                      

   double top = poiss_reg(slc->blocks[0]);
   double mid = poiss_reg(slc->blocks[1]);
   double bot = poiss_reg(slc->blocks[2]);

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
  ml_slice *slc
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
//   '*slc' : variable to store the slice.                              
//                                                                      
// RETURN:                                                              
//   'void'                                                             
//                                                                      
// SIDE-EFFECTS:                                                        
//   Update 'slc' in place.                                             
//                                                                      

   int l;
   int pos;
   int row;
   int col;

   // Initialize block sizes to 0.
   for (l = 0 ; l < 3 ; l++) {
      slc->blocks[l]->size = 0;
   }

   // Iterate over data vertically.
   for (col = start ; col < end+1 ; col++) {
      for (row = 0 ; row < n ; row++) {

         // Skip NAs in the data.
         if (isnan(obs[row+col*n]))
            continue;


         //TODO: test those lines.

         // Find which block to update ('l').
         //     0: top                       
         //     1: middle                    
         //     2: bottom                    
         
         if        (row < start)  l = 0;
         else if   (row < col)    l = 1;
         else if   (row > end)    l = 2;
         else                     continue;

         pos = slc->blocks[l]->size;

         slc->blocks[l]->counts[pos] = log_gamma[row+col*n];
         slc->blocks[l]->counts[pos] = obs[row+col*n];
         slc->blocks[l]->dist[pos] = dist[row+col*n];
         // The weight is the square root of the product of the
         // diagonal terms (the product of potencies).
         slc->blocks[l]->weights[pos] = \
                  sqrt(obs[row+row*n]*obs[col+col*n]);

         slc->blocks[l]->size++;

      }
   }

   return;

}

void
mlwalk(
  /* input */
  const double *llik,
  const int n,
  const int m,
  /* output */
  double *mllik,
  int *breakpoints
){
// SYNOPSIS:                                                            
//   Dynamic programming algorithm to compute the most likely position  
//   of breakpoints given a matrix of slice maximum log-likelihood.     
//                                                                      
// PARAMETERS:                                                          
//   '*llik' : matrix of maximum log-likelihood values.                 
//   'n' : row/col number of 'llik'.                                    
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

   // The max number of segments is 1/4 row/col number.
   const int max_n_breaks = n/4;

   int i;
   int j;
   int nbreaks;

   double new_llik[n];
   double old_llik[n];

   // Breakpoint lists. The first index (line) is the end of the slice,
   // the second is 1 if this position is an end (breakpoint).
   // These must be allocated from the heap because 'n' can be large.
   int *new_bkpt_list = (int *) malloc(n*n * sizeof(int));
   int *old_bkpt_list = (int *) malloc(n*n * sizeof(int));

   // Initializations.
   memset(breakpoints, 0, n*max_n_breaks);

   for (i = 0 ; i < n*n ; i++) {
      new_bkpt_list[i] = 0;
      old_bkpt_list[i] = 0;
   }

   for (i = 0 ; i < max_n_breaks ; i++) {
      mllik[i] = NAN;
   }

   // Initialize 'old_llik' to the first line of 'llik' containing the
   // log-likelihood of segments starting at index 0.
   for (i = 0 ; i < n ; i++) {
      old_llik[i] = llik[i*n];
      new_llik[i] = -INFINITY;
   }


   // Dynamic programming.
   for (nbreaks = 1 ; nbreaks < max_n_breaks ; nbreaks++) {
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
            double tmp = old_llik[i-1] + llik[i+j*n];
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


void
tadbit(
  /* input */
  double **obs,
  int n,
  const int m,
  double max_tad_size,
  int n_threads,
  const int verbose,
  const int heuristic,
  /* output */
  tadbit_output *seg
){


   // Get absolute max tad size.
   max_tad_size = max_tad_size > 1 ? max_tad_size : max_tad_size * n;

   int i, j, k, l;
   const int init_n = n;

   // Allocate memory and initialize variables. The distance
   // matrix 'dis' is the distance to the main diagonal. Every
   // element of coordinate (i,j) is on a diagonal; the distance
   // is the log-shift to the main diagonal 'i-j'.

   double *init_dis = (double *) malloc(init_n*init_n * sizeof(double));
   double *llikmat = (double *) malloc(init_n*init_n * sizeof(double));

   for (l = 0, i = 0; i < init_n ; i++) {
   for (j = 0; j < init_n ; j++) {
      init_dis[l] = log(abs(i-j));
      llikmat[l] = NAN;
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
   double **log_gamma = (double **) malloc(m * sizeof(double *));
   double **new_obs = (double **) malloc(m * sizeof(double *));
   double *dis = (double *) malloc(n*n * sizeof(double));
   for (k = 0 ; k < m ; k++) {
      l = 0;
      log_gamma[k] = (double *) malloc(n*n * sizeof(double));
      new_obs[k] = (double *) malloc(n*n * sizeof(double));
      for (j = 0 ; j < init_n ; j++) {
      for (i = 0 ; i < init_n ; i++) {
         if (remove[i] || remove[j]) {
            continue;
         }
         log_gamma[k][l] = lgamma(obs[k][i+j*init_n]);
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


   //FIXME: re-implement heurisitc.
   char *to_include = (char*) malloc(n*n * sizeof(char));
   memset(to_include, 1, n*n);


   for (i = 0 ; i < n-3 ; i++) {
   for (j = i+3 ; j < n ; j++) {
      //if (j-i < max_tad_size) {
      if (j-i < max_tad_size && to_include[i+j*n]) {
         assignment[i+j*n] = to_process % n_threads;
         to_process++;
      }
      else {
         // Do not process if slice is larger than 'max_tad_size'.
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
      ml_block *top = (ml_block *) malloc(sizeof(ml_block));
      ml_block *mid = (ml_block *) malloc(sizeof(ml_block));
      ml_block *bot = (ml_block *) malloc(sizeof(ml_block));
      slc->blocks[0] = top;
      slc->blocks[1] = mid;
      slc->blocks[2] = bot;
      // Readability variable.
      int nmax = (n+1)*(n+1)/4;

      // TODO: Reduce max size to minimum actually used.
      // My last attempts caused segmentation faults.
      top->lgamma  = (double *) malloc(nmax * sizeof(double));
      top->counts  = (double *) malloc(nmax * sizeof(double));
      top->dist    = (double *) malloc(nmax * sizeof(double));
      top->weights = (double *) malloc(nmax * sizeof(double));

      mid->lgamma  = (double *) malloc(2 * nmax * sizeof(double));
      mid->counts  = (double *) malloc(2 * nmax * sizeof(double));
      mid->dist    = (double *) malloc(2 * nmax * sizeof(double));
      mid->weights = (double *) malloc(2 * nmax * sizeof(double));

      bot->lgamma  = (double *) malloc(nmax * sizeof(double));
      bot->counts  = (double *) malloc(nmax * sizeof(double));
      bot->dist    = (double *) malloc(nmax * sizeof(double));
      bot->weights = (double *) malloc(nmax * sizeof(double));
      
      
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
            slice(log_gamma[k], obs[k], dis, n, i, j, slc);
            // Get the likelihood and sum (see macro definition).
            llik[i+j*n] += fit_slice(slc);
         }
         processed++;
         if (verbose) {
            fprintf(stderr, "computing likelihood (%0.f%% done)\r",
               99 * processed / (float) to_process);
         }
      }
      } // End of the (i,j) for loop.

      // Free allocated memory.
      for (i = 0 ; i < 3 ; i++) {
         free(slc->blocks[i]->lgamma);
         free(slc->blocks[i]->counts);
         free(slc->blocks[i]->dist);
         free(slc->blocks[i]->weights);
         free(slc->blocks[i]);
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
      fprintf(stderr, "computing likelihood (100%% done)\n");
   }

   free(assignment);

   // The matrix 'llik' now contains the log-likelihood of the
   // segments. The breakpoints are found by dynamic
   // programming.

   //int *all_breakpoints = (int *) malloc(n*n * sizeof(int));
   double *mllik = (double *) malloc(n/4 * sizeof(double));
   int *bkpts = (int *) malloc(n*n/4 * sizeof(int));

   mlwalk(llik, n, m, mllik, bkpts);

   // Resize output to match original.
   int *resized_bkpts = (int *) malloc(init_n*n/4 * sizeof(int));
   memset(resized_bkpts, 0, init_n*n/4);

   for (l = 0, i = 0 ; i < init_n ; i++) {
      if (!remove[i]) {
         for (j = 0 ; j < n/4 ; j++) {
            resized_bkpts[i+j*init_n] = bkpts[l+j*n];
         }
         l++;
      }
   }

   free(bkpts);

   for (l = 0, i = 0 ; i < init_n ; i++) {
      if (remove[i]) continue;
      for (k = 0, j = 0 ; j < init_n ; j++) {
         if (remove[j]) continue;
         llikmat[i+j*init_n] = llik[l+k*n];
         k++;
      }
      l++;
   }

   for (k = 0 ; k < m ; k++) {
      free(new_obs[k]);
   }
   free(to_include);
   free(new_obs);
   free(dis);
   free(llik);

   // Update output struct.
   seg->nbreaks_opt = 15;
   seg->llikmat = llikmat;
   seg->mllik = mllik;
   seg->bkpts = resized_bkpts;

   return;

}

void
free_tadbit_output(
  tadbit_output *seg
){
  free(seg->llikmat);
  free(seg->mllik);
  free(seg->bkpts);
  free(seg);
}
