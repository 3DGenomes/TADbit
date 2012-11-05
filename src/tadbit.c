#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include <ctype.h>
#include <unistd.h>

#include "tadbit.h"


// Globals variables. //

int *_cache_index;
int n_processed;              // Number of slices processed so far.
int n_to_process;             // Total number of slices to process.
int taskQ_i;                  // Index used for task queue.
pthread_mutex_t tadbit_lock;  // Mutex to access task queue.


void
calcfg(
  /* input */
  const int    n,
  const int    i_,
  const int    _i,
  const int    j_,
  const int    _j,
  const int    diag,
  const double *k,
  const double *d,
  const double *w,
  const double a,
  const double b,
  const double da,
  const double db,
        double *c,
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
// TODO Update arguments.
//                                                                      
// RETURN:                                                              
//   'void'                                                             
//                                                                      
// SIDE-EFFECTS:                                                        
//   Update 'f' and 'g' in place.                                       
//                                                                      

   // 'tmp' is a computation intermediate that will be the return
   // value of 'exp'. This can call '__slowexp' which on 64-bit machines
   // can return a long double (causing segmentation fault if 'tmp' is
   // declared as long).
   long double tmp;
   int i;
   int j;

   *f = 0.0; *g = 0.0;
   for (i = 0 ; i < n ; i++) c[i] = NAN;
   // Only the boundaries change for the diagonal half-block.
   // The computations are the same in both cases.
   if (diag) {
      for (j = j_+1 ; j < _j+1 ; j++) {
      for (i = i_ ; i < j ; i++) {
         // Retrive value of the exponential from cache.
         if (isnan(c[_cache_index[i+j*n]])) {
            c[_cache_index[i+j*n]] = exp(a+da+(b+db)*d[i+j*n]);
         }
         tmp  =  w[i+j*n] * c[_cache_index[i+j*n]] - k[i+j*n];
         *f  +=  tmp;
         *g  +=  tmp * d[i+j*n];
      }
      }
   }
   else {
      for (j = j_ ; j < _j+1 ; j++) {
      for (i = i_ ; i < _i+1 ; i++) {
         // Retrive value of the exponential from cache.
         if (isnan(c[_cache_index[i+j*n]])) {
            c[_cache_index[i+j*n]] = exp(a+da+(b+db)*d[i+j*n]);
         }
         tmp  =  w[i+j*n] * c[_cache_index[i+j*n]] - k[i+j*n];
         *f  +=  tmp;
         *g  +=  tmp * d[i+j*n];
      }
      }
   }

   return;

}

double
ll(
  const int    n,
  const int    i_,
  const int    _i,
  const int    j_,
  const int    _j,
  const int    diag,
  const double *k,
  const double *d,
  const double *w,
  const double *lg,
        double *c
){
// SYNOPSIS:                                                            
//   The fitted model (by maximum likelihood) is Poisson with lambda    
//   paramter such that lambda = w * exp(a + b*d). So the full          
//   log-likelihood of the model is the sum of terms                    
//                                                                      
//      - w_i exp(a + b*d_i) + k_i(log(w_i) + a + b*d_i) - log(k_i!)    
//                                                                      
// ARGUMENTS:                                                           
// TODO: Describe the arguments.
//                                                                      
// RETURN:                                                              
//   The maximum log-likelihood of a block of hiC data.                 
//                                                                      

   if ((i_ >= _i) || (j_ >= _j)) return 0.0;

   int i;
   int j;
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
   // See the comment about 'tmp' in 'calcfg'.
   long double tmp; 

   calcfg(n, i_, _i, j_, _j, diag, k, d, w, a, b, da, db, c, &f, &g);

   // Newton-Raphson until gradient function is less than TOLERANCE.
   // The gradient function is the square norm 'f*f + g*g'.
   while ((oldgrad = f*f + g*g) > TOLERANCE && iter++ < MAXITER) {

      for (i = 0 ; i < n ; i++) c[i] = NAN;
      // Compute the derivatives.
      dfda = dfdb = dgda = dgdb = 0.0;
      // Only the boundaries change for the diagonal half-block.
      // The computations are the same in both cases.
      if (diag) {
         for (j = j_+1 ; j < _j+1 ; j++) {
         for (i = i_ ; i < j ; i++) {
            // Retrive value of the exponential from cache.
            if (isnan(c[_cache_index[i+j*n]])) {
               c[_cache_index[i+j*n]] = exp(a+b*d[i+j*n]);
            }
            tmp   =   w[i+j*n] * c[_cache_index[i+j*n]];
            dfda +=   tmp;
            dgda +=   tmp * d[i+j*n];
            dgdb +=   tmp * d[i+j*n]*d[i+j*n];
         }
         }
      }
      else {
         for (j = j_ ; j < _j+1 ; j++) {
         for (i = i_ ; i < _i+1 ; i++) {
            // Retrive value of the exponential from cache.
            if (isnan(c[_cache_index[i+j*n]])) {
               c[_cache_index[i+j*n]] = exp(a+b*d[i+j*n]);
            }
            tmp   =   w[i+j*n] * c[_cache_index[i+j*n]];
            dfda +=   tmp;
            dgda +=   tmp * d[i+j*n];
            dgdb +=   tmp * d[i+j*n]*d[i+j*n];
         }
         }
      }
      dfdb = dgda;

      denom = dfdb*dgda - dfda*dgdb;
      da = (f*dgdb - g*dfdb) / denom;
      db = (g*dfda - f*dgda) / denom;

      calcfg(n, i_, _i, j_, _j, diag, k, d, w, a, b, da, db, c, &f, &g);

      // Traceback if we are not going down the gradient. Cut the
      // length of the steps in half until this step goes down
      // the gradient.
      while (f*f + g*g > oldgrad) {
         da /= 2;
         db /= 2;
         calcfg(n, i_, _i, j_, _j, diag, k, d, w, a, b, da, db, c, &f, &g);
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
   if (diag) {
      for (j = j_+1 ; j < _j+1 ; j++) {
      for (i = i_ ; i < j ; i++) {
         llik += exp(a+b*d[i+j*n]) + k[i+j*n]*(a+b*d[i+j*n]) - lg[i+j*n];
      }
      }
   } 
   else {
      for (j = j_ ; j < _j+1 ; j++) {
      for (i = i_ ; i < _i+1 ; i++) {
         llik += exp(a+b*d[i+j*n]) + k[i+j*n]*(a+b*d[i+j*n]) - lg[i+j*n];
      }
      }
   } 

   return llik;

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
// TODO Describe the parameters.
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
   const double **k = (const double **) myargs->k;
   const double *d = (const double*) myargs->d;
   const double **w = (const double **) myargs->w;
   const double **lg= (const double **) myargs->lg;
   const int *skip = (const int *) myargs->skip;
   double *llikmat = myargs->llikmat;
   const int verbose = myargs->verbose;

   int i;
   int j;
   int l;

   // Cache to speed up computation.
   double *c= (double *) malloc(n * sizeof(double));
   memset(c, 0.0, n * sizeof(double));

   int job_index;
   
   // Break out of the loop when task queue is empty.
   while (1) {

      pthread_mutex_lock(&tadbit_lock);
      while ((taskQ_i < n*n) && (skip[taskQ_i] > 0)) {
         // Fast forward to the next job.
         taskQ_i++;
      }
      if (taskQ_i >= n*n) {
         // Task queue is empty. Exit loop and return
         pthread_mutex_unlock(&tadbit_lock);
         break;
      }
      job_index = taskQ_i;
      taskQ_i++;
      pthread_mutex_unlock(&tadbit_lock);

      // Compute the log-likelihood of slice '(i,j)'.
      i = job_index % n;
      j = job_index / n;

      // Distinct parts of the array, no lock needed.
      llikmat[i+j*n] = 0.0;
      for (l = 0 ; l < m ; l++) {
         llikmat[i+j*n] += 
            ll(n, 0, i-1, i, j, 0, k[l], d, w[l], lg[l], c) / 2 +
            ll(n, i,   j, i, j, 1, k[l], d, w[l], lg[l], c) +
            ll(n, j+1, n, i, j, 0, k[l], d, w[l], lg[l], c) / 2;
      }
      n_processed++;
      if (verbose) {
         fprintf(stderr, "computing likelihood (%0.f%% done)\r",
            99 * n_processed / (float) n_to_process);
      }
   }

   free(c);
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

   _cache_index = (int *) malloc(n*n * sizeof(int));
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
         _cache_index[l] = i > j ? i-j : j-i;
         log_gamma[k][l] = lgamma(obs[k][i+j*N]+1);
         new_obs[k][l] = obs[k][i+j*N];
         dist[l] = init_dist[i+j*N];
         l++;
      }
      }
   }
   fprintf(stderr, "\n");

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
   for (l = 0 ; l < m ; l++) {
      weights[l] = (double *) malloc(n*n * sizeof(double));
      memset(weights[l], 0.0, n*n * sizeof(double));
      // Compute scalar product.
      for (j = 0 ; j < n ; j++) {
      for (i = 0 ; i < n ; i++) {
         // TODO implement the following weights.
         // weights[l][i+j*n] = sqrt(rowsums[l][i]*rowsums[l][j]);
         weights[l][i+j*n] = sqrt(obs[l][i+i*n]*obs[l][j+j*n]);
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

   } // End of heuristic 'max_tad_size' definition.

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
      .k = obs,
      .d = dist,
      .w = weights,
      .lg = log_gamma,
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
