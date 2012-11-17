#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <pthread.h>
#include <ctype.h>
#include <unistd.h>
#include <float.h>

#include "tadbit.h"


// Global variables. //

int *_cache_index;
int _max_cache_index;
int n_processed;              // Number of slices processed so far.
int n_to_process;             // Total number of slices to process.
int taskQ_i;                  // Index used for task queue.
pthread_mutex_t tadbit_lock;  // Mutex to access task queue.


void
fg(
  // input //
  const int    n,
  const int    i_,
  const int    _i,
  const int    j_,
  const int    _j,
  const int    diag,
  const int *k,
  const double *d,
  const double *w,
  const double a,
  const double b,
  const double da,
  const double db,
        double *c,
  // output //
        double *f,
        double *g
){
// SYNOPSIS:                                                            
//   Subfroutine of 'll' that computes 'f' and 'g' for Newton-Raphson   
//   cycles.                                                            
//                                                                      
// ARGUMENTS:                                                           
//   See the function 'll' for the description of 'n', 'i_', '_i',      
//      'j_', '_j', 'diag', 'k', 'd', and 'w'.                          
//   'a': parameter 'a' of the Poisson regression (see 'poiss_reg').    
//   'b': parameter 'b' of the Poisson regression (see 'poiss_reg').    
//   'da': computed differential of 'a' (see 'poiss_reg').              
//   'db': computed differential of 'b' (see 'poiss_reg').              
//        -- output arguments --                                        
//   'f': first function to zero, recomputed by the routine.            
//   'g': second function to zero, recomputed by the routine.           
//   'c': address of an array of double for caching.                    
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
   int i_low = i_;
   int i_high = -1;
   int j_low = diag ? j_+1 : j_;
   int j_high = _j+1;
   int index;

   *f = 0.0; *g = 0.0;
   // Initialize cache.
   for (index = 0 ; index < _max_cache_index ; index++) c[index] = NAN;

   for (j = j_low ; j < j_high ; j++) {
      i_high = diag ? j : _i+1;
      for (i = i_low ; i < i_high ; i++) {
         // Retrieve value of the exponential from cache.
         index = _cache_index[i+j*n];
         if (c[index] != c[index]) {
            c[index] = exp(a+da+(b+db)*d[i+j*n]);
         }
         tmp  =  w[i+j*n] * c[index] - k[i+j*n];
         *f  +=  tmp;
         *g  +=  tmp * d[i+j*n];
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
  const int    *k,
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
//   'n': row/column number of the counts.                              
//   'i_': first value of index i (row).                                
//   '_i': last value of index i (row).                                 
//   'j_': first value of index j (column).                             
//   '_j': last value of index j (column).                              
//   'diag': whether the block is half-diagonal (middle block).         
//   'k': raw hiC counts.                                               
//   'd': distances from diagonal in log.                               
//   'w': weights measuring hiC bias.                                   
//   'lg': log-gamma terms.                                             
//   'c': address of an array of double for caching.                    
//                                                                      
// RETURN:                                                              
//   The maximum log-likelihood of a block of hiC data.                 
//                                                                      

   // For slices at the border of the hiC matrix, the top or bottom
   // blocks have 0 height. Returning 0.0 makes the summation at
   // the line labelled "slice ll summation" still valid.
   if ((i_ >= _i) || (j_ >= _j)) return 0.0;
   // For slices of length 2, the diagonal block has only 1 value,
   // which creates an infinite loop (because there are two parameters
   // to fit). Return NAN because estimation is impossible.
   if ((_i < i_+2) || (_j < j_+2)) return NAN;

   int i;
   int j;
   int i_low = i_;
   int i_high = -1;
   int j_low = diag ? j_+1 : j_;
   int j_high = _j+1;
   int index;
   int iter = 0;
   double denom;
   double oldgrad;
   double f = INFINITY;
   double g = INFINITY;
   double a = 0.0;
   double b = 0.0;
   double da = 0.0;
   double db = 0.0;
   double dfda = 0.0;
   double dfdb = 0.0;
   double dgda = 0.0;
   double dgdb = 0.0;
   // See the comment about 'tmp' in 'fg'.
   long double tmp; 

   fg(n, i_, _i, j_, _j, diag, k, d, w, a, b, da, db, c, &f, &g);

   // Newton-Raphson until gradient function is less than TOLERANCE.
   // The gradient function is the square norm 'f*f + g*g'.
   while ((oldgrad = f*f + g*g) > TOLERANCE && iter++ < MAXITER) {

      for (index = 0 ; index < _max_cache_index ; index++) c[index] = NAN;
      // Compute the derivatives.
      dfda = dfdb = dgda = dgdb = 0.0;

      for (j = j_low ; j < j_high ; j++) {
         i_high = diag ? j : _i+1;
         for (i = i_low ; i < i_high ; i++) {
            index = _cache_index[i+j*n];
            // Retrive value of the exponential from cache.
            if (c[index] != c[index]) { // ERROR.
               c[index] = exp(a+b*d[i+j*n]);
            }
            tmp   =   w[i+j*n] * c[index];
            dfda +=   tmp;
            tmp  *=   d[i+j*n];
            dgda +=   tmp;
            tmp  *=   d[i+j*n];
            dgdb +=   tmp;
         }
      }
      dfdb = dgda;

      denom = dfdb*dgda - dfda*dgdb;
      da = (f*dgdb - g*dfdb) / denom;
      db = (g*dfda - f*dgda) / denom;

      fg(n, i_, _i, j_, _j, diag, k, d, w, a, b, da, db, c, &f, &g);

      // Traceback if we are not going down the gradient. Cut the
      // length of the steps in half until this step goes down
      // the gradient.
      for (i = 0 ; (i < 20) && (f*f + g*g > oldgrad) ; i++) {
         da /= 2;
         db /= 2;
         fg(n, i_, _i, j_, _j, diag, k, d, w, a, b, da, db, c, &f, &g);
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

   //////////////////////////////////////////////////////////////////////
   // The last call to 'fg' has set the cache to the right values. //
   // No need to reset the cache.                                      //
   //////////////////////////////////////////////////////////////////////
   for (j = j_low ; j < j_high ; j++) {
      i_high = diag ? j : _i+1;
      for (i = i_low ; i < i_high ; i++) {
         index = _cache_index[i+j*n];
         // Retrive value of the exponential from cache.
         llik += c[index] + k[i+j*n]*(a+b*d[i+j*n]) - lg[i+j*n];
      }
   }

   return llik;

}

void *
fill_DP(
  void *arg
){
// SYNOPSIS:                                                            
//   Thread function to compute the values of the dynamic programming   
//   'DPwalk' for a given number of breaks (value of 'nbreaks').        
//                                                                      
// PARAMETERS:                                                          
//   'arg': arguments for a thread (see header file).                   
//                                                                      
// RETURN:                                                              
//   'void *'                                                           
//                                                                      
// SIDE-EFFECTS:                                                        
//   Update 'old_llik', 'new_llik' and 'new_bkpt_list' in place.        
//                                                                      

   dpworker_arg *myargs = (dpworker_arg *) arg;
   const int n = myargs->n;
   const double *llikmat = (const double *) myargs->llikmat;
   double *old_llik = (double *) myargs->old_llik;
   double *new_llik = (double *) myargs->new_llik;
   const int nbreaks = myargs->nbreaks;
   int *new_bkpt_list = (int *) myargs->new_bkpt_list;
   const int *old_bkpt_list = (const int *) myargs->old_bkpt_list;

   int i;

   while (1) {
      pthread_mutex_lock(&tadbit_lock);
      if (taskQ_i > n-1) {
         // Task queue is empty. Exit loop and return
         pthread_mutex_unlock(&tadbit_lock);
         break;
      }
      // A task gives an end point 'j'.
      int j = taskQ_i;
      taskQ_i++;
      pthread_mutex_unlock(&tadbit_lock);

      new_llik[j] = -INFINITY;
      int new_bkpt = -1;

      // Cycle over start point 'i'.
      for (i = 3 * nbreaks ; i < j-3 ; i++) {

         // If NAN the following condition evaluates to false.
         double tmp = old_llik[i-1] + llikmat[i+j*n];
         if (tmp > new_llik[j]) {
            new_llik[j] = tmp;
            new_bkpt = i-1;
         }
      }

      // Update breakpoint list (skip if log-lik is undefined).
      // No need to use mutex because 'j' is different for every thread.
      if (new_llik[j] > -INFINITY) {
         for (i = 0 ; i < n ; i++) {
            new_bkpt_list[j+i*n] = old_bkpt_list[new_bkpt+i*n];
         }
         new_bkpt_list[j+new_bkpt*n] = 1;
      }
   }

   return NULL;

}

void
DPwalk(
  // input //
  const double *llikmat,
  const int n,
  const int MAXBREAKS,
  int n_threads,
  // output //
  double *mllik,
  int *breakpoints
){
// SYNOPSIS:                                                            
//   Dynamic programming algorithm to compute the most likely position  
//   of breakpoints given a matrix of slice maximum log-likelihood.     
//                                                                      
// PARAMETERS:                                                          
//   '*llikmat': matrix of maximum log-likelihood values.               
//   'n': row/col number of 'llikmat'.                                  
//   'MAXBREAKS': The maximum number of breakpoints.                    
//        -- output arguments --                                        
//   '*mllik': maximum log-likelihood of the segmentations.             
//   '*breakpoints': optimal breakpoints per number of breaks.          
//                                                                      
// RETURN:                                                              
//   'void'                                                             
//                                                                      
// SIDE-EFFECTS:                                                        
//   Update 'breakpoints' in place.                                     
//                                                                      

   int i;
   int nbreaks;

   double new_llik[n];
   double old_llik[n];

   // Breakpoint lists. The first index (row) is the end of the slice,
   // the second is 1 if this position is an end (breakpoint).
   // These must be allocated from the heap because 'n' can be large.
   int *new_bkpt_list = (int *) malloc(n*n * sizeof(int));
   int *old_bkpt_list = (int *) malloc(n*n * sizeof(int));

   // Initializations.
   // 'breakpoints' is a 'n' x 'MAXBREAKS' array. The first index (row)
   // is 1 if there is a breakpoint at that location, the second index
   // (column) is the number of breakpoints.
   for (i = 0 ; i < n*MAXBREAKS ; i++) breakpoints[i] = 0;

   for (i = 0 ; i < n*n ; i++) {
      new_bkpt_list[i] = 0;
      old_bkpt_list[i] = 0;
   }

   for (i = 0 ; i < MAXBREAKS ; i++) {
      mllik[i] = NAN;
   }

   // Initialize 'old_llik' to the first line of 'llikmat' containing
   // the log-likelihood of segments starting at index 0.
   for (i = 0 ; i < n ; i++) {
      old_llik[i] = llikmat[i*n];
      new_llik[i] = -INFINITY;
   }

   int err = pthread_mutex_init(&tadbit_lock, NULL);
   if (err) {
      fprintf(stderr, "error initializing mutex (%d)\n", err);
      return;
   }

   dpworker_arg arg = {
      .n = n,
      .llikmat = llikmat,
      .old_llik = old_llik,
      .new_llik = new_llik,
      .nbreaks = 1,
      .new_bkpt_list = new_bkpt_list,
      .old_bkpt_list = old_bkpt_list,
   };

   pthread_t *tid = (pthread_t *) malloc(n_threads * sizeof(pthread_t));
   for (i = 0 ; i < n_threads ; i++) tid[i] = 0;
   for (i = 0 ; i < n_threads ; i++) {
      err = pthread_create(&(tid[i]), NULL, &fill_DP, &arg);
      if (err) {
         fprintf(stderr, "error creating thread (%d)\n", err);
         return;
      }
   }

   // Dynamic programming.
   for (nbreaks = 1 ; nbreaks < MAXBREAKS ; nbreaks++) {

      arg.nbreaks = nbreaks;
      // Update breakpoint lists.
      for (i = 0 ; i < n*n ; i++) {
         old_bkpt_list[i] = new_bkpt_list[i];
      }
      taskQ_i = 3 * nbreaks + 2;

      for (i = 0 ; i < n_threads ; i++) tid[i] = 0;
      for (i = 0 ; i < n_threads ; i++) {
         err = pthread_create(&(tid[i]), NULL, &fill_DP, &arg);
         if (err) {
            fprintf(stderr, "error creating thread (%d)\n", err);
            return;
         }
      }

      // Wait for threads to return.
      for (i = 0 ; i < n_threads ; i++) {
         pthread_join(tid[i], NULL);
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
//   'arg': thread arguments (see header file for definition).          
//                                                                      
// RETURN:                                                              
//   'void'                                                             
//                                                                      
// SIDE-EFFECTS:                                                        
//   Update 'llikmat' in place.                                         
//                                                                      

   llworker_arg *myargs = (llworker_arg *) arg;
   const int n = myargs->n;
   const int m = myargs->m;
   const int **k = (const int **) myargs->k;
   const double *d = (const double*) myargs->d;
   const double **w = (const double **) myargs->w;
   const double **lg= (const double **) myargs->lg;
   const char *skip = (const char *) myargs->skip;
   double *llikmat = myargs->llikmat;
   const int verbose = myargs->verbose;

   int i;
   int j;
   int l;

   // Cache to speed up computation. Get the max of '_cache_index'
   // in order to allocate the right size.
   double *c= (double *) malloc(_max_cache_index * sizeof(double));
   for (i = 0 ; i < _max_cache_index ; i++) c[i] = 0.0;

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
         // This is the slice ll summation.
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
allocate_heur_job(
  char *skip,
  const int i0,
  const int j0,
  const int n
){
// SYNOPSIS:                                                            
//   Create or update thread jobs (used in pre-heuristic).
//                                                                      
// PARAMETERS:                                                          
//   'skip': the job matrix to update in place.                         
//   'i0': start position of the approximate TAD.                       
//   'j0': end position of the approximate TAD.                         
//   'n': number of rows/columns of the hiC matrix (or 'skip').         
//                                                                      
// RETURN:                                                              
//   'void'                                                             
//                                                                      
// SIDE-EFFECTS:                                                        
//   Update 'skip' in place.                                            
//                                                                      

   int i;
   int j;

   for (j = j0-2 ; j < j0+3 ; j++)
   for (i = i0-2 ; i < i0+3 ; i++)
      if ((i+j*n > 0) && (i+j*n < n*n)) skip[i+j*n] = 0;

}

void
allocate_new_jobs(
  char *skip,
  const int *bkpts,
  const int MAXBREAKS,
  const int nbreaks_opt,
  const int n
){
// SYNOPSIS:                                                            
//   Create or update thread jobs. For an approximate TAD defined by    
//   ('i0', 'j0'), create jobs that will compute the log-likelihood of  
//   this TAD, and all single splits of this TAD.                       
//                                                                      
// PARAMETERS:                                                          
// TODO Update parameters
//   'skip': the job matrix to update in place.                         
//   'n': number of rows/columns of the hiC matrix (or 'skip').         
//                                                                      
// RETURN:                                                              
//   'void'                                                             
//                                                                      
// SIDE-EFFECTS:                                                        
//   Update 'skip' in place.                                            
//                                                                      

   int i;
   int j;
   int i0;
   int j0;
   int shift;
   char *starts = (char*) malloc(n* sizeof(char));
   char *ends = (char*) malloc(n* sizeof(char));
   for (i = 0 ; i < n ; i++) {
      starts[i] = 0;
      ends[i] = 0;
   }
   
   for (shift = -10 ; shift < 11 ; shift++) {
      for (i0 = 0, j0 = 0 ; j0 < n ; j0++) {
         if (bkpts[j0+(shift+nbreaks_opt)*n]) {

            // Jobs for splitting the TAD.
            for (j = i0 ; j < j0 ; j++)
               skip[i0+j*n] = 0;
            for (i = i0 ; i < j0 ; i++)
               skip[i+j0*n] = 0;

            starts[i0] = 1;
            ends[j0] = 1;
            i0 = j0+1;
         }
      }
   }

   // Jobs for merging the TADs.
   for (i = 0 ; i < n ; i++)
   for (j = 0 ; j < n ; j++)
      if (starts[i] && ends[j] && (j-i < 500) && (i < j))
         skip[i+j*n] = 0;

   free(starts);
   free(ends);

}

void
tadbit(
  // input //
  int **obs,
  int n,
  const int m,
  int n_threads,
  const int verbose,
  int max_tad_size,
  const int do_not_use_heuristic,
  // output //
  tadbit_output *seg
){

   // Get thread number if set to 0 (automatic).
   if (n_threads < 1) {
      #ifdef _SC_NPROCESSORS_ONLN
         n_threads = (int) sysconf(_SC_NPROCESSORS_ONLN);
      #else
         n_threads = 1;
      #endif
   }

   const int N = n; // Original size.
   int err;       // Used for error checking.

   int i;
   int j;
   int k;
   int l;
   int i0;


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

   const int MAXBREAKS = n/5;

   _cache_index = (int *) malloc(n*n * sizeof(int));
   _max_cache_index = N+1;
   // Allocate and copy.
   double **log_gamma = (double **) malloc(m * sizeof(double *));
   int **new_obs = (int **) malloc(m * sizeof(int *));
   double *dist = (double *) malloc(n*n * sizeof(double));
   for (k = 0 ; k < m ; k++) {
      l = 0;
      log_gamma[k] = (double *) malloc(n*n * sizeof(double));
      new_obs[k] = (int *) malloc(n*n * sizeof(int));
      for (j = 0 ; j < N ; j++) {
      for (i = 0 ; i < N ; i++) {
         if (remove[i] || remove[j]) continue;
         _cache_index[l] = i > j ? i-j : j-i;
         log_gamma[k][l] = lgamma(obs[k][i+j*N]+1);
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
      for (i = 0 ; i < n ; i++) rowsums[k][i] = 0.0;
   }
   for (i = 0 ; i < n ; i++)
   for (k = 0 ; k < m ; k++)
   for (l = 0 ; l < n ; l++)
      rowsums[k][i] += obs[k][i+l*n];

   // Compute the weights.
   double **weights = (double **) malloc(m * sizeof(double *));
   for (l = 0 ; l < m ; l++) {
      weights[l] = (double *) malloc(n*n * sizeof(double));
      for (i = 0 ; i < n*n ; i++) weights[l][i] = 0.0;
      // Compute scalar product.
      for (j = 0 ; j < n ; j++)
      for (i = 0 ; i < n ; i++)
         weights[l][i+j*n] = sqrt(rowsums[l][i]*rowsums[l][j]);
   }

   // We don't need the row/column sums any more.
   for (l = 0 ; l < m ; l++) free(rowsums[l]);
   free(rowsums);

   double *mllik = (double *) malloc(MAXBREAKS * sizeof(double));
   int *bkpts = (int *) malloc(MAXBREAKS*n * sizeof(int));
   double *llikmat = (double *) malloc(n*n * sizeof(double));
   for (i = 0 ; i < n*n ; i++)
      llikmat[i] = NAN;

   // 'skip' will contain only 0 or 1 and can be stored as 'char'.
   char *skip = (char *) malloc(n*n * sizeof(char));

   // Use the heuristic by default (hence the name of the parameter).
   // The parameter 'max_tad_size' is needed only in case the heuristic
   // is not used.
   if (do_not_use_heuristic) {
      for (j = 0 ; j < n ; j++)
      for (i = 0 ; i < n ; i++)
         skip[i+j*n] = (j-i) > max_tad_size ? 1 : 0;
   }
   else {
      if (verbose) {
         fprintf(stderr, "running pre-heuristic\n");
      }

      for (i = 0 ; i < n*n ; i++) skip[i] = 1;

      // 'S[i+j*n]' is the weighted sum of reads within the triangle
      // defined by ('i','j') in the upper triangular matrix of
      // observations.
      double *S = (double *) malloc(n*n * sizeof(double));
      for (i = 0 ; i < n*n ; i++) S[i] = 0.0;
      for (j = 1 ; j < n ; j++) {
      for (i = 0 ; i < n-j ; i++) {
         double weighted_value = 0.0;
         for (l = 0 ; l < m ; l++)
            weighted_value += obs[l][i+(i+j)*n]/weights[l][i+(i+j)*n];
         S[i+(i+j)*n] = S[i+(i+j-1)*n] + S[i+1+(i+j)*n] -
            S[i+1+(i+j-1)*n] + weighted_value;
      }
      }

      double *heur_score = (double *) malloc(n*n * sizeof(double));
      for (i = 0 ; i < n*n ; i++) heur_score[i] = NAN;
      for (j = 1 ; j < n ; j++)
      for (i = 0 ; i < j ; i++)
        heur_score[i+j*n] = log(S[i+j*n]);

      // Use dynamic programming to find approximate break points.
      // The matrix 'mllik' is used only to make the function call valid
      // (it is updated in place, but the value is disregarded), and
      // the heuristic score 'heur_score' plays the role of the
      // log-likelihood 'llikmat'.
      DPwalk(heur_score, n, MAXBREAKS, n_threads, mllik, bkpts);

      free(heur_score);
      free(S);

      // Create a thread job for each approximate TAD.
      for (i = 0 ; i < n*n ; i++) skip[i] = 1;
      for (j = 1 ; j < MAXBREAKS ; j++) {
         i0 = 0;
         for (i = 0 ; i < n ; i++) {
            if (bkpts[i+j*n]) {
               allocate_heur_job(skip, i0, i, n);
               i0 = i+1;
            }
         }
      }

      // Erase the lower triangular part of 'skip'.
      for (j = 0 ; j < n ; j++)
      for (i = j ; i < n ; i++)
         skip[i+j*n] = 1;

      // Allocate estimation of the log likelihood for all small
      // TADs (less than 3 bins).
      for (j = 6 ; j < n ; j++)
      for (i = j-6 ; i < j-3 ; i++)
         skip[i+j*n] = 0;

      // Allocate jobs at the ends of the chromosomes/units because
      // these regions are a bit noisier.
      for (j = 1 ; j < 51 ; j++)
      for (i = 0 ; i < j-3 ; i++)
         skip[i+j*n] = 0;
      for (j = n-51 ; j < n ; j++)
      for (i = n-51 ; i < j-3 ; i++)
         skip[i+j*n] = 0;

      // Reset lower triangular part of 'skip'.
      for (j = 0 ; j < n ; j++)
      for (i = j ; i < n ; i++)
         skip[i+j*n] = 1;

   } // End of pre-heuristic.


   // Allocate 'tid'.
   pthread_t *tid = (pthread_t *) malloc(n_threads * sizeof(pthread_t));

   llworker_arg arg = {
      .n = n,
      .m = m,
      .k = (const int **) obs,
      .d = dist,
      .w = (const double **) weights,
      .lg = (const double **) log_gamma,
      .skip = skip,
      .llikmat = llikmat,
      .verbose = verbose,
   };

   err = pthread_mutex_init(&tadbit_lock, NULL);
   if (err) {
      fprintf(stderr, "error initializing mutex (%d)\n", err);
      return;
   }

   int n_params;
   int nbreaks_opt = 0;
   double AIC = -INFINITY;
   double newAIC = -DBL_MAX;

   while (newAIC > AIC) {

      if (verbose) {
         fprintf(stderr, "starting new cycle\n");
      }

      AIC = newAIC;

      // Initialize task queue.
      n_to_process = 0;
      for (i = 0 ; i < n*n ; i++) {
         // Skip all computation done in previous cycles.
         if (!isnan(llikmat[i])) skip[i] = 1;
         n_to_process += (1-skip[i]);
      }
      n_processed = 0;
      taskQ_i = 0;
      
      // Instantiate threads and start running jobs.
      for (i = 0 ; i < n_threads ; i++) tid[i] = 0;
      for (i = 0 ; i < n_threads ; i++) {
         err = pthread_create(&(tid[i]), NULL, &fill_llikmat, &arg);
         if (err) {
            fprintf(stderr, "error creating thread (%d)\n", err);
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

      // The matrix 'llikmat' now contains the log-likelihood of the
      // segments. The breakpoints are found by dynamic programming.
      int maxbreaks = nbreaks_opt ? nbreaks_opt + 11 : MAXBREAKS;
      DPwalk(llikmat, n, maxbreaks, n_threads, mllik, bkpts);

      // Get optimal number of breaks by AIC.
      newAIC = -INFINITY;
      for (nbreaks_opt = 1 ; nbreaks_opt < MAXBREAKS ; nbreaks_opt++) {
         n_params = nbreaks_opt + m*(8 + nbreaks_opt*6);
         if (newAIC > mllik[nbreaks_opt] - n_params) {
            break;
         }
         else {
            newAIC = mllik[nbreaks_opt] - n_params;
         }
      }
      nbreaks_opt -= 1;

      allocate_new_jobs(skip, bkpts, MAXBREAKS, nbreaks_opt, n);

   }

   AIC = newAIC;

   pthread_mutex_destroy(&tadbit_lock);
   free(skip);
   free(tid);
   free(_cache_index);

   // Compute breakpoint confidence by penalized dynamic progamming.
   double *llikmatcpy = (double *) malloc (n*n * sizeof(double));
   double *mllikcpy = (double *) malloc(MAXBREAKS * sizeof(double));
   int *bkptscpy = (int *) malloc(n*MAXBREAKS * sizeof(int));
   int *passages = (int *) malloc(n * sizeof(int));
   for (i = 0 ; i < n ; i++) passages[i] = 0;
   for (i = 0 ; i < n*MAXBREAKS ; i++) bkptscpy[i] = bkpts[i];
   for (i = 0 ; i < n*n ; i++) llikmatcpy[i] = llikmat[i];

   for (l = 0 ; l < 10 ; l++) {
      i = 0;
      for (j = 0 ; j < n ; j++) {
         if (bkptscpy[j+nbreaks_opt*n]) {
            // Apply a constant penalty every time a TAD is present
            // in the final decomposition. The penalty is set to
            // '6*m' because it is the expected log-likelihood gain
            // for adding a new TAD around the optimum log-likelihood.
            llikmatcpy[i+j*n] -= m*6;
            passages[j] += bkpts[j+nbreaks_opt*n];
            i = j+1;
         }
      }
      if (i < n) llikmatcpy[i+(n-1)*n] -= m*6;
      DPwalk(llikmatcpy, n, nbreaks_opt+1, n_threads, mllikcpy, bkptscpy);
   }
   free(llikmatcpy);
   free(mllikcpy);
   free(bkptscpy);
   
   // Resize output to match original.
   int *resized_bkpts = (int *) malloc(N*MAXBREAKS * sizeof(int));
   int *resized_passages = (int *) malloc(N * sizeof(int));
   for (i = 0 ; i < N*MAXBREAKS ; i++) resized_bkpts[i] = 0;
   for (i = 0 ; i < N ; i++) resized_passages[i] = 0;

   for (l = 0, i = 0 ; i < N ; i++) {
      if (!remove[i]) {
         resized_passages[i] = passages[l];
         for (j = 0 ; j < MAXBREAKS ; j++) {
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
   seg->maxbreaks = MAXBREAKS;
   seg->nbreaks_opt = nbreaks_opt;
   seg->passages = resized_passages;
   seg->llikmat = resized_llikmat;
   seg->mllik = mllik;
   seg->bkpts = resized_bkpts;

   return;

}
