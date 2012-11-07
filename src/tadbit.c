#include <stdlib.h>
#include <stdio.h>
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
   for (index = 0 ; index < n ; index++) c[index] = NAN;

   for (j = j_low ; j < j_high ; j++) {
      i_high = diag ? j : _i+1;
      for (i = i_low ; i < i_high ; i++) {
         // Retrieve value of the exponential from cache.
         index = _cache_index[i+j*n];
         if (isnan(c[index])) {
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
//   'n': row/column number of the counts.                              
//   'i_': first value of index i (row).                                
//   '_i': last value of index i (row).                                 
//   'j_': first value of index j (column).                             
//   '_j': last value of index j (column).                              
//   'diag': whether the block is half-diagonal (middle block).         
//   'k': raw hiC counts as double.                                     
//   'd': distances from diagonal in log.                               
//   'w': weights measuring hiC bias.                                   
//   'lg': log-gamma terms.                                             
//   'c': address of an array of double for caching.                    
//                                                                      
// RETURN:                                                              
//   The maximum log-likelihood of a block of hiC data.                 
//                                                                      

   if ((i_ >= _i) || (j_ >= _j)) return 0.0;

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

      for (index = 0 ; index < n ; index++) c[index] = NAN;
      // Compute the derivatives.
      dfda = dfdb = dgda = dgdb = 0.0;

      for (j = j_low ; j < j_high ; j++) {
         i_high = diag ? j : _i+1;
         for (i = i_low ; i < i_high ; i++) {
            index = _cache_index[i+j*n];
            // Retrive value of the exponential from cache.
            if (isnan(c[index])) {
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

   //////////////////////////////////////////////////////////////////////
   // The last call to 'calcfg' has set the cache to the right values. //
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

void
DPwalk(
  /* input */
  const double *llik_mat,
  const int n,
  const int MAXBREAKS,
  /* output */
  double *mllik,
  int *breakpoints
){
// SYNOPSIS:                                                            
//   Dynamic programming algorithm to compute the most likely position  
//   of breakpoints given a matrix of slice maximum log-likelihood.     
//                                                                      
// PARAMETERS:                                                          
//   '*llik_mat': matrix of maximum log-likelihood values.              
//   'n': row/col number of 'llik_mat'.                                 
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

   // Initialize 'old_llik' to the first line of 'llik_mat' containing
   // the log-likelihood of segments starting at index 0.
   for (i = 0 ; i < n ; i++) {
      old_llik[i] = llik_mat[i*n];
      new_llik[i] = -INFINITY;
   }


   // Dynamic programming.
   for (nbreaks = 1 ; nbreaks < MAXBREAKS ; nbreaks++) {
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
//   'arg': thread arguments (see header file for definition).          
//                                                                      
// RETURN:                                                              
//   'void'                                                             
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
   for (i = 0 ; i < n ; i++) c[i] = 0.0;

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
allocate_task(
  int *skip,
  const int i0,
  const int j0,
  const int n
){
// SYNOPSIS:                                                            
//   Create or update thread jobs. For an approximate TAD defined by    
//   ('i0', 'j0'), create jobs that will compute the log-likelihood of  
//   this TAD, and all single splits of this TAD.                       
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

   for (j = i0 ; j < j0 ; j++)
      skip[i0+j*n] = 0;
   for (i = i0 ; i < j0 ; i++)
      skip[i+j0*n] = 0;
}

void
tadbit(
  /* input */
  double **obs,
  int n,
  const int m,
  int n_threads,
  const int verbose,
  int max_tad_size,
  const int do_not_use_heuristic,
  /* output */
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

   const int MAXBREAKS = n/5;

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

   int *skip = (int *) malloc(n*n *sizeof(int));

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
         fprintf(stderr, "applying heuristic pruning\n");
      }
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
      DPwalk(heur_score, n, MAXBREAKS, mllik, bkpts);

      free(heur_score);
      free(S);

      int i0;
      int j0;
      // Create a thread job for each approximate TAD. (see
      // 'allocate_task').
      for (i = 0 ; i < n*n ; i++) skip[i] = 1;
      for (j = 1 ; j < MAXBREAKS ; j++) {
         i0 = 0;
         for (i = 0 ; i < n ; i++) {
            if (bkpts[i+j*n]) {
               j0 = i;
               allocate_task(skip, i0, j0, n);
               i0 = i+1;
            }
         }
      }

      // Allocate estimation of the log likelihood for all small
      // TADs (less than 3 bins).
      for (j = 6 ; j < n ; j++)
      for (i = j-6 ; i < j-3 ; i++)
         skip[i+j*n] = 0;

      // Allocate jobs at the ends of the chromosomes/units because
      // these regions are a bit noisier.
      for (j = 1 ; j < 101 ; j++)
      for (i = 0 ; i < j-3 ; i++)
         skip[i+j*n] = 0;
      for (j = n-101 ; j < n ; j++)
      for (i = n-101 ; i < j-3 ; i++)
         skip[i+j*n] = 0;

   } // End of heuristic task allocation.


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

   // Initialize task queue.
   n_processed = 0;
   n_to_process = 0;
   for (i = 0 ; i < n*n ; i++) n_to_process += (1-skip[i]);
   taskQ_i = 0;

   // Instantiate threads and start running jobs.
   for (i = 0 ; i < n_threads ; i++) tid[i] = 0;
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
   free(_cache_index);

   // The matrix 'llikmat' now contains the log-likelihood of the
   // segments. The breakpoints are found by dynamic programming.
   DPwalk(llikmat, n, MAXBREAKS, mllik, bkpts);

   // Get optimal number of breaks by AIC.
   double AIC = -INFINITY;
   int n_params;
   int nbreaks_opt;
   for (nbreaks_opt = 1 ; nbreaks_opt  < MAXBREAKS ; nbreaks_opt++) {
      n_params = nbreaks_opt + m*(8 + nbreaks_opt*6);
      if (AIC > mllik[nbreaks_opt] - n_params) {
         break;
      }
      else {
         AIC = mllik[nbreaks_opt] - n_params;
      }
   }
   nbreaks_opt -= 1;

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
      DPwalk(llikmatcpy, n, MAXBREAKS, mllikcpy, bkptscpy);
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
