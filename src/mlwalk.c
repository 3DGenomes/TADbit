#include <stdlib.h>

int **
mlwalk(
  /* input */
  const double *llik,
  const int n,
  const int m,
  /* output */
  double *mllk
){

   // The max number of segments is 1/4 row/col number.
   const int max_n_breaks = n/4;

   int i;
   int j;
   int nbreaks;

   double new_llik[n];
   double old_llik[n];

   // Return value, the optimal breakpoints for very number
   // of breakpoint.
   int **breakpoints = (int **) malloc(max_n_breaks * sizeof(int *));

   // Breakpoint lists. The first index (line) is the end of the slice,
   // the second is 1 if this position is an end (breakpoint).
   // These must be allocated from the heap because 'n' can be large.
   int *new_bkpt_list = (int *) malloc(n*max_n_breaks * sizeof(int));
   int *old_bkpt_list = (int *) malloc(n*max_n_breaks * sizeof(int));

   // Initializations.
   for (i = 0 ; i < n*n ; i++) {
      new_bkpt_list[i] = old_bkpt_list[i] = 0;
   }

   for (i = 1 ; i < max_n_breaks ; i++) {
      mllik[i] = NAN;
      breakpoints[i] = (int *) malloc(i * sizeof(int));
      for (j = 0 ; j < i ; j++) {
         breakpoints[i][j] = -1;
      }
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
      for (j = 0, i = 0 ; i < n ; i++) {
         old_llik[i] = new_llik[i];
         if (new_bkpt_list[n-1+i*n]) {
            breakpoints[nbreaks][j] = i;
            j++;
         }
      }

   }

   free(new_bkpt_list);
   free(old_bkpt_list);

   return breakpoints;

}
