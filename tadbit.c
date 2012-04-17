#include <math.h>
#include <stdlib.h>
#include <R.h>

#define TOLERANCE 1e-6


double ml_ab(double *k, double *d, double *ab, int n) {
/*
 * The 2-double-array 'ab' is upated in place and the log-likelihood
 * is returned.
 * The fitted model (by maximum likelihood) is Poisson with lambda
 * paramter such that lambda = exp(a + b*d). So the full log-likelihood
 * of the model is Sigma -exp(a + b*d_i) + k_i(a + b*d_i).
 */

   int i, j;
   double a = ab[0], b = ab[1], llik;
   double f = 0.0, g = 0.0, dfda, dfdb, dgda, dgdb;
   double denom, tmp; // 'tmp' is used as computation intermediate.

   // Initialize 'f' and 'g'.
   for (i = 0 ; i < n ; i++) {
      tmp  =  exp(a+b*d[i])-k[i];
      f   +=  tmp;
      g   +=  tmp * d[i];
   }

   // Newton-Raphson until gradient function < TOLERANCE.
   while ((f*f + g*g) > TOLERANCE) {

      // Compute the derivatives.
      dfda = dfdb = dgda = dgdb = 0.0;
      for (i = 0 ; i < n ; i++) {
         tmp   =   exp(a+b*d[i]);
         dfda +=   tmp;
         dgda +=   tmp * d[i];
         dgdb +=   tmp * d[i]*d[i];
      }
      dfdb = dgda;

      // NB: skip gradient test.
      denom = dfdb*dgda - dfda*dgdb;
      a += (f*dgdb - g*dfdb) / denom;
      b += (g*dfda - f*dgda) / denom;

      // Update 'f' and 'g' for testing.
      f = g = 0.0;
      for (i = 0 ; i < n ; i++) {
         tmp  =  exp(a+b*d[i])-k[i];
         f   +=  tmp;
         g   +=  tmp * d[i];
      }
   }

   // Update 'ab' in place.
   ab[0] = a; ab[1] = b;

   // Compute log-likelihood (starting from 'dfda').
   llik = -dfda;
   for (i = 0 ; i < n ; i++) {
      llik += k[i] * (a + b*d[i]);
   }

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

void fill_ml_matrix(double *obs, double *llik, int *dim) {
/*
 * Fill in the upper half of 'llik' with log-likelihood values of
 * segments starting at index i and ending at index j.
 */

   int n = (int) *dim;
   int i, j, k;

   double ab[3][2] = {{0.0,0.0},{0.0,0.0}, {0.0,0.0}};


   // Build the (n,n) distance matrix.
   double *dis = (double *) malloc(n*n * sizeof(double));

   k = 0;
   for (i = 0; i < n ; i++) {
      for (j = 0; j < n ; j++) {
         //dis[k++] = log(1+abs(i-j));
         dis[k++] = abs(i-j);
      }
   }

   // Allocate max possible size to blocks matrices.
   double **d_blk = (double **) malloc(3 * sizeof(double *));
   d_blk[0] = (double *) malloc(n*n/4 * sizeof(double));
   d_blk[1] = (double *) malloc(n*(n-1)/2 * sizeof(double));
   d_blk[2] = (double *) malloc(n*n/4 * sizeof(double));

   double **k_blk = (double **) malloc(3 * sizeof(double *));
   k_blk[0] = (double *) malloc(n*n/4 * sizeof(double));
   k_blk[1] = (double *) malloc(n*(n-1)/2 * sizeof(double));
   k_blk[2] = (double *) malloc(n*n/4 * sizeof(double));


   for (i = 0 ; i < n-2 ; i++) {
      for (j = i+2 ; j < n ; j++) {
         // Cut the (i,j)-blocks.
         k_blk = break_in_blocks(obs, n, i, j, k_blk);
         d_blk = break_in_blocks(dis, n, i, j, d_blk);

         // Get the likelihood per block and sum.
         llik[i+j*n] =
               ml_ab(k_blk[0], d_blk[0], ab[0], i*(j-i+1))       / 2  +
               ml_ab(k_blk[1], d_blk[1], ab[1], (j-i)*(j-i+1)/2)      +
               ml_ab(k_blk[2], d_blk[2], ab[2], (n-j-1)*(j-i+1)) / 2;
      }
   }

}

// TODO: write the maximization function.
