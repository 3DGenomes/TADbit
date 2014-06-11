#include "visibility.h"

double **visibility(int **obs, int m, int n){
  // Compute row/column sums (identical by symmetry).
  double **rowsums = (double **) malloc(m * sizeof(double *));
  int k;
  int i;
  int j;
  int l;
  for (k = 0 ; k < m ; k++) {
    rowsums[k] = (double *) malloc(n * sizeof(double));
    for (i = 0 ; i < n ; i++) 
      rowsums[k][i] = 0.0;
  }

   for (k = 0 ; k < m ; k++)
     for (i = 0 ; i < n ; i++)
       for (j = 0 ; j < n ; j++)
	 rowsums[k][i] += obs[k][i+j*n];

   // compute the weights.
   double **weights = (double **) malloc(m * sizeof(double *));
   for (k = 0 ; k < m ; k++) {
      weights[k] = (double *) malloc(n*n * sizeof(double));
      for (i = 0 ; i < n*n ; i++) weights[k][i] = 0.0;
      // Compute product.
      for (j = 0 ; j < n ; j++) {
	for (i = 0 ; i < n ; i++) {
	  weights[k][i+j*n] = rowsums[k][i]*rowsums[k][j];
	}
      }
   }

   // We don't need the row/column sums any more.
   for (l = 0 ; l < m ; l++) free(rowsums[l]);
   free(rowsums);
   return weights;
}
