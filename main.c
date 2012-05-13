#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "tadbit.h"

int *read_int_matrix (const char *fname, int *n_ptr) {
/*
   * Read integers from a file and return a line-matrix (a pointer
   * of pointers to int). Assume that the first line is a header
   * containing as many space characters as columns, that the number
   * of rows is the same as the number of columns, and that each
   * row starts with the row name.
*/

   int i, j = 0, check_row = 0, check_col, n = 0;
   char line[65536], row_name[256], *a, *z;

   int *array;
   long u;

   // Garbage collection in case of early return.
   int *failure_exit() {
      free(array);
      return NULL;
   }

   FILE *f = fopen(fname, "r");
   if (f == NULL) {
      return NULL;
   }

   // Use first line to get 'n' by counting tabs.
   fgets(line, sizeof(line), f);
   for (i = 0 ; i < strlen(line)-1 ; i++) {
      if (isspace(line[i])) {
         n++;
      }
   }

   array = (int *) malloc (n*n * sizeof(int));

   while (fgets(line, sizeof(line), f)) {
      if (++check_row > n) {
         printf("check_row = %d\n", check_row);
         return failure_exit();
      }

      a = line;
      z = NULL;

      // Skip row name.
      if (sscanf(line, "%s", row_name)) {
         a += strlen(row_name);
      }

      check_col = 0;
      u = strtol(a, &z, 10);
      while (a != z) {
         check_col++;
         array[j++] = (int) u;
         a = z;
         u = strtol(a, &z, 10);
      }

      if (check_col != n) {
         printf("n = %d, check_col = %d\n", n, check_col);
         return failure_exit();
      }

   }

   if (check_row != n) {
      return failure_exit();
   }

   fclose(f);
   *n_ptr = n;
   return array;

}


int main (int argc, const char* argv[]) {
   int i, j, n, *read = read_int_matrix(argv[1], &n);
   if (read == NULL) {
      fprintf(stderr, "read error");
      exit(1);
   }

   double *counts = (double *) malloc(n*n * sizeof(double));
   for (i = 0 ; i < n*n ; i++) {
	   counts[i] = (double) read[i];
   }

   // Set NAs on rows and cols that have 0 diagonal term.
   for (i = 0 ; i < n ; i++) {
	   if (counts[i+i*n] == 0.0) {
		   for (j = 0 ; j < n ; j++) {
			   counts[j+i*n] = 0/0.0;
			   counts[i+j*n] = 0/0.0;
		   }
	   }
   }

   double **obs = &counts;
   int *val = (int *) malloc (n * sizeof(int));
   tadbit((const double **) obs, n, 1, 50, 0, 1, val);

   for (i = 0 ; i < n ; i++) {
      printf("%d ", val[i]);
   }
   printf("\n");
}
