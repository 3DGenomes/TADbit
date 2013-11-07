#include "test.h"

char error_buffer[1024];
int backup;


void
redirect_stderr_to
(char buffer[])
{
   // Flush stderr, redirect to /dev/null and set buffer.
   fflush(stderr);
   int temp = open("/dev/null", O_WRONLY);
   dup2(temp, STDERR_FILENO);
   memset(buffer, '\0', 1024 * sizeof(char));
   setvbuf(stderr, buffer, _IOFBF, 1024);
   close(temp);
   // Fill the buffer (needed for reset).
   fprintf(stderr, "fill the buffer");
   fflush(stderr);
}

void
unredirect_sderr
(void)
{ 
   fflush(stderr);
   dup2(backup, STDERR_FILENO);
   setvbuf(stderr, NULL, _IONBF, 0);
}

void
test_enforce_symmetry
(void)
{

   int case1[16] = {
      1,2,3,4,
      1,2,3,4,
      1,2,3,4,
      1,2,3,4,
   };
   int case2[16] = {
      1,1,1,1,
      2,2,2,2,
      3,3,3,3,
      4,4,4,4,
   };

   int **obs = malloc(2 * sizeof(int *));
   obs[0] = case1;
   obs[1] = case2;

   int expected[16] = {
      1,3,4,5,
      3,2,5,6,
      4,5,3,7,
      5,6,7,4,
   };
  
   redirect_stderr_to(error_buffer);
   int asymmetry = enforce_symmetry(obs, 4, 2);
   unredirect_sderr();   

   // Check the output.
   g_assert_cmpint(asymmetry, ==, 1);
   for (int i = 0 ; i < 16 ; i++) {
      g_assert_cmpint(expected[i], ==, obs[0][i]);
      g_assert_cmpint(expected[i], ==, obs[1][i]);
   }

   // Check the error message.
   g_assert_cmpstr(error_buffer, ==,
         "input matrix not symmetric: symmetrizing\n");


   // Now the matrices are symmetric.
   asymmetry = enforce_symmetry(obs, 4, 2);
   g_assert_cmpint(asymmetry, ==, 0);

   free(obs);
   return;

}

void
test_tadbit
(void)
{

   // -- INPUT -- //
   int obs1[400] = {
    99, 50, 35, 28, 25, 22, 20, 18, 17, 16,   7,  7,  6,  6,  6,  5,  5,  5,  4,  4,
    50, 99, 50, 35, 28, 25, 22, 20, 18, 17,   8,  7,  7,  6,  6,  6,  5,  5,  5,  4,
    35, 50, 99, 50, 35, 28, 25, 22, 20, 18,   9,  8,  7,  7,  6,  6,  6,  5,  5,  5,
    28, 35, 50, 99, 50, 35, 28, 25, 22, 20,  10,  9,  8,  7,  7,  6,  6,  6,  5,  5,
    25, 28, 35, 50, 99, 50, 35, 28, 25, 22,  11, 10,  9,  8,  7,  7,  6,  6,  6,  5,
    22, 25, 28, 35, 50, 99, 50, 35, 28, 25,  13, 11, 10,  9,  8,  7,  7,  6,  6,  6,
    20, 22, 25, 28, 35, 50, 99, 50, 35, 28,  16, 13, 11, 10,  9,  8,  7,  7,  6,  6,
    18, 20, 22, 25, 28, 35, 50, 99, 50, 35,  20, 16, 13, 11, 10,  9,  8,  7,  7,  6,
    17, 18, 20, 22, 25, 28, 35, 50, 99, 50,  28, 20, 16, 13, 11, 10,  9,  8,  7,  7,
    16, 17, 18, 20, 22, 25, 28, 35, 50, 99,  50, 28, 20, 16, 13, 11, 10,  9,  8,  7,

     7,  8,  9, 10, 11, 13, 16, 20, 28, 50,   99, 50, 35, 28, 25, 22, 20, 18, 17, 16,
     7,  7,  8,  9, 10, 11, 13, 16, 20, 28,   50, 99, 50, 35, 28, 25, 22, 20, 18, 17,
     6,  7,  7,  8,  9, 10, 11, 13, 16, 20,   35, 50, 99, 50, 35, 28, 25, 22, 20, 18,
     6,  6,  7,  7,  8,  9, 10, 11, 13, 16,   28, 35, 50, 99, 50, 35, 28, 25, 22, 20,
     6,  6,  6,  7,  7,  8,  9, 10, 11, 13,   25, 28, 35, 50, 99, 50, 35, 28, 25, 22,
     5,  6,  6,  6,  7,  7,  8,  9, 10, 11,   22, 25, 28, 35, 50, 99, 50, 35, 28, 25,
     5,  5,  6,  6,  6,  7,  7,  8,  9, 10,   20, 22, 25, 28, 35, 50, 99, 50, 35, 28,
     5,  5,  5,  6,  6,  6,  7,  7,  8,  9,   18, 20, 22, 25, 28, 35, 50, 99, 50, 35,
     4,  5,  5,  5,  6,  6,  6,  7,  7,  8,   17, 18, 20, 22, 25, 28, 35, 50, 99, 50,
     4,  4,  5,  5,  5,  6,  6,  6,  7,  7,   16, 17, 18, 20, 22, 25, 28, 35, 50, 99,
   };


   int **obs = malloc(2 * sizeof(int *));
   obs[0] = obs1;
   obs[1] = obs1;

   // -- OUTPUT -- //
   tadbit_output *seg = malloc(sizeof(tadbit_output));

   double expected_weights[400] = {
      16, 17.6, 18.5, 19.1, 19.6, 19.9, 20.1, 20.3, 20.5, 20.8, 20.8, 20.5, 20.3, 20.1, 19.9, 19.6, 19.1, 18.5, 17.6, 16,
      17.6, 19.3, 20.3, 21, 21.5, 21.9, 22.1, 22.3, 22.5, 22.9, 22.9, 22.5, 22.3, 22.1, 21.9, 21.5, 21, 20.3, 19.3, 17.6,
      18.5, 20.3, 21.5, 22.2, 22.7, 23.1, 23.3, 23.5, 23.7, 24.1, 24.1, 23.7, 23.5, 23.3, 23.1, 22.7, 22.2, 21.5, 20.3, 18.5,
      19.1, 21, 22.2, 22.9, 23.5, 23.9, 24.1, 24.3, 24.5, 25, 25, 24.5, 24.3, 24.1, 23.9, 23.5, 22.9, 22.2, 21, 19.1,
      19.6, 21.5, 22.7, 23.5, 24, 24.4, 24.7, 24.9, 25.1, 25.6, 25.6, 25.1, 24.9, 24.7, 24.4, 24, 23.5, 22.7, 21.5, 19.6,
      19.9, 21.9, 23.1, 23.9, 24.4, 24.8, 25.1, 25.3, 25.5, 26, 26, 25.5, 25.3, 25.1, 24.8, 24.4, 23.9, 23.1, 21.9, 19.9,
      20.1, 22.1, 23.3, 24.1, 24.7, 25.1, 25.4, 25.6, 25.8, 26.3, 26.3, 25.8, 25.6, 25.4, 25.1, 24.7, 24.1, 23.3, 22.1, 20.1,
      20.3, 22.3, 23.5, 24.3, 24.9, 25.3, 25.6, 25.8, 26, 26.5, 26.5, 26, 25.8, 25.6, 25.3, 24.9, 24.3, 23.5, 22.3, 20.3,
      20.5, 22.5, 23.7, 24.5, 25.1, 25.5, 25.8, 26, 26.2, 26.7, 26.7, 26.2, 26, 25.8, 25.5, 25.1, 24.5, 23.7, 22.5, 20.5,
      20.8, 22.9, 24.1, 25, 25.6, 26, 26.3, 26.5, 26.7, 27.2, 27.2, 26.7, 26.5, 26.3, 26, 25.6, 25, 24.1, 22.9, 20.8,
      20.8, 22.9, 24.1, 25, 25.6, 26, 26.3, 26.5, 26.7, 27.2, 27.2, 26.7, 26.5, 26.3, 26, 25.6, 25, 24.1, 22.9, 20.8,
      20.5, 22.5, 23.7, 24.5, 25.1, 25.5, 25.8, 26, 26.2, 26.7, 26.7, 26.2, 26, 25.8, 25.5, 25.1, 24.5, 23.7, 22.5, 20.5,
      20.3, 22.3, 23.5, 24.3, 24.9, 25.3, 25.6, 25.8, 26, 26.5, 26.5, 26, 25.8, 25.6, 25.3, 24.9, 24.3, 23.5, 22.3, 20.3,
      20.1, 22.1, 23.3, 24.1, 24.7, 25.1, 25.4, 25.6, 25.8, 26.3, 26.3, 25.8, 25.6, 25.4, 25.1, 24.7, 24.1, 23.3, 22.1, 20.1,
      19.9, 21.9, 23.1, 23.9, 24.4, 24.8, 25.1, 25.3, 25.5, 26, 26, 25.5, 25.3, 25.1, 24.8, 24.4, 23.9, 23.1, 21.9, 19.9,
      19.6, 21.5, 22.7, 23.5, 24, 24.4, 24.7, 24.9, 25.1, 25.6, 25.6, 25.1, 24.9, 24.7, 24.4, 24, 23.5, 22.7, 21.5, 19.6,
      19.1, 21, 22.2, 22.9, 23.5, 23.9, 24.1, 24.3, 24.5, 25, 25, 24.5, 24.3, 24.1, 23.9, 23.5, 22.9, 22.2, 21, 19.1,
      18.5, 20.3, 21.5, 22.2, 22.7, 23.1, 23.3, 23.5, 23.7, 24.1, 24.1, 23.7, 23.5, 23.3, 23.1, 22.7, 22.2, 21.5, 20.3, 18.5,
      17.6, 19.3, 20.3, 21, 21.5, 21.9, 22.1, 22.3, 22.5, 22.9, 22.9, 22.5, 22.3, 22.1, 21.9, 21.5, 21, 20.3, 19.3, 17.6,
      16, 17.6, 18.5, 19.1, 19.6, 19.9, 20.1, 20.3, 20.5, 20.8, 20.8, 20.5, 20.3, 20.1, 19.9, 19.6, 19.1, 18.5, 17.6, 16,
   };

   tadbit(obs, 20, 2, 1, 0, 20, 1, 1, seg);

   // Check max breaks and optimal number of breaks.
   g_assert_cmpint(seg->maxbreaks, ==, 4);
   g_assert_cmpint(seg->nbreaks_opt, ==, 1);
   // Check the position of the optimal break.
   for (int i = 0 ; i < 20 ; i++) {
      g_assert_cmpint(seg->bkpts[i+1*20], == , i == 9);
   }


   // Check the computed weights.
   for (int j = 0 ; j < 20 ; j++) {
      for (int i = j ; i < 20 ; i++) {
         g_assert_cmpfloat(
            abs(seg->weights[0][i+j*20] / 9272 -
               expected_weights[i+j*20]), <, .1
         );
      }
   }

   free(obs);
   destroy_tadbit_output(seg);
   
}

void
test_ll
(void)
{

   // -- INPUT -- //
   int obs[400] = {
    99, 50, 35, 28, 25, 22, 20, 18, 17, 16,   7,  7,  6,  6,  6,  5,  5,  5,  4,  4,
    50, 99, 50, 35, 28, 25, 22, 20, 18, 17,   8,  7,  7,  6,  6,  6,  5,  5,  5,  4,
    35, 50, 99, 50, 35, 28, 25, 22, 20, 18,   9,  8,  7,  7,  6,  6,  6,  5,  5,  5,
    28, 35, 50, 99, 50, 35, 28, 25, 22, 20,  10,  9,  8,  7,  7,  6,  6,  6,  5,  5,
    25, 28, 35, 50, 99, 50, 35, 28, 25, 22,  11, 10,  9,  8,  7,  7,  6,  6,  6,  5,
    22, 25, 28, 35, 50, 99, 50, 35, 28, 25,  13, 11, 10,  9,  8,  7,  7,  6,  6,  6,
    20, 22, 25, 28, 35, 50, 99, 50, 35, 28,  16, 13, 11, 10,  9,  8,  7,  7,  6,  6,
    18, 20, 22, 25, 28, 35, 50, 99, 50, 35,  20, 16, 13, 11, 10,  9,  8,  7,  7,  6,
    17, 18, 20, 22, 25, 28, 35, 50, 99, 50,  28, 20, 16, 13, 11, 10,  9,  8,  7,  7,
    16, 17, 18, 20, 22, 25, 28, 35, 50, 99,  50, 28, 20, 16, 13, 11, 10,  9,  8,  7,

     7,  8,  9, 10, 11, 13, 16, 20, 28, 50,   99, 50, 35, 28, 25, 22, 20, 18, 17, 16,
     7,  7,  8,  9, 10, 11, 13, 16, 20, 28,   50, 99, 50, 35, 28, 25, 22, 20, 18, 17,
     6,  7,  7,  8,  9, 10, 11, 13, 16, 20,   35, 50, 99, 50, 35, 28, 25, 22, 20, 18,
     6,  6,  7,  7,  8,  9, 10, 11, 13, 16,   28, 35, 50, 99, 50, 35, 28, 25, 22, 20,
     6,  6,  6,  7,  7,  8,  9, 10, 11, 13,   25, 28, 35, 50, 99, 50, 35, 28, 25, 22,
     5,  6,  6,  6,  7,  7,  8,  9, 10, 11,   22, 25, 28, 35, 50, 99, 50, 35, 28, 25,
     5,  5,  6,  6,  6,  7,  7,  8,  9, 10,   20, 22, 25, 28, 35, 50, 99, 50, 35, 28,
     5,  5,  5,  6,  6,  6,  7,  7,  8,  9,   18, 20, 22, 25, 28, 35, 50, 99, 50, 35,
     4,  5,  5,  5,  6,  6,  6,  7,  7,  8,   17, 18, 20, 22, 25, 28, 35, 50, 99, 50,
     4,  4,  5,  5,  5,  6,  6,  6,  7,  7,   16, 17, 18, 20, 22, 25, 28, 35, 50, 99,
   };

   double lg[400] = {0};
   double *c = malloc(21 * sizeof(double));

   double w[400] = {[0 ... 399] = 1.0};
   double d[400];
   int C[400];

   double rowsums[20];
   double total = 0;
   for (int i = 0 ; i < 400 ; i++) total += obs[i];

   _cache_index = C;
   _max_cache_index = 21;

   for (int j = 0 ; j < 20 ; j++) {
      rowsums[j] = 0.0;
      for (int i = 0 ; i < 20 ; i++) {
         rowsums[j] += obs[i+j*20];
         d[i+j*20] = log(abs(j-i));
         C[i+j*20] = abs(j-i);
      }
   }

   double loglik = ll(20, 0, 9, 0, 9, 1, obs, d, w, lg, c);

   free(c);

}


void
test_tadbit_on_real_input
(void)
{
   // -- Open input file (the code is hopelessly not portable) -- //
   char *pch; 
   int i;

   int *obs[2];
   obs[0] = malloc(3191*3191 * sizeof(int));
   obs[1] = malloc(3191*3191 * sizeof(int));

   FILE *f[2];
   f[0] = fopen("data/hESC_chr19-rep1.txt", "r");
   f[1] = fopen("data/hESC_chr19-rep2.txt", "r");

   g_assert(f[0] != NULL);
   g_assert(f[1] != NULL);

   // `getline` is available because we define _GNU_SOURCE.
   char *line = NULL;
   size_t len = 0;
   ssize_t read;

   // Read both files the same way.
   for (int j = 0 ; j < 2 ; j++) {
      // Discard header.
      //read = getline(&line, &len, f[j]);
      i = 0;
      while ((read = getline(&line, &len, f[j])) != -1) {
         pch = strtok(line, "\t");
         // Discard row name.
         //pch = strtok(NULL, "\t");
         while (pch != NULL) {
            obs[j][i++] = atoi(pch);
            pch = strtok(NULL, "\t");
         }
         // Check the number of items in each row.
         g_assert_cmpint(i % 3191, == , 0);
      }
      // Check the total number of items read.
      g_assert_cmpint(i, == , 3191*3191);
   }

   fclose(f[0]);
   fclose(f[1]);
   free(line);

   tadbit_output *seg = malloc(sizeof(tadbit_output));
   redirect_stderr_to(error_buffer);
   tadbit(obs, 3191, 2, 8, 1, 200, 0, 0, seg);
   unredirect_sderr();

   int nTADs1 = seg->nbreaks_opt;

   //tadbit(obs, 1413, 2, 16, 1, 200, 0, 0, seg);
   //int nTADs2 = seg->nbreaks_opt;

   //fprintf(stderr, "%d, %d\n", nTADs1, nTADs2);
   destroy_tadbit_output(seg);

   free(obs[0]);
   free(obs[1]);

}


int
main(
   int argc,
   char **argv
)
{

   // Store 'stderr', file descriptor.
   backup = dup(STDERR_FILENO);

   g_test_init(&argc, &argv, NULL);
   //g_test_add_func("/ll", test_ll);
   g_test_add_func("/enforce_symmetry", test_enforce_symmetry);
   g_test_add_func("/tadbit", test_tadbit);
   g_test_add_func("/tadbit_on_real_input", test_tadbit_on_real_input);

   int g_test_result = g_test_run();
   close(backup);

   return g_test_result;

}
