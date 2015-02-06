#include "test.h"

char error_buffer[1024];
int backup;


const int ideal_matrix_20x20[400] = {
// Follows (50 * d^-0.5) on the two diagonal blocks, and
// (50 * d-.08) on the two off-diagonal blocks.
99,50,35,28,25,22,20,18,17,16,   7, 7, 6, 6, 6, 5, 5, 5, 4, 4,
50,99,50,35,28,25,22,20,18,17,   8, 7, 7, 6, 6, 6, 5, 5, 5, 4,
35,50,99,50,35,28,25,22,20,18,   9, 8, 7, 7, 6, 6, 6, 5, 5, 5,
28,35,50,99,50,35,28,25,22,20,  10, 9, 8, 7, 7, 6, 6, 6, 5, 5,
25,28,35,50,99,50,35,28,25,22,  11,10, 9, 8, 7, 7, 6, 6, 6, 5,
22,25,28,35,50,99,50,35,28,25,  13,11,10, 9, 8, 7, 7, 6, 6, 6,
20,22,25,28,35,50,99,50,35,28,  16,13,11,10, 9, 8, 7, 7, 6, 6,
18,20,22,25,28,35,50,99,50,35,  20,16,13,11,10, 9, 8, 7, 7, 6,
17,18,20,22,25,28,35,50,99,50,  28,20,16,13,11,10, 9, 8, 7, 7,
16,17,18,20,22,25,28,35,50,99,  50,28,20,16,13,11,10, 9, 8, 7,

 7, 8, 9,10,11,13,16,20,28,50,  99,50,35,28,25,22,20,18,17,16,
 7, 7, 8, 9,10,11,13,16,20,28,  50,99,50,35,28,25,22,20,18,17,
 6, 7, 7, 8, 9,10,11,13,16,20,  35,50,99,50,35,28,25,22,20,18,
 6, 6, 7, 7, 8, 9,10,11,13,16,  28,35,50,99,50,35,28,25,22,20,
 6, 6, 6, 7, 7, 8, 9,10,11,13,  25,28,35,50,99,50,35,28,25,22,
 5, 6, 6, 6, 7, 7, 8, 9,10,11,  22,25,28,35,50,99,50,35,28,25,
 5, 5, 6, 6, 6, 7, 7, 8, 9,10,  20,22,25,28,35,50,99,50,35,28,
 5, 5, 5, 6, 6, 6, 7, 7, 8, 9,  18,20,22,25,28,35,50,99,50,35,
 4, 5, 5, 5, 6, 6, 6, 7, 7, 8,  17,18,20,22,25,28,35,50,99,50,
 4, 4, 5, 5, 5, 6, 6, 6, 7, 7,  16,17,18,20,22,25,28,35,50,99,
};

const double expected_weights[400] = {
// Weights computed manually from the matrix above with R.
148225,162855,171710,177485,181720,184800,186725,188265,189805,193270,
193270,189805,188265,186725,184800,181720,177485,171710,162855,148225,
162855,178929,188658,195003,199656,203040,205155,206847,208539,212346,
212346,208539,206847,205155,203040,199656,195003,188658,178929,162855,
171710,188658,198916,205606,210512,214080,216310,218094,219878,223892,
223892,219878,218094,216310,214080,210512,205606,198916,188658,171710,
177485,195003,205606,212521,217592,221280,223585,225429,227273,231422,
231422,227273,225429,223585,221280,217592,212521,205606,195003,177485,
181720,199656,210512,217592,222784,226560,228920,230808,232696,236944,
236944,232696,230808,228920,226560,222784,217592,210512,199656,181720,
184800,203040,214080,221280,226560,230400,232800,234720,236640,240960,
240960,236640,234720,232800,230400,226560,221280,214080,203040,184800,
186725,205155,216310,223585,228920,232800,235225,237165,239105,243470,
243470,239105,237165,235225,232800,228920,223585,216310,205155,186725,
188265,206847,218094,225429,230808,234720,237165,239121,241077,245478,
245478,241077,239121,237165,234720,230808,225429,218094,206847,188265,
189805,208539,219878,227273,232696,236640,239105,241077,243049,247486,
247486,243049,241077,239105,236640,232696,227273,219878,208539,189805,
193270,212346,223892,231422,236944,240960,243470,245478,247486,252004,
252004,247486,245478,243470,240960,236944,231422,223892,212346,193270,
193270,212346,223892,231422,236944,240960,243470,245478,247486,252004,
252004,247486,245478,243470,240960,236944,231422,223892,212346,193270,
189805,208539,219878,227273,232696,236640,239105,241077,243049,247486,
247486,243049,241077,239105,236640,232696,227273,219878,208539,189805,
188265,206847,218094,225429,230808,234720,237165,239121,241077,245478,
245478,241077,239121,237165,234720,230808,225429,218094,206847,188265,
186725,205155,216310,223585,228920,232800,235225,237165,239105,243470,
243470,239105,237165,235225,232800,228920,223585,216310,205155,186725,
184800,203040,214080,221280,226560,230400,232800,234720,236640,240960,
240960,236640,234720,232800,230400,226560,221280,214080,203040,184800,
181720,199656,210512,217592,222784,226560,228920,230808,232696,236944,
236944,232696,230808,228920,226560,222784,217592,210512,199656,181720,
177485,195003,205606,212521,217592,221280,223585,225429,227273,231422,
231422,227273,225429,223585,221280,217592,212521,205606,195003,177485,
171710,188658,198916,205606,210512,214080,216310,218094,219878,223892,
223892,219878,218094,216310,214080,210512,205606,198916,188658,171710,
162855,178929,188658,195003,199656,203040,205155,206847,208539,212346,
212346,208539,206847,205155,203040,199656,195003,188658,178929,162855,
148225,162855,171710,177485,181720,184800,186725,188265,189805,193270,
193270,189805,188265,186725,184800,181720,177485,171710,162855,148225,
};


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

   int *obs[2];
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

   return;

}

void
test_tadbit
(void)
{

   // -- INPUT -- //
   int *obs[2];
   int obs0[400];
   int obs1[400];
   memcpy(obs0, ideal_matrix_20x20, 400 * sizeof(int));
   memcpy(obs1, ideal_matrix_20x20, 400 * sizeof(int));
   obs[0] = obs0;
   obs[1] = obs1;

   // -- OUTPUT -- //
   tadbit_output *seg = malloc(sizeof(tadbit_output));
   char *remove = (char *) malloc (400 * sizeof(char));
   for (int j = 0 ; j < 400 ; j++){
    remove[j] = 0; // automatic casting into char
  }

   tadbit(obs, remove, 20, 2, 1, 0, 20, 0, 1, seg);

   // Check max breaks and optimal number of breaks.
   g_assert_cmpint(seg->maxbreaks, ==, 4);
   g_assert_cmpint(seg->nbreaks_opt, ==, 1);
   // Check the position of the optimal break.
   for (int i = 0 ; i < 20 ; i++) {
      g_assert_cmpint(seg->bkpts[i+1*20], == , i == 9);
   }

   // Check the computed weights.
   /* for (int j = 0 ; j < 20 ; j++) { */
   /*   for (int i = j ; i < 20 ; i++) { */
   /*    g_assert_cmpfloat( */
   /*       seg->weights[0][i+j*20], ==,  expected_weights[i+j*20] */
   /*    ); */
   /* } */
   /* } */

   destroy_tadbit_output(seg);
   
}

void
test_ll
(void)
{

   double lg[400] = {0};
   double *c = malloc(21 * sizeof(double));

   double w[400] = {[0 ... 399] = 1.0};
   double d[400];
   int C[400];

   _cache_index = C;
   _max_cache_index = 21;

   for (int j = 0 ; j < 20 ; j++) {
      for (int i = 0 ; i < 20 ; i++) {
         d[i+j*20] = log(abs(j-i));
         C[i+j*20] = abs(j-i);
      }
   }

   double loglik1 = ll(20, 0, 9, 0, 9, 1, ideal_matrix_20x20, d, w, lg, c);
   // Value checked manually with R. The value is sensitive to
   // the value of the estimates, which is why the  precision
   // cannot be higher than 0.1.
   g_assert_cmpfloat(abs(loglik1-6138.2), <, 1e-1);

   // Check symmetry/reproducibility.
   double loglik2 = ll(20, 10, 19, 10, 19, 1, ideal_matrix_20x20, d, w, lg, c);
   g_assert_cmpfloat(abs(loglik1-loglik2), <, 1e-12);

   // Same as above, checked manually with R.
   loglik1 = ll(20, 0, 9, 10, 19, 0, ideal_matrix_20x20, d, w, lg, c);
   g_assert_cmpfloat(abs(loglik1-3036.8), <, 1e-1);

   // Check symmetry/reproducibility again.
   loglik2 = ll(20, 10, 19, 0, 9, 0, ideal_matrix_20x20, d, w, lg, c);
   g_assert_cmpfloat(abs(loglik1-loglik2), <, 1e-12);

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
   char *remove = (char *) malloc (400 * sizeof(char));
   tadbit(obs, remove, 3191, 2, 8, 1, 200, 0, 0, seg);
   unredirect_sderr();

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
   g_test_add_func("/ll", test_ll);
   g_test_add_func("/enforce_symmetry", test_enforce_symmetry);
   g_test_add_func("/tadbit", test_tadbit);
   if (g_test_thorough()) {
      g_test_add_func("/tadbit_on_real_input", test_tadbit_on_real_input);
   }

   int g_test_result = g_test_run();
   close(backup);

   return g_test_result;

}
