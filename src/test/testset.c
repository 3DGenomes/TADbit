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
test_tadbit
(void)
{

   // -- INPUT -- //
   int obs1[36] = {
      32, 16,  8,  4,  2,  1,
      16, 32, 16,  8,  4,  2,
       8, 16, 32, 16,  8,  4,
       4,  8, 16, 32, 16,  8,
       2,  4,  8, 16, 32, 16,
       1,  2,  4,  8, 16, 32,
   };

   int obs2[36] = {
      32, 16,  8,  4,  2,  1,
      16, 32, 16,  8,  4,  2,
       8, 16, 32, 16,  8,  4,
       4,  8, 16, 32, 16,  8,
       2,  4,  8, 16, 32, 16,
       1,  2,  4,  8, 16, 32,
   };


   int **obs = malloc(2 * sizeof(int *));
   obs[0] = obs1;
   obs[1] = obs2;

   // -- OUTPUT -- //
   tadbit_output *seg = malloc(sizeof(tadbit_output));

   tadbit(obs, 6, 2, 1, 0, 6, 0, 0, seg);

   g_assert_cmpint(seg->maxbreaks, ==, 1);

   free(obs);
   destroy_tadbit_output(seg);
   
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
   g_test_add_func("/tadbit", test_tadbit);

   int g_test_result = g_test_run();
   close(backup);

   return g_test_result;

}
