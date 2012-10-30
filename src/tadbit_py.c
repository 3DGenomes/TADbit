/* @(#)tadbit_py.c
 */

#include "tadbit.c"

main(){
  double obs[1][10]={0.0,1.0,2.0,3.0,4.0,4.0,6.0,2.0,8.0,1.0};
  int n;
  const int m=1;
  int n_threads;
  const int verbose;
  const int speed;
  const int heuristic;
  /* output */
  tadbit_output *seg;
  int i;
  n=10;
  
  for (i = 0; i < 10; ++i) 
    printf(" %f", obs[0][i]);
  printf("\nhola\n");
  tadbit(obs, n, m, n_threads, verbose, speed, heuristic, seg);
}
