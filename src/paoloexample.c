

#include <stdio.h>
#include <stdlib.h>




#define GETTIME



#ifdef GETTIME
#include <sys/time.h>
struct timeval _t1,_t2;
double duration;

#define START_CLOCK   gettimeofday(&_t1,NULL ); 
#define END_CLOCK   gettimeofday(&_t2,NULL);   duration = (_t2.tv_sec-_t1.tv_sec)+ (double)(_t2.tv_usec-_t1.tv_usec)/1000000;    printf("----------> get time %e sec<------\n",duration); 
#endif /*  GETTIME */


#ifdef CLOCK
#include <time.h>
clock_t _t1,_t2;
double duration;

#define START_CLOCK   _t1 = clock(); 
#define END_CLOCK     _t2 = clock(); duration =  ((double)(_t2-_t1))/CLOCKS_PER_SEC; \
  printf("clock time %e s \n", duration); 
#endif

typedef int Mat ;
static int DEBUG=0;

#include <SparseBLAS.h>
#include <sorting.h>
#include <parsparsecoo.h>




int main() {
  int k, D,Ps;
  COO M,MT;
  COO temp;
  
  printf("K "); scanf("%d", &k);
  printf("D "); scanf("%d", &D);
  printf("P "); scanf("%d", &Ps);
  
  M = buildrandom_coo_list(k, D);

  printf("M %d %d %d \n", M.length, M.M, M.N); 

  MT = M;

  MT.data = (COOE*) malloc(M.length*sizeof(COOE));
  
  assert(MT.data);
    
  for (int i=0;i<M.length;i++) MT.data[i] = M.data[i];

  columnsort(&MT); // we really transpose the data so that the order is

  printf("MT %d %d %d \n", MT.length, MT.M, MT.N); 
  
  if (DEBUG) print_coo(M);

  
  START_CLOCK;
  if (Ps==0) 
    temp = matmul_coo(M,M,MT);
  else
    temp = matmul_coo_par(M,M,MT,Ps);
  END_CLOCK;


  if (DEBUG) print_coo(temp);
  
  free(MT.data);
  free(M.data);
  free(temp.data);

  return 0;
}

