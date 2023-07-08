

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

//#define GRAPH_PATH 1
static int DEBUG=0;

#include <SparseBLAS.h>
#include <sorting.h>
#include <parsparsecoo.h>
#include <parsorting.h>



int main(int argc, char **argv) {
  int k, D,Ps;
  COOMB M,MT;
  COOMB temp;
  int ver;
  int res;
  Mat *C, *A, *B;

  
  printf("Size/dimension of a square matix K "); res= scanf("%d", &k);
  printf("Degree or integer part of a percentage as density D "); res=scanf("%d", &D);
  printf("Number of parallel processes P "); res=scanf("%d", &Ps);


  // We build a sparse matrix with blocks 
  
  M = buildrandom_coomb_list(k, D);
  ver = validate_b(M);

  printf("\n BM %d BN %d \n", BM_, BN_) ;
  
  if (!ver) {
    return 2;
  }
      
  
  // We create a copy of the sparse matrix 
  printf("M %ld %d %d \n", M.length, M.M, M.N); 
  MT = M;
  MT.data = (COOB*) malloc(M.length*sizeof(COOB));
  assert(MT.data);
  for (int i=0;i<M.length;i++) MT.data[i] = M.data[i];



  // We sort the matrix M so that is row friendly. The sorting is
  // parallel: we need to specify the number of threads
  //print_coo(MT);
  START_CLOCK;
  rowsort_bp(M,Ps); // we really transpose the data so that the order is
  END_CLOCK;
  //print_coo_c(MT);
  printf("ROW SORT\n"); 

  //print_coomb(MT);
  START_CLOCK;
  columnsort_b(&MT); // we really transpose the data so that the order is
  END_CLOCK;
  //
  printf("COL SORT\n"); 
  //print_coomb(MT);
  ver = validate_b(M);
  if (!ver) {
    return 2;
  }
  
  
  printf("MT %ld %d %d \n", MT.length, MT.M, MT.N); 
  
  if (DEBUG) print_coomb(M);

  // We can run the M*MT either as a single computation or as a
  // pthreaded algorithm
  
  START_CLOCK;
  if (Ps<=1 ) { 
    temp = matmul_coo_b(M,MT);
  } 
  else { 
    temp = matmul_coo_par_b(M,M,MT,Ps);
  }
  END_CLOCK;
  printf("## %2d %3d %6d %ld %fMFLOPS\n",  1, D,MT.M,temp.ops,
	 temp.ops/duration/1000000);
  
  
  //print_coo(temp);


  if (temp.length*BN_*BN_<1000000)  { 
    A = build_densemb(M,1);
    
    B = build_densemb(MT,1);
    if (Ps>1) C = build_densemb(M,1);
    else      C = build_densemb(temp,0);

    
    START_CLOCK;
    matmul_f(C, temp.M*BM_, temp.N*BN_,
	     A, M.M*BM_, M.N*BN_,
	     B, MT.M*BM_, MT.N*BN_);
    END_CLOCK;
    printf("DIFF %f \n",compare_dense_mb(temp,C));
    
    free(A);
    free(B);
    free(C);
  }

  // Clean up after your self Python's coders
  free(MT.data);
  free(M.data);
  free(temp.data);

  return 0;
}

