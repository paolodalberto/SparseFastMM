
#define GRAPH_PATH 1
typedef int Mat ;
static int DEBUG = 0;

#include <SparseBLAS.h>

#define add(a,b) (((a)<(b))?(a):(b))
#define e_a  INT_MAX
#define mul(a,b) ((a)+(b))
#define e_m  0



#include <stdio.h>
#include <stdlib.h>
#include <sorting.h>

#include <pthread.h>

#define _GNU_SOURCE
#include <sched.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <assert.h>

/*
#define _GNU_SOURCE
#include <sched.h>
*/
#define CPU_SETSIZE __CPU_SETSIZE
# define CPU_ZERO(cpusetp) __CPU_ZERO_S (sizeof (cpu_set_t), cpusetp)
# define CPU_SET(cpu, cpusetp) __CPU_SET_S (cpu, sizeof (cpu_set_t), cpusetp)
# define CPU_ISSET(cpu, cpusetp) __CPU_ISSET_S (cpu, sizeof (cpu_set_t), cpusetp)
# define CPU_COUNT(cpusetp)      __CPU_COUNT_S (sizeof (cpu_set_t), cpusetp)


typedef void  (*MatrixComputation)(COO *C,  COO A, COO B);
typedef struct operands_addition TAddOperands;

struct operands_addition { 
  int  pi;
  MatrixComputation m;  // C = A*B 
  COO   *c;
  COO   a;
  COO   b;
} ;


void *basicComputation( void *s) {
  TAddOperands mc = *(TAddOperands *)s;
  int p1;
  cpu_set_t mask;
  if (mc.pi > 0)  {

    CPU_ZERO(&mask);
    CPU_SET(mc.pi, &mask);

    //p1 = sched_setaffinity(0,sizeof(mc.pi),&(mc.pi));
    p1 = sched_setaffinity(0,sizeof(mask),&(mask));
    if (p1<0) { 
      printf(" Fail processor setting pt %d \n",mc.pi);
    }
  }
  if (DEBUG) printf(" A %d x %d x %d \n",mc.c->M,mc.a.N,mc.b.N);
  mc.m(mc.c, mc.a,mc.b);
  return 0;
}


void MatrixComputations(TAddOperands *args, int len)  {
  
  pthread_t*  p_thread; /* thread's structure */
  pthread_attr_t attr;
  int* thr_id;
  int i;
  int k=len;
  
  thr_id = malloc(k * sizeof(int) );
  p_thread = malloc(k * sizeof(pthread_t) );
  pthread_attr_init(&attr);


  pthread_attr_setscope(&attr, PTHREAD_SCOPE_SYSTEM);
  
  for (i = 0; i<k-1; i++){
    //printf("k %d \n",k);
    thr_id[i] = pthread_create(&p_thread[i], 
			       &attr, 
			       basicComputation, 
			       (void *)(args+i));
  }

  basicComputation((void *)(args+i));
  
 //START_CLOCK;
 /* wait for the threads to complete */
 for (i = 0; i<k-1; i++){
   pthread_join(p_thread[i], NULL);
 }
 if (DEBUG) printf(" Done pthreading \n");

 free(thr_id);
 free(p_thread);


}




COO merge( COO C, COO T) {  // R = C+T both sparse

  COOTemporary TR= {NULL, 0, C.M, C.N};  
  COO R;
  int i,j;
  COOE c,t;
  initialize_coot(&TR);          

  
  for (i=0, j =0; i<C.length || j<T.length; ){
    c = C.data[i];
    t = T.data[j];
    if ((c.m*C.N+c.n)<(t.m*T.N+t.n)) {
      append_coot(&TR,c);
      i++;
    } else if ((c.m*C.N+c.n)>(t.m*T.N+t.n)) {
      append_coot(&TR,t);
      j++;
    } else {
      c.value = add(c.value,t.value);
      append_coot(&TR,c);
      i++;
      j++;
    }
  }
  R.length = TR.length;
  R.data = (COOE*) malloc(TR.length*sizeof(COOE)); 
  for (int t=0; t<TR.length;t++)	{
    R.data[t] = index_coot(&TR,t);
  }
  if (DEBUG) printf("Compressed  \n");
	
  free_coot(&TR);
  
  return R;
}

static inline int count_rows(COO A, int *b) {
  int count; 
  
  if (A.length<=0 ) return 0;


  count = 0;
  b[count] = 0; 
  for (int i=1; i<A.length; i++)
    if (A.data[i-1].m != A.data[i].m) {
      count ++;
      b[count] = i;
    }

  return count+1;
}




COO *split_rows(COO A, int Ps) {

  int L = A.M;
  int *b;
  COO *Rows;
  int r;
  int K, RK;
  int k=0;
  int i;

  Rows = (COO*) malloc(Ps*sizeof(COO));
  b = (int*) malloc(L*sizeof(int));

  r = count_rows(A,b);


  RK = r%Ps;
  
  K = r/(Ps) + ((RK>0)?1:0) ;

  if (DEBUG) printf("Rows = %d K =%d RK =%d Ps=%d \n",r,K, RK,Ps);     
  for (i=0; k<Ps-1;i+=K,k++) {
    if (DEBUG) printf("k < %d i =%d \n",k,i);     
    Rows[k].data = A.data +b[i];
    Rows[k].length = b[i+K]-b[i];
    Rows[k].M = A.M;
    Rows[k].N = A.N;
  }
  if (1) {
    if (DEBUG) printf("k > %d i =%d r =%d \n",k,i,r);     
    Rows[k].data = A.data +b[i];
    Rows[k].length = b[r-1]-b[i];
    Rows[k].M = A.M;
    Rows[k].N = A.N;
  }  
  free(b);
  if (DEBUG) printf("#Rows => %d \n",k);     

  
  return Rows;
}




COO matmul_coo_par(COO C,COO A,COO B,
	       int Ps /* number of threads */
	       ) {

  int i, j, k;
  COO *Rows = split_rows(A,Ps);
  COO *Ts   = (COO*) malloc(Ps*sizeof(COO)); 

  COO TR = { NULL, 0, A.M, A.N};
  COO R = { NULL, 0, A.M, A.N};
  TAddOperands *args = (TAddOperands*) malloc(Ps*sizeof(TAddOperands));
  

    // This will be parallelized 
  for (i=0;i<Ps;i++) {
    args[i].pi = i;
    args[i].m = matmul_coo_AB;
    args[i].c = Ts+i;
    args[i].a = Rows[i];
    args[i].b = B;
  }
  MatrixComputations(args,Ps);
  free(args);

  // collecting theresults 
  j = 0;
  for (k=0,i=0;i<Ps;i++) 
    j += Ts[i].length;
  
  TR.data = (COOE *) malloc(j*sizeof(COOE));
  TR.length = j;
  
  for (k=0,i=0;i<Ps;i++) 
    for (j=0; j< Ts[i].length; j++)
      TR.data[k++] = Ts[i].data[j];
  
  
  for (k=0,i=0;i<Ps;i++) 
    free(Ts[i].data);
  free(Ts);
  free(Rows);

  R = merge(C,TR);

  free(TR.data);
  
  return R;
}
  
  
  
  
  
