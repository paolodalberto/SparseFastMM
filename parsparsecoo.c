
//#define GRAPH_PATH 1

static int DEBUG = 0;
static int DEBUG2=0;
#include <SparseBLAS.h>



 
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


typedef COO  (*MatrixComputation)(COO A, COO B);
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
  if (mc.pi >= 0)  {

    CPU_ZERO(&mask);
    CPU_SET(mc.pi, &mask);

    //p1 = sched_setaffinity(0,sizeof(mc.pi),&(mc.pi));
    p1 = sched_setaffinity(0,sizeof(mask),&(mask));
    if (p1<0) { 
      printf(" Fail processor setting pt %d \n",mc.pi);
    }
  }
  
  *mc.c = mc.m(mc.a,mc.b);
  printf("C =%2d  Ops %lu %ld x %ld x %ld \n",
	 mc.pi,
	 mc.c->ops,
	 mc.c->length,
	 mc.a.length,
	 mc.b.length);
  print_coo(mc.a);
  print_coo_c(mc.b);
  print_coo(mc.c);
    
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
 if (DEBUG2) printf(" Done pthreading \n");

 free(thr_id);
 free(p_thread);


}




COO merge( COO C, COO T) {  // R = C+T both sparse

  COOTemporary TR= {NULL, 0, C.M, C.N};  
  COO R = initialize_COO(NULL, 0, C.M, C.N);  
  int i,j;
  COOE c,t;
  initialize_coot(&TR);          
  long unsigned int ops=0;
  
  for (i=0, j =0; i<C.length && j<T.length; ){
    c = C.data[i];
    t = T.data[j];
    if ((c.m*C.N+c.n)<(t.m*T.N+t.n)) {
      append_coot(&TR,c);
      i++;
    } else if ((c.m*C.N+c.n)>(t.m*T.N+t.n)) {
      append_coot(&TR,t);
      j++;
    } else {
      ops++;
      c.value = add(c.value,t.value);
      append_coot(&TR,c);
      i++;
      j++;
    }
  }
  for (; i<C.length ;i++ )
    append_coot(&TR,C.data[i]);
  
  for (;  j<T.length; j++)
    append_coot(&TR,T.data[j]);
  

  R.length = TR.length;
  R.data = (COOE*) malloc(TR.length*sizeof(COOE));
  R.ops = ops + T.ops;
  for (int t=0; t<TR.length;t++)	{
    R.data[t] = index_coot(&TR,t);
  }
  if (DEBUG) printf("Compressed  \n");
  
  free_coot(&TR);
  
  return R;
}

COO merge_alt( COO C, COO T) {  // R = C+T both sparse

  COO TR= initialize_COO(NULL, 0, C.M, C.N);  
  COO R;
  int i,j,k;
  COOE c,t;

  TR.data = (COOE*) malloc((T.length+C.length)*sizeof(COO));
  
  k =0;
  for (i=0, j =0; i<C.length && j<T.length; k++){
    c = C.data[i];
    t = T.data[j];
    if ((c.m*C.N+c.n)<(t.m*T.N+t.n)) {
      TR.data[k]=c;
      i++;
    } else if ((c.m*C.N+c.n)>(t.m*T.N+t.n)) {
      TR.data[k]=t;
      j++;
    } else {
      c.value = add(c.value,t.value);
      TR.data[k]=c;
      i++;
      j++;
    }
  }
  for (; i<C.length ;i++,k++ )
    TR.data[k] = C.data[i];
  
  for (;  j<T.length; j++,k++)
    TR.data[k] = T.data[j];
  
  
  R.length = k;
  R.data = (COOE*) malloc(TR.length*sizeof(COOE)); 
  for (int t=0; t<k;t++)	{
    R.data[t] = TR.data[t];
  }
  if (DEBUG) printf("Compressed merge all \n");
	
  free(TR.data);
  
  return R;
}

static inline int count_rows(COO A, int *b) {
  int count; 
  int i;
  if (A.length<=0 ) return 0;


  count = 0;
  b[count] = 0; 
  for (i=1; i<A.length; i++)
    if (A.data[i-1].m != A.data[i].m) {
      count ++;
      b[count] = i;
    }
  count++;
  b[count] = i;
  return count;
}




COO *split_rows(COO A, int Ps) {

  int *b;
  COO *Rows;
  int L = A.M;
  int r;
  int K, RK;
  int k=0;
  int i;

  b = (int*) malloc((L+1)*sizeof(int));
  Rows = (COO*) malloc(Ps*sizeof(COO));
  
  //if (DEBUG) printf("%d\n",(int)b); 
  
  r = count_rows(A,b);
  
  if (DEBUG) printf("L = %d-%dx%d r =%d \n",L,A.M,A.N,r); 
  RK = r%Ps;
  
  K = r/(Ps) + ((RK>0)?1:0) ;

  if (DEBUG) printf("Rows = %d K =%d RK =%d Ps=%d \n",r,K, RK,Ps);     
  for (i=0; k<Ps-1;i+=K,k++) {
    if (DEBUG) printf("k < %d i =%d b[i] = %d L = %d \n",k,i,b[i],b[i+K]-b[i]-1);     
    Rows[k].data = A.data +b[i];
    Rows[k].length = b[i+K]-b[i];
    Rows[k].M = A.M;
    Rows[k].N = A.N;
  }
  if (1) {
    if (DEBUG) printf("k = %d i =%d b[i] = %d r =%d L = %d\n",k,i,b[i],r,b[r-1]-b[i]);     
    Rows[k].data = A.data +b[i];
    Rows[k].length = b[r]-b[i];
    Rows[k].M = A.M;
    Rows[k].N = A.N;
  }  
  if (DEBUG) printf("#Rows => %d \n",k);
  //printf("%u\n",(unsigned int)b); 
  free(b);


  
  return Rows;
}




COO matmul_coo_par(COO C,COO A,COO B,
	       int Ps /* number of threads */
	       ) {
  long unsigned int ops = 0;
  int i, j, k;
  COO *Rows;
  COO *Ts; 

  COO TR = initialize_COO( NULL, 0, A.M, B.N);
  COO R = initialize_COO( NULL, 0, C.M, C.N);
  TAddOperands *args = (TAddOperands*) malloc(Ps*sizeof(TAddOperands));
  
  if (DEBUG2) printf("Parallel %d\n",Ps);

  Rows = split_rows(A,Ps);
  if (DEBUG2) printf("Rows \n");
  Ts   = (COO*) calloc(Ps,sizeof(COO));
  // This will be parallelized 
  for (i=0;i<Ps;i++) {
    
    Ts[i].M = A.M; Ts[i].N = B.N; 
    args[i].pi = i;
    args[i].m = matmul_coo;
    args[i].c = Ts + i;
    args[i].a = Rows[i];
    args[i].b = B;
  }

  MatrixComputations(args,Ps);
  free(args);

  // collecting theresults 
  j = 0;
  for (k=0,i=0;i<Ps;i++)  {
    j += Ts[i].length;
    ops += Ts[i].ops;
    if (DEBUG) printf("Ts[%i] %d %d \n",i,Ts[i].M,Ts[i].N); 
    if (DEBUG && !validate(Ts[i])) {
      printf("Problems with T[%d]\n",i);
    } 
  }
  TR.data = (COOE *) malloc(j*sizeof(COOE));
  TR.length = j;
  TR.ops = ops;
  if (DEBUG2) printf("Combining %d\n",Ps);
  for (k=0,i=0;i<Ps;i++) 
    for (j=0; j< Ts[i].length; j++)
      TR.data[k++] = Ts[i].data[j];

  if (DEBUG) printf("TR %d %d %ld  \n",TR.M,TR.N,TR.length); 
  if (DEBUG && !validate(TR)) {
    
    printf("Problems with TR\n");
  }
  if (DEBUG) printf("C %d %d \n",C.M,C.N); 
  if (DEBUG && !validate(C)) {
    printf("Problems with C before merge with TR\n");
  } 
  
  for (k=0,i=0;i<Ps;i++) 
    free(Ts[i].data);
  free(Ts);
  free(Rows);
  
  if ( DEBUG2) printf("Merging C and T  %d\n",Ps);
  R = merge(C,TR);

  if (DEBUG &&  !validate(R)) {
    printf("Problems with R\n");
  } 
  

  
  if (DEBUG2) printf("last free \n");
  free(TR.data);
  
  return R;
}

COO
matmul_coo_par_basic(
		     int *_CX, int *_CY, Mat *_CV,
		     long unsigned int LC,int MC, int NC,
		     int *_AX,int *_AY,  Mat *_AV,
		     long unsigned int LA, int MA,int NA,
		     int *_BX, int *_BY, Mat *_BV,
		     long unsigned int LB, int MB, int NB,
		     int Ps /* number of threads */
		     ) {
  printf("matmul_coo_par_basic \n");
  if (DEBUG2) printf("C %lu \n", LC);
  COOE *_C = from_three_to_one(_CX,_CY,_CV,LC);
  if (DEBUG2) printf("A %lu \n", LA);
  COOE *_A = from_three_to_one(_AX,_AY,_AV,LA);;
  if (DEBUG2) printf("B %lu \n",LB);
  COOE *_B = from_three_to_one(_BX,_BY,_BV,LB);; 
  COO C = initialize_COO( _C, LC, MC, NC);
  COO A = initialize_COO( _A, LA, MA, NA);
  COO B = initialize_COO( _B, LB, MB, NB);
  COO R;

  columnsort(&B); // we really transpose the data so that the order is  


  if (!validate(C)) {
    printf("Problems with C\n");
  } 
  if (!validate(A)) {
    printf("Problems with A\n");
  }
  if (!validateT(B)) {
    printf("Problems with B\n");
  }
  
  R = matmul_coo_par(C,A,B,Ps);

  if (!validate(R)) {
    printf("Problems with R\n");
  } 


  if (DEBUG2) printf("Done matmul_coo_par\n");
  free(_C);
  free(_A);
  free(_B);
  // This has to be free by who calls this 
  return R;
}
  
  
  
  
  
