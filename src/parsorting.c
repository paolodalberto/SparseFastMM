#include<stdio.h>

static int DEBUG=0;



//#define GRAPH_PATH 1
#include <SparseBLAS.h>


#define SORT_DEFINITION_PAR 1
#include <sorting.h>
#include <pthread.h>

typedef void (*in_place_sort)    (COO arr,        Comparing comp, Ordering *order);
typedef void (*in_place_sort_par)(COO arr, int T, int Q, Comparing comp, Ordering *order);

int stages(int L, int T) {
  int K=T ;
  int j;
  for (j=1; j<T; j++)   {
    K >>=1;
    if (K==1)  break;
  }
  return j;
}    


static int ceil_(int a, int b) {
  int b1 = a/b;
  int r = (a%b>0)?1:0;
  return b1+ r;
}
static int floor_(int a, int b) {
  int b1 = a/b;
  return b1;
}





struct operands_addition_a { 
  int  pi;
  in_place_sort_par m;  // C = A*B 
  COO       c;
  Comparing cmp;
  Ordering  *order;
  int       T;
  int       Q;
} ;

typedef struct operands_addition_a TOperands;


void *bComputation( void *s) {
  TOperands mc = *(TOperands *)s;
  
  if (mc.pi >= 0)  {
    int p1;
    cpu_set_t mask;
    CPU_ZERO(&mask);
    CPU_SET(mc.pi+mc.Q, &mask);

    //p1 = sched_setaffinity(0,sizeof(mc.pi),&(mc.pi));
    p1 = sched_setaffinity(0,sizeof(mask),&(mask));
    if (p1<0) { 
      printf(" Fail processor setting pt %d \n",mc.pi);
    }
  }
  
  mc.m(mc.c,mc.T,mc.Q,mc.cmp,mc.order);
  printf("S =%2d  Ops %lu %ld  \n",
	 mc.pi+mc.Q,
	 mc.c.ops,
	 mc.c.length);

  return 0;
}


void MComputations(TOperands *args, int len)  {
  
  pthread_t*  p_thread; /* thread's structure */
  pthread_attr_t attr;
  int* thr_id;
  int i;
  int k=len;
  
  thr_id = malloc(k * sizeof(int) );
  p_thread = malloc(k * sizeof(pthread_t) );
  pthread_attr_init(&attr);


  pthread_attr_setscope(&attr, PTHREAD_SCOPE_SYSTEM);
  //printf("Tk %d \n",k);
  for (i = 0; i<k; i++){
    thr_id[i] = pthread_create(&p_thread[i], 
			       &attr, 
			       bComputation, 
			       (void *)(args+i));
  }

  // bComputation((void *)(args+i));
  
 //START_CLOCK;
 /* wait for the threads to complete */
 for (i = 0; i<k; i++){
   pthread_join(p_thread[i], NULL);
 }
 if (DEBUG) printf(" Done pthreading \n");

 free(thr_id);
 free(p_thread);


}


COO MERGECOO( COO C, COO T, Comparing comp, Ordering *order) {  // R = C+T both sparse

  COO TR= initialize_COO(NULL, 0, C.M, C.N);  
  COO R;
  int i,j,k;
  COOE c,t;

  TR.data = (COOE*) malloc((T.length+C.length)*sizeof(COO));
  
  k =0;
  for (i=0, j =0; i<C.length && j<T.length; k++){
    c = C.data[i];
    t = T.data[j];
    if (comp(&c,&t,order)<0) { //(c.m*C.N+c.n)<(t.m*T.N+t.n)) {
      TR.data[k]=c;
      i++;
    } else if (comp(&c,&t,order)>0){ // (c.m*C.N+c.n)>(t.m*T.N+t.n)) {
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
    
  TR.length = k;
  
  return TR;
}




void   sort_merge_coo  (COO   A, int T, int Q, Comparing comp, Ordering *order ) {
  
  if (T>=2) {
    TOperands P[2]; 

    COO A0 = A;
    A0.length  = ceil_(A.length,2);
    A0.data = A.data;
    COO A1 = A;
    A1.length  = floor_(A.length,2);
    A1.data = A.data + A0.length;
    P[0].pi = 0 ; P[0].m = sort_merge_coo; P[0].c = A0; P[0].cmp = comp; P[0].order = order; P[0].T = ceil_(T,2);   P[0].Q = Q;
    P[1].pi = 0;  P[1].m = sort_merge_coo; P[1].c = A1; P[1].cmp = comp; P[1].order = order; P[1].T = floor_(T,2);  P[1].Q = P[0].T+Q;

    MComputations(P,2);
    
    COO TR = MERGECOO(A0, A1, comp, order);
    for (int t=0; t<A.length;t++)	{
      A.data[t] = TR.data[t];
    }
    if (DEBUG) printf("Compressed merge all \n");
      
    free(TR.data);
  }
  else  {

    //printf("T %d Q %d\n", T,Q);
    quickSort(A.data, 0, A.length-1, comp,order);
  }
}

 
void columnsort_p(COO M, int T) {
  Ordering order = {M.M, M.N, 1 } ;
  sort_merge_coo(M,T,0, colorder, &order);
}

void rowsort_p(COO M, int T) {
  Ordering order = {M.M, M.N, 1 } ;
  sort_merge_coo(M,T,0, roworder, &order);
}
