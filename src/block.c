
//#define GRAPH_PATH 1

static int DEBUG = 0;
static int DEBUG2=0;

//#define SPARSEBLASMBDEF 1

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


typedef COOMB  (*MatrixComputationMB)(COOMB A, COOMB B);
typedef struct operands_addition_b TAddOperandsB;

struct operands_addition_b { 
  int  pi;
  MatrixComputationMB m;  // C = A*B 
  COOMB   *c;
  COOMB   a;
  COOMB   b;
} ;


void *basicComputation( void *s) {
  TAddOperandsB mc = *(TAddOperandsB *)s;
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
  //print_coo(mc.a);
  //print_coo_c(mc.b);
  //print_coo(mc.c);
    
  return 0;
}


void MatrixComputationsB(TAddOperandsB *args, int len)  {
  
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




COOMB merge( COOMB C, COOMB T) {  // R = C+T both sparse

  COOBTemporary TR= {NULL, 0, C.M, C.N};  
  COOMB R = initialize_COOMB(NULL, 0, C.M, C.N);  
  int i,j;
  COOB c,t;
  initialize_coot_b(&TR);          
  long unsigned int ops=0;
  
  for (i=0, j =0; i<C.length && j<T.length; ){
    c = C.data[i];
    t = T.data[j];
    if ((c.m*C.N+c.n)<(t.m*T.N+t.n)) {
      append_coot_b(&TR,c);
      i++;
    } else if ((c.m*C.N+c.n)>(t.m*T.N+t.n)) {
      append_coot_b(&TR,t);
      j++;
    } else {
      ops++;
      add_b(&c, &c,&t);
      append_coot_b(&TR,c);
      i++;
      j++;
    }
  }
  for (; i<C.length ;i++ )
    append_coot_b(&TR,C.data[i]);
  
  for (;  j<T.length; j++)
    append_coot_b(&TR,T.data[j]);
  

  R.length = TR.length;
  R.data = (COOB*) malloc(TR.length*sizeof(COOB));
  R.ops = ops + T.ops;
  for (int t=0; t<TR.length;t++)	{
    R.data[t] = index_coot_b(&TR,t);
  }
  if (DEBUG) printf("Compressed  \n");
  
  free_coot_b(&TR);
  
  return R;
}


static inline int count_rows_b(COOMB A, int *b) {
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




COOMB *split_rows(COOMB A, int Ps) {

  int *b;
  COOMB *Rows;
  int L = A.M;
  int r;
  int K, RK;
  int k=0;
  int i;

  b = (int*) malloc((L+1)*sizeof(int));
  Rows = (COOMB*) malloc(Ps*sizeof(COOMB));
  
  //if (DEBUG) printf("%d\n",(int)b); 
  
  r = count_rows_b(A,b);
  
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


/***************************************
 * Sparse Matrix COO = COO * COO
 * COO  = [*COOE, nnz, ops, M, N] 
 * COOE =[ m,n, val ]
 *  
 **************************************/
// in row format
long unsigned int
collectrow_b(
	   COOB *array,
	   long unsigned int len,
	   long unsigned int i,
	   int row) {
  long unsigned int j;

  for (j=i; j<len && array[j].m ==row; j++);

  return j;
}

// in column format
long unsigned int
collectcol_b(
	   COOB *array,
	   long unsigned int len,
	   long unsigned int i,
	   int col) {

  long unsigned int j;

  for (j=i; j<len && array[j].n ==col; j++);

  return j;
}


COOMB matmul_coo_b(COOMB A,COOMB B) {

  long unsigned int i, j, t, l,row, col;
  COOBTemporary T = { NULL, 0, A.M, B.N}; 
  COOMB CT = initialize_COOMB( NULL,  0,A.M, B.N ); 
  initialize_coot_b(&T);
  long unsigned int ops = 0;
  COOB temp_m = { 0, 0, EMPTY_BLOCK};

 
  l = 0; // C and T runner 
  i = 0; // A runner

  while (i<A.length) {
    // from i to iii there is the A[row] vector 
    long unsigned int iii = collectrow_b(A.data, A.length,i,A.data[i].m); 
    long unsigned int ii=i;
    row = A.data[i].m;
    if (DEBUG) printf("i = %lu Row %lu ii=%lu iii=%lu \n",i,row,ii,iii);
    // filling entire rows from C till we have the first row of A
    // if (kk) { DEBUG=1; kk =0; }
    //else DEBUG =0;
	      
    j = 0;  // B runner
    while (j<B.length) {
      // temporary to hold the product
      COOB temp = { row, B.data[j].n, EMPTY_BLOCK}; 
      // from j to jjj there is the B[col] vector 
      long unsigned int jjj = collectcol_b(B.data, B.length,j,B.data[j].n); 
      long unsigned int jj=j;
      col = B.data[j].n;
      if (DEBUG)
	printf("%lu\t j=%lu Col %lu jj=%lu jjj=%lu (%d,%d)\n",
	       l,j,col,jj,jjj,B.data[j].m,B.data[j].n);
      ii= i;
      // a_row * b_col is like a merge
      while (ii<iii && jj<jjj) {
	if (A.data[ii].n == B.data[jj].m)  {
	  ops += 2*8*8*8;
	  mul_b(&temp_m, &A.data[ii],&B.data[jj]);
	  add_b(&temp, &temp, &temp_m);

	  if (DEBUG) {
	    printf("\nA:");print_block(A.data[ii]); 
	    printf("\nB:");print_block(B.data[jj]);
	    printf("\nM:");print_block(temp_m);
	    printf("\nC:");print_block(temp);
	  }
	  ii ++;  jj ++;
	}
	else { 
	  if (A.data[ii].n < B.data[jj].m)   ii++;
	  else                               jj++;
	}
      }
      //Done because if either is empty nothing to do e_a*w = e_a
      
      // if temp is not e_a (identity for +)  
      if (!e_ab(&temp)) {
	int res = append_coot_b(&T, temp);
	if (DEBUG)
	  printf("\t\t append CT %d temp (%d,%d,%d) \n",
		 res,temp.m,temp.n,(int)temp.value[0]);
      }
      
      j = jjj;  // next column 
      if (0 && DEBUG) printf("\t end jjj %lu \n",jjj);
    }
    if (0 && DEBUG) printf("end iii %lu \n",iii);
    i =iii; // next row
  }
    
  // we copy the temporary result as a sparse and contiguous
  // matrix and deallocate the temporary file.
  if (DEBUG) printf("Compressing %lu \n", T.length);
  CT.length = T.length;
  CT.ops = ops;
  CT.data = (COOB*) malloc(T.length*sizeof(COOB)); 
  for (t=0; t<T.length;t++)	{
    CT.data[t] = index_coot_b(&T,t);
  }

  //printf("====================================\n");
  //print_coomb(CT);
  if (DEBUG) printf("Compressed  CT %d %d %ld \n",CT.M, CT.N, CT.length);
  
  free_coot_b(&T);
  if (DEBUG) printf("free TEMP \n");
  
  if (!validate_b(CT)) {
    printf("Problems with CT\n");
  } 
  
  
  return CT;
    
}


COOMB matmul_coo_par_b(COOMB C,COOMB A,COOMB B,
	       int Ps /* number of threads */
	       ) {
  long unsigned int ops = 0;
  int i, j, k;
  COOMB *Rows;
  COOMB *Ts; 

  COOMB TR = initialize_COOMB( NULL, 0, A.M, B.N);
  COOMB R = initialize_COOMB( NULL, 0, C.M, C.N);
  TAddOperandsB *args = (TAddOperandsB*) malloc(Ps*sizeof(TAddOperandsB));
  
  if (DEBUG2) printf("Parallel %d\n",Ps);

  Rows = split_rows(A,Ps);
  if (DEBUG2) printf("Rows \n");
  Ts   = (COOMB*) calloc(Ps,sizeof(COOMB));
  // This will be parallelized 
  for (i=0;i<Ps;i++) {
    
    Ts[i].M = A.M; Ts[i].N = B.N; 
    args[i].pi = i;
    args[i].m = matmul_coo_b;
    args[i].c = Ts + i;
    args[i].a = Rows[i];
    args[i].b = B;
  }

  MatrixComputationsB(args,Ps);
  free(args);

  // collecting theresults 
  j = 0;
  for (k=0,i=0;i<Ps;i++)  {
    j += Ts[i].length;
    ops += Ts[i].ops;
    if (DEBUG) printf("Ts[%i] %d %d \n",i,Ts[i].M,Ts[i].N); 
    if (DEBUG && !validate_b(Ts[i])) {
      printf("Problems with T[%d]\n",i);
    } 
  }
  TR.data = (COOB *) malloc(j*sizeof(COOB));
  TR.length = j;
  TR.ops = ops;
  if (DEBUG2) printf("Combining %d\n",Ps);
  for (k=0,i=0;i<Ps;i++) 
    for (j=0; j< Ts[i].length; j++)
      TR.data[k++] = Ts[i].data[j];

  if (DEBUG) printf("TR %d %d %ld  \n",TR.M,TR.N,TR.length); 
  if (DEBUG && !validate_b(TR)) {
    
    printf("Problems with TR\n");
  }
  if (DEBUG) printf("C %d %d \n",C.M,C.N); 
  if (DEBUG && !validate_b(C)) {
    printf("Problems with C before merge with TR\n");
  } 
  
  for (k=0,i=0;i<Ps;i++) 
    free(Ts[i].data);
  free(Ts);
  free(Rows);
  
  if ( DEBUG2) printf("Merging C and T  %d\n",Ps);
  R = merge(C,TR);

  if (DEBUG &&  !validate_b(R)) {
    printf("Problems with R\n");
  } 
  

  
  if (DEBUG2) printf("last free \n");
  free(TR.data);
  
  return R;
}

  
  
  
  
  
