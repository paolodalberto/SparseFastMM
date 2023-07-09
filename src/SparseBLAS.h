//#include <Python.h>
//#include <arrayobject.h>
#include <limits.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#ifdef GRAPH_PATH
typedef int Mat ;
#define add(a,b) (((a)<(b))?(a):(b))
#define add(a,b) (((a)<(b))?(a):(b))
#define e_a  INT_MAX
#define mul(a,b) ((a)+(b))
#define e_m  0
#else 
typedef float Mat ;
#define add(a,b) (a+b)
#define e_a  0
#define mul(a,b) ((a)*(b))
#define e_m  1
#endif
#include <block.h>




/*
#define _GNU_SOURCE
#include <sched.h>
*/
#define CPU_SETSIZE __CPU_SETSIZE
# define CPU_ZERO(cpusetp) __CPU_ZERO_S (sizeof (cpu_set_t), cpusetp)
# define CPU_SET(cpu, cpusetp) __CPU_SET_S (cpu, sizeof (cpu_set_t), cpusetp)
# define CPU_ISSET(cpu, cpusetp) __CPU_ISSET_S (cpu, sizeof (cpu_set_t), cpusetp)
# define CPU_COUNT(cpusetp)      __CPU_COUNT_S (sizeof (cpu_set_t), cpusetp)









struct degree_type {
  int d;  // minimum 
  int D;  // maximum 
  int ad; // average degree
} ;
typedef struct degree_type Degree;




// sparse matrix is a composition of coordinate (m,n,value)
struct coo_type {
  int m;  // < M
  int n;  // < N
  Mat value;
} ;
typedef struct coo_type COOE;








struct ordering {
  int M;
  int N;
  int order;
};
typedef struct ordering Ordering;



static COOE* from_three_to_one(int *x, int *y, Mat *v, long unsigned int len) {

  COOE * t = (COOE*) calloc(len,sizeof(COOE));
  for (long unsigned int i=0; i<len; i++) {
    t[i].n =x[i];
    t[i].m =y[i];
    t[i].value = v[i];
  }
  return t;

}



  
/* 
   This is to abstract the comparison function used in the sorting.
   comp(A,B,o) ~ A<B The coordinates are not enough we need to have a
   dimension so that we can express an order by row or column (M, N)
   and decreasing or increasing order.
   
*/
typedef int (*Comparing)(COOE *a, COOE *b, Ordering *c);
typedef int (*ComparingB)(COOB *a, COOB *b, Ordering *c);

/* 
   this is a real matrix a sequence of coordinate elements: length
   expresses the number of elements and M and N specify the dense
   number of row and columns.
 */
struct coo_matrix {
  COOE *data;
  long unsigned int length; // number of COOE
  long unsigned int ops;    // number of operations
  int M;      // Dimensions Row
  int N;      // Dimensions Column
  
};

typedef struct coo_matrix COO;
static COO initialize_COO(COOE *e, long unsigned int L, int M, int N) {
  COO C = { e, L, 0, M, N}; 
  return C;
}


/*
  In C we must allocate space in advance. We create chunks of memory
  that can be used and then solidified into a single block
 */
struct cootemp_matrix {
  COOE **data;   // data[(L/(3*M)][(L/3 % N)*3]
  long unsigned int length; // number of COOE 
  int M;      // Dimensions Row
  int N;      // Dimensions Column
};
typedef struct cootemp_matrix COOTemporary;






struct adj_matrix {
  Mat *data;
  int m;
  int n;
};
typedef struct adj_matrix Adjacent ;



static int initialize_coot(COOTemporary *T) { 
  T->data = (COOE**) calloc(T->M,sizeof(COOE*));
  assert(T->data);
  for (int i=0;i<T->M;i++) {
    T->data[i] = 0;
  }
  T->data[0] = (COOE*) malloc(T->M*sizeof(COOE));
  assert(T->data[0]);
  return 1;
}



static inline COOE index_coot(COOTemporary *T, long unsigned int L) {
  int i = L/T->M;
  int j = L % T->M;
  
  if (DEBUG && L<T->length) printf("T length %lu i %d j %d \n", T->length,i,j);
  //assert((L>T->length)?1:0) ;
  return T->data[i][j];
}
static inline COOB index_coot_b(COOBTemporary *T, long unsigned int L) {
  int i = L/T->M;
  int j = L % T->M;
  
  if (DEBUG && L<T->length) printf("T length %lu i %d j %d \n", T->length,i,j);
  //assert((L>T->length)?1:0) ;
  return T->data[i][j];
}


static inline int free_coot(COOTemporary *T) {
  long unsigned int L = T->length;
  int i = L/T->M;


  for (;i>=0; i--) {
    free(T->data[i]);
    T->data[i] = NULL;
  }
  T-> length = 0;
  free(T->data);
  return 1;
}


static inline int append_coot(COOTemporary *T, COOE val) {
  long unsigned int L = (T->length);
  int i = L/T->M;
  int j = L % T->M;
  int add=0;
  if (DEBUG && !T->data[i])  printf(" append_coot L =%lu i=%i j=%d  \n",L, i, j);
  if (!T->data[i]) {
    T->data[i] = (COOE*) malloc(T->M*sizeof(COOE));
    assert(T->data[i] );
    add = 1;
  }
  T->data[i][j] = val;
  T->length ++;
  return add;
}



#ifndef SPARSEBLASDEF
#define SPARSEBLASDEF

extern COO matmul_coo(COO A,COO B);
extern COOMB matmul_coo_b(COOMB A,COOMB B);
extern COO buildrandom_coo_list(int k, int D);
extern COOMB buildrandom_coomb_list(int k, int D);
extern void matmul_f(
		     Mat *C, int cm, int cn,
		     Mat *A, int am, int an,
		     Mat *B, int bm, int bn);

extern void print_block(COOB B);
extern void print_coo(COO B);
extern void print_coomb(COOMB B);
extern void print_coo_c(COO B);
extern void print_dense(Mat *A, int M, int N);
extern int validate(COO B);
extern int validate_b(COOMB B);
extern int validateT(COO B);
extern Mat *build_dense(COO A, int def);
extern Mat *build_densemb(COOMB A, int def);
extern double compare_dense(COO A, Mat *d);
extern double compare_dense_mb(COOMB A, Mat *d);


#endif



//#include<sorting.h>



