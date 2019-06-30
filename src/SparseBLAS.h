//#include <Python.h>
//#include <arrayobject.h>
#include <limits.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h> 



#ifdef GRAPH_PATH
typedef int Mat ;
#define add(a,b) (((a)<(b))?(a):(b))
#define e_a  INT_MAX
#define mul(a,b) ((a)+(b))
#define e_m  0
#endif

#ifdef ALGEBRA_PATH
typedef double Mat ;
#define add(a,b) (a+b)
#define e_a  0
#define mul(a,b) ((a)*(b))
#define e_m  1
#endif


struct degree_type {
  int d;  // minimum 
  int D;  // maximum 
  int ad; // average degree
} ;



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

  
/* 
   This is to abstract the comparison function used in the sorting.
   comp(A,B,o) ~ A<B The coordinates are not enough we need to have a
   dimension so that we can express an order by row or column (M, N)
   and decreasing or increasing order.
   
*/
typedef int (*Comparing)(COOE *a, COOE *b, Ordering *c);


/* 
   this is a real matrix a sequence of coordinate elements: length
   expresses the number of elements and M and N specify the dense
   number of row and columns.
 */
struct coo_matrix {
  COOE *data;
  long unsigned int length; // number of COOE
  int M;      // Dimensions Row
  int N;      // Dimensions Column
};

typedef struct coo_matrix COO;

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




typedef struct degree_type Degree;
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
  T->data[0] = (COOE*) malloc(T->N*sizeof(COOE));
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


static inline int free_coot(COOTemporary *T) {
  long unsigned int L = T->length;
  int i = L/T->N;


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
    T->data[i] = (COOE*) malloc(T->N*sizeof(COOE));
    assert(T->data[i] );
    add = 1;
  }
  T->data[i][j] = val;
  T->length ++;
  return add;
}


#ifndef SPARSEBLASDEF
#define SPARSEBLASDEF

extern COO matmul_coo(COO C,COO A,COO B);
extern void matmul_coo_AB(COO *C, COO A,COO B);
extern COO buildrandom_coo_list(int k, int D);
extern void matmul_f(
		     Mat *C, int cm, int cn,
		     Mat *A, int am, int an,
		     Mat *B, int bm, int bn);

extern void print_coo(COO B);
extern void print_coo_c(COO B);
extern int validate(COO B);
extern int validateT(COO B);

#endif



//#include<sorting.h>



