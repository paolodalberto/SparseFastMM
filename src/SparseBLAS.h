//#include <Python.h>
//#include <arrayobject.h>
#include <limits.h>
#include <assert.h>




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
  int length; // number of COOE
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
  int length; // number of COOE 
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



//#include<sorting.h>



