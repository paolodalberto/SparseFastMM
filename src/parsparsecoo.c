
typedef int Mat ;
#include <SparseBLAS.h>

#define add(a,b) (((a)<(b))?(a):(b))
#define e_a  INT_MAX
#define mul(a,b) ((a)+(b))
#define e_m  0

static int DEBUG = 0;

#include <stdio.h>
#include <stdlib.h>
#include <sorting.h>




COO merge( COO C, COO T) {  // R = C+T both sparse

  COOTemporary T= {NULL, 0, C.M, C.N};  
  COO R; 
  return R;
}

inline int count_rows(COO A, int *b) {
  int row;
  int count 
  
  if (A.length<=0 ) return 0;

  row  = A.data[0].m;
  count = 0;
  b[count] = 0; 
  for (int i=1; i<A.length; i++)
    if (A.data[i-1].m != A.data[i].m) {
      count ++;
      b[count] = i;
    }

  return count;
}




COO *split_rows(COO A, int Ps) {

  int L = A.M;
  int *b;
  COO *Rows;
  int r;
  int K, RK;
  int k=0;
  int i,j;

  Rows = (COO*) malloc(Ps*sizeof(COO));
  b = (int*) malloc(L*sizeof(int));

  r = count_rows(A,b);
  
  K = r/Ps;
  RK = r%Ps;
  
  
  for (i=0; i<r-K;i+=K,k++) {
    Rows[k] = {A.data +b[i], b[i+K]-b[i],A.M, A.N};
  }

  Rows[k] = {A.data +b[i], b[r-1]-b[i],A.M, A.N};
  
  
  free(b);
  
  return Rows;
}








COO matmul_coo(COO C,COO A,COO B,
	       int Ps /* number of threads */
	       ) {

  int i, j, t, l,row, col;
  COO *Rows = split_rows(A,Ps);
  COO *Ts   = (COO*) malloc(P*sizeof(COO)); 

  COO TR = { NULL, 0, A.M, A.N};
  COO R = { NULL, 0, A.M, A.N};


  // This will be parallelized 
  for (i=0;i<Ps;i++) {
    Ts[i] = matmul_coo_AB(Rows[i],B);
    j += Ts[i].length;
  }

  TR.data = (COOE *) malloc(j*sizeof(COOE));
  TR.lemgth = j;
  
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
  
  
  
  
  
