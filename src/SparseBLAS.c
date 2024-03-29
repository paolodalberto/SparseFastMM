


//#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION


static int DEBUG = 0;
//#define GRAPH_PATH 1



#include <stdio.h>
#include <stdlib.h>

#include <SparseBLAS.h>



#include <sorting.h>


int validate(COO A) {
  Ordering order = {A.M, A.N, 1 } ;
  
  for (long unsigned int i=1; i< A.length; i++) { 
    COOE l = A.data[i-1];
    COOE c = A.data[i];
    if (roworder(&c,&l,&order)<0 || A.data[i].m >A.M || A.data[i].n > A.N || A.data[i].n<0||A.data[i].m<0 ) {
      printf(" %d <%d location %lu A.M %d A.N %d L = %lu Current (%d,%d,%d) < previous (%d,%d,%d ) \n",
	     c.m*A.N+c.n, l.m*A.N+l.n,
	     i,A.M, A.N, A.length,
	     c.m,c.n,(int)c.value,l.m,l.n,(int)l.value);
      return 0;
    }
    
  }
  return 1;
}
int validate_b(COOMB A) {
  Ordering order = {A.M, A.N, 1 } ;
  
  for (long unsigned int i=1; i< A.length; i++) { 
    COOB l = A.data[i-1];
    COOB c = A.data[i];
    if (roworder_b(&c,&l,&order)<0 || A.data[i].m >A.M || A.data[i].n > A.N || A.data[i].n<0||A.data[i].m<0 ) {
      printf(" %d <%d location %lu A.M %d A.N %d L = %lu Current (%d,%d,%d) < previous (%d,%d,%d ) \n",
	     c.m*A.N+c.n, l.m*A.N+l.n,
	     i,A.M, A.N, A.length,
	     c.m,c.n,(int)c.value[0],l.m,l.n,(int)l.value[0]);
      return 0;
    }
    
  }
  return 1;
}

int validateT(COO A) {
  Ordering order = {A.M, A.N, 1 } ;
  
  for (long unsigned int i=1; i< A.length; i++) { 
    COOE l = A.data[i-1];
    COOE c = A.data[i];
    if (colorder(&c, &l,&order)<0  || A.data[i].m >A.M || A.data[i].n > A.N || A.data[i].n<0||A.data[i].m<0) {
      printf(" %d < %d location %lu A.M %d A.N %d L = %lu Current (%d,%d,%d) < previous (%d,%d,%d ) \n",
	     c.m+A.M*c.n ,  l.m+A.M*l.n,
	     i, A.M, A.N, A.length,
	     c.m,c.n,(int)c.value,
	     l.m,l.n,(int)l.value);
      return 0;
    }
    
  }
  return 1;
}


COO buildrandom_coo_list(int k, int D){
  long unsigned int range = 0;
  COO res = initialize_COO( NULL, 0, k,k);
  COOE *array ;
  array =  (COOE *) malloc(sizeof(COOE)*k*((10*D<k)?(10*D):k));
  int r[D];
  int q[D];
  int min, imin=0;
  static int FFF = 0;
  
  res.data = array;

  if (DEBUG)
    printf("Degree %d \n", D);

  if (FFF) printf("\n");
  for (int i=0; i < k; i++) {
    r[0] = i;
    for (int j=1; j<D; j++ )  {
      r[j] = (r[j-1] + 1+random() % (k-1)/D) %k; 
    }
    min = k;
    for (int j=0; j<D; j++ )  {
      if (r[j] <min) { min=r[j]; imin=j;}
    }
    if (FFF) {   printf("%d %d r ",i,imin);
      for (int j=0; j<D; j++ )
	printf("%2d ",r[j]);
      printf(" \tq ");
    }
    for (int j=0; j<D; j++ ) { 
      q[j] = r[(j+imin)%D];
    }
    if (FFF)  {
      for (int j=0; j<D; j++ )
	printf("%2d ",q[j]);
      printf("\n");
    }
    for (int j=0; j<D; j++ )  {
      array[range].m = i;
      array[range].n = q[j];
      array[range].value = 1;
      range ++;
      }
  }
  res.length = range;

  if (FFF) print_coo(res);
  
  if (DEBUG)
    printf("Size allocated %d and used %lu \n", k*k,range);
  return res;
}



COOMB buildrandom_coomb_list(int k, int D){
  long unsigned int range = 0;
  COOMB res = initialize_COOMB( NULL, 0, k,k);
  COOB *array ;
  array =  (COOB *) malloc(sizeof(COOB)*k*((10*D<k)?(10*D):k));
  int r[D];
  int q[D];
  int min, imin=0;
  static int FFF = 0;
  
  res.data = array;

  if (DEBUG)
    printf("Degree %d \n", D);

  if (FFF) printf("\n");
  for (int i=0; i < k; i++) {
    r[0] = i;
    for (int j=1; j<D; j++ )  {
      r[j] = (r[j-1] + 1+random() % (k-1)/D) %k; 
    }
    min = k;
    for (int j=0; j<D; j++ )  {
      if (r[j] <min) { min=r[j]; imin=j;}
    }
    if (FFF) {   printf("%d %d r ",i,imin);
      for (int j=0; j<D; j++ )
	printf("%2d ",r[j]);
      printf(" \tq ");
    }
    for (int j=0; j<D; j++ ) { 
      q[j] = r[(j+imin)%D];
    }
    if (FFF)  {
      for (int j=0; j<D; j++ )
	printf("%2d ",q[j]);
      printf("\n");
    }
    for (int j=0; j<D; j++ )  {
      array[range].m = i;
      array[range].n = q[j];
      for (int r=0; r<BM_; r++)
	for (int c=0; c< BN_; c++)
	  array[range].value[r*BN_+c] = 1;
      range ++;
      }
  }
  res.length = range;

  if (FFF) print_coomb(res);
  
  if (DEBUG)
    printf("Size allocated %d and used %lu \n", k*k,range);
  return res;
}



// in row format
long unsigned int
collectrow(
	   COOE *array,
	   long unsigned int len,
	   long unsigned int i,
	   int row) {
  long unsigned int j;

  for (j=i; j<len && array[j].m ==row; j++);

  return j;
}

// in column format
long unsigned int
collectcol(
	   COOE *array,
	   long unsigned int len,
	   long unsigned int i,
	   int col) {

  long unsigned int j;

  for (j=i; j<len && array[j].n ==col; j++);

  return j;
}





/*  C is row-column ordered
 *  A is row-column ordered 
 *  B is column-row ordered
 */



/***
    Matrix vector product - basic version                                   
    183.equake in SPEC2000 
*/
void smvp(int nodes,
	  double ***A,int *Acol,
	  int *Aindex, double **v, double **w) {
  int i;
  int Anext, Alast, col;
  double sum0, sum1, sum2;

  for (i = 0; i < nodes; i++) {
    Anext = Aindex[i];
    Alast = Aindex[i + 1];

    sum0 = A[Anext][0][0]*v[i][0] + A[Anext][0][1]*v[i][1] + A[Anext][0][2]*v[i][2];
    sum1 = A[Anext][1][0]*v[i][0] + A[Anext][1][1]*v[i][1] + A[Anext][1][2]*v[i][2];
    sum2 = A[Anext][2][0]*v[i][0] + A[Anext][2][1]*v[i][1] + A[Anext][2][2]*v[i][2];

    Anext++;
    while (Anext < Alast) {
      col = Acol[Anext];

      sum0 += A[Anext][0][0]*v[col][0] + A[Anext][0][1]*v[col][1] + A[Anext][0][2]*v[col][2];
      sum1 += A[Anext][1][0]*v[col][0] + A[Anext][1][1]*v[col][1] + A[Anext][1][2]*v[col][2];
      sum2 += A[Anext][2][0]*v[col][0] + A[Anext][2][1]*v[col][1] + A[Anext][2][2]*v[col][2];


      w[col][0] += A[Anext][0][0]*v[i][0] + A[Anext][1][0]*v[i][1] + A[Anext][2][0]*v[i][2];
      w[col][1] += A[Anext][0][1]*v[i][0] + A[Anext][1][1]*v[i][1] + A[Anext][2][1]*v[i][2];
      w[col][2] += A[Anext][0][2]*v[i][0] + A[Anext][1][2]*v[i][1] + A[Anext][2][2]*v[i][2];
      Anext++;
    }
    w[i][0] += sum0;
    w[i][1] += sum1;
    w[i][2] += sum2;
  }
}


/***************************************
 * Sparse Matrix COO = COO * COO
 * COO  = [*COOE, nnz, ops, M, N] 
 * COOE =[ m,n, val ]
 *  
 **************************************/


COO matmul_coo(COO A,COO B) {

  long unsigned int i, j, t, l,row, col;
  COOTemporary T = { NULL, 0, A.M, B.N}; 
  COO CT = initialize_COO( NULL,  0,A.M, B.N ); 
  initialize_coot(&T);
  long unsigned int ops = 0;
 
  l = 0; // C and T runner 
  i = 0; // A runner

  while (i<A.length) {
    // from i to iii there is the A[row] vector 
    long unsigned int iii = collectrow(A.data, A.length,i,A.data[i].m); 
    long unsigned int ii=i;
    row = A.data[i].m;
    if (DEBUG) printf("i = %lu Row %lu ii=%lu iii=%lu \n",i,row,ii,iii);
    // filling entire rows from C till we have the first row of A
    // if (kk) { DEBUG=1; kk =0; }
    //else DEBUG =0;
	      
    j = 0;  // B runner
    while (j<B.length) {
      // temporary to hold the product
      COOE temp = { row, B.data[j].n, e_a}; 
      // from j to jjj there is the B[col] vector 
      long unsigned int jjj = collectcol(B.data, B.length,j,B.data[j].n); 
      long unsigned int jj=j;
      col = B.data[j].n;
      if (DEBUG)
	printf("%lu\t j=%lu Col %lu jj=%lu jjj=%lu (%d,%d)\n",
	       l,j,col,jj,jjj,B.data[j].m,B.data[j].n);
      ii= i;
      // a_row * b_col is like a merge
      while (ii<iii && jj<jjj) {
	if (A.data[ii].n == B.data[jj].m)  {
	  ops +=2;
	  temp.value = add(temp.value, mul(A.data[ii].value,B.data[jj].value));
	  if (DEBUG) printf("\t\t Merge  (%d,%d,%d) \n", temp.m, temp.n,(int)temp.value);
	  ii ++;  jj ++;
	}
	else { 
	  if (A.data[ii].n < B.data[jj].m)   ii++;
	  else                               jj++;
	}
      }
      //Done because if either is empty nothing to do e_a*w = e_a
      
      // if temp is not e_a (identity for +)  
      if (temp.value != e_a) {
	int res = append_coot(&T, temp);
	if (DEBUG)
	  printf("\t\t append CT %d temp (%d,%d,%d) \n",
		 res,temp.m,temp.n,(int)temp.value);
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
  CT.data = (COOE*) malloc(T.length*sizeof(COOE)); 
  for (t=0; t<T.length;t++)	{
    CT.data[t] = index_coot(&T,t);
  }
  if (DEBUG) printf("Compressed  CT %d %d %ld \n",CT.M, CT.N, CT.length);
  //printf("====================================\n");
  //print_coo(CT);
  
  free_coot(&T);
  if (DEBUG) printf("free TEMP \n");
  
  if (!validate(CT)) {
    printf("Problems with CT\n");
  } 
  
  
  return CT;
    
}





void print_coo(COO B) {

  int rows =B.data[0].m;
  printf("L=%lu M=%d N=%d S=%lu \n",B.length,B.M, B.N, sizeof(COOE));
  for (long unsigned int ktemp=0; ktemp<B.length; ktemp++) {
    if (B.data[ktemp].m!= rows) {
      printf("\n");
      rows = B.data[ktemp].m;
    }
    printf("(%d,%d,%d)", B.data[ktemp].m,B.data[ktemp].n,(int)B.data[ktemp].value  );
  }
  printf("\n");

}


void print_block(COOB B) {
  
  for (int r=0; r<BM_; r++)
    for (int c=0; c< BN_; c++)
      printf("(%d)", (int)B.value[r*BN_+c]);
}

void print_coomb(COOMB B) {

  int rows =B.data[0].m;
  printf("L=%lu M=%d N=%d S=%lu \n",B.length,B.M, B.N, sizeof(COOE));
  for (long unsigned int ktemp=0; ktemp<B.length; ktemp++) {
    if (B.data[ktemp].m!= rows) {
      printf("\n");
      rows = B.data[ktemp].m;
    }
    printf("(%d,%d)", B.data[ktemp].m,B.data[ktemp].n);
    print_block(B.data[ktemp]);
  }
  printf("\n");

}

void print_dense(Mat *a, int M, int N) {

  for (int r=0; r<M; r++){ 
    for (int c=0; c< N; c++)
      printf("(%d)", (int) a[r*N+c]);
    printf("\n");
  }
}


  

void print_coo_c(COO B) {

  int cols =B.data[0].n;
  printf("L=%lu M=%d N=%d S=%lu \n",B.length,B.M, B.N, sizeof(COOE));
  for (long unsigned int ktemp=0; ktemp<B.length; ktemp++) {
    if (B.data[ktemp].n!= cols) {
      printf("\n");
      cols = B.data[ktemp].n;
    }
    printf("(%d,%d,%d)", B.data[ktemp].m,B.data[ktemp].n,(int)B.data[ktemp].value  );
  }
  printf("\n");

}


Mat *build_dense(COO A, int def) {
  
  Mat *c =  (Mat *) calloc(A.M*A.N,sizeof(Mat));
  if (def==0) return c;
  
  for (int l=0; l<A.length; l++) {
    COOE d = A.data[l];
    c[d.m*A.N + d.n] = d.value*def;
  }

  return c;
}
Mat *build_densemb(COOMB A, int def) {
  
  Mat *c =  (Mat *) calloc(A.M*A.N*BN_*BM_,sizeof(Mat));
  if (def==0) return c;
  
  for (int l=0; l<A.length; l++) {
    COOB d = A.data[l];
    Mat * C = c+(d.m)*A.N*BN_*BM_ +d.n*BN_;

    for (int r=0; r<BM_; r++)
      for (int cc=0; cc< BN_; cc++)
	
	C[r*A.N*BN_ + cc] = d.value[r*BN_+cc]*def;
    //print_dense(c,A.M*BM_, A.N*BN_);
  }

  return c;
}

double C(Mat *A, Mat *B, int M, int N) {

  double res = 0;
  for (int i=0; i<M; i++)
    for (int j=0; j<N;j++) {
#if (COMPLEX_FLOAT_ || COMPLEX_DOUBLE_)
      double dif = creal((A[i*N+j]-B[i*N+j])*(A[i*N+j]-B[i*N+j]));
#else
      double dif = (A[i*N+j]-B[i*N+j])*(A[i*N+j]-B[i*N+j]);
#endif
      if (dif) printf("(%d,%d,%f)\n", i,j,dif);
      res += A[i*N+j]-B[i*N+j];
    }
  return res;
}

double compare_dense_(void *B, Mat *def, int alg) {
  double res = 0; 
  if (alg==0) {
    COO BB = *(COO*) B;
    Mat *BM = build_dense( BB,1);
    res= C(BM,def, BB.M, BB.N);
    free(BM);
    
  }
  else if (alg==1) {
    COOMB BB = *(COOMB*) B;
    Mat *BM = build_densemb(BB,1);
    res= C(BM,def, BB.M, BB.N);
    free(BM);
  }
  else { res = 1;}
  return res ;
}



double compare_dense(COO B, Mat *def) {


  double res = 0;
  Mat diff_m = 0;
  double dif;
  int cols =B.data[0].m;
  if (B.length<10) printf("L=%lu M=%d N=%d S=%lu \n",B.length,B.M, B.N, sizeof(COOE));
  for (long unsigned int ktemp=0; ktemp<B.length; ktemp++) {
    if (B.data[ktemp].m!= cols) {
      if (B.length< 100) printf("\n");
      cols = B.data[ktemp].m;
    }
    diff_m = (B.data[ktemp].value-def[B.data[ktemp].m*B.N+B.data[ktemp].n])*(B.data[ktemp].value-def[B.data[ktemp].m*B.N+B.data[ktemp].n]);
#if (COMPLEX_FLOAT_ || COMPLEX_DOUBLE_)
    dif = creal(diff_m);
#else
    dif = (diff_m);
#endif
    
    res += dif;
    if (B.length<10)
      printf("(%d,%d,%f)", B.data[ktemp].m,B.data[ktemp].n,
	     dif );
  }
  if (B.length<10) printf("\n");
  return res;
}
double compare_dense_mb(COOMB B, Mat *def) {

  double dif;
  double res = 0;
  int cols =B.data[0].m;
  if (B.length<10) printf("L=%lu M=%d N=%d S=%lu \n",B.length,B.M, B.N, sizeof(COOE));
  for (long unsigned int ktemp=0; ktemp<B.length; ktemp++) {
    double resb = 0;
    COOB d = B.data[ktemp];
    Mat * C = def+(d.m)*B.N*BN_*BM_ +d.n*BN_;
    if (d.m!= cols) {
      if (B.length< 10) printf("\n");
      cols = d.m;
    }
    

    for (int r=0; r<BM_; r++)
      for (int c=0; c< BN_; c++) {
#if (COMPLEX_FLOAT_ || COMPLEX_DOUBLE_)
	dif = creal((d.value[r*BN_+c]-  C[r*B.N*BN_+c])*(d.value[r*BN_+c]-  C[r*B.N*BN_+c]));
#else
	dif = (d.value[r*BN_+c]-  C[r*B.N*BN_+c])*(d.value[r*BN_+c]-  C[r*B.N*BN_+c]);
#endif
	resb += dif;
      }
    if (B.length<10) printf("(%d,%d,%f)", d.m,d.n, resb );
    res += resb;
  }
  if (B.length<10) printf("\n");
  return res;
}





void matmul_f(
	      Mat *C, int cm, int cn,
	      Mat *A, int am, int an,
	      Mat *B, int bm, int bn, float beta) {

  if (DEBUG) {
    printf("MAT C <%d, %d > ", cm,cn);
    printf("MAT A <%d, %d > ", am,an);
    printf("MAT B <%d, %d > \n", bm,bn);
  }




  for (int i=0; i<am; i++) 
    for (int j=0; j< bn; j++) 
      for (int k=0;k< an; k++) 
	
	C[i*cn + j] = add(C[i*cn + j],
			  mul(A[i*an + k],B[k*bn + j])
			  );
  

  

  
}


#include <cblas.h>


void matmul_f_par(
		  Mat *C, int cm, int cn,
		  Mat *A, int am, int an,
		  Mat *B, int bm, int bn, float beta) {
  
  //void cblas_sgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, const int M, const int N, const int K, const float alpha, const float *A, const int lda, const float *B, const int ldb, const float beta, float *C, const int ldc)

#if (COMPLEX_FLOAT_ || COMPLEX_DOUBLE_)
  Mat one = e_m;
  Mat zero = e_a;
#endif

  
#if (FLOAT_) 
  cblas_sgemm(CblasRowMajor,CblasNoTrans, CblasNoTrans,
	      am, bn, an, 1.0, A, an, B, bn, beta, C, cn);
#endif
#if (DOUBLE_) 
  cblas_dgemm(CblasRowMajor,CblasNoTrans, CblasNoTrans,
	      am, bn, an, 1.0, A, an, B, bn, beta, C, cn);

#endif
#if (COMPLEX_FLOAT_) 
  if (beta==1.0) zero = one;
  cblas_cgemm(CblasRowMajor,CblasNoTrans, CblasNoTrans,
	      am, bn, an, &one, A, an, B, bn, &zero, C, cn);

#endif
#if (COMPLEX_DOUBLE_) 
  if (beta==1.0) zero = one;
  cblas_zgemm(CblasRowMajor,CblasNoTrans, CblasNoTrans,
	      am, bn, an, &one, A, an, B, bn, &zero, C, cn);

#endif
  
}


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






