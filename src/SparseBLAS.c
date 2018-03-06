


//#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

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



COO buildrandom_coo_list(int k, int D){
  int range = 0;
  COO res = { NULL, 0, k,k};
  COOE *array ;
  array =  (COOE *) malloc(sizeof(COOE)*k*((10*D<k)?(10*D):k));
  
  res.data = array;

  if (DEBUG)
    printf("Degree %d \n", D);



  for (int i=0; i < k; i++) 
    for (int j=0; j<k; j++ ) 
      if ((random() % k) <= D || i ==j) {
	array[range].m = i;
	array[range].n = j;
	array[range].value = 1;
	range ++;
      }
  
  res.length = range;

  if (DEBUG)
    printf("Size allocated %d and used %d \n", k*k,range);
  return res;
}



// in row format
int collectrow(COOE *array, int len, int i, int row) {
  int j;

  for (j=i; j<len && array[j].m ==row; j++);

  return j;
}

// in column format
int collectcol(COOE *array, int len, int i, int col) {

  int j;

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


int initialize_coot(COOTemporary *T) { 
  T->data = (COOE**) calloc(T->M,sizeof(COOE*));
  assert(T->data);
  for (int i=0;i<T->M;i++) {
    T->data[i] = 0;
  }
  T->data[0] = (COOE*) malloc(T->N*sizeof(COOE));
  assert(T->data[0]);
  return 1;
}

static inline COOE index_coot(COOTemporary *T, int L) {
  int i = L/T->M;
  int j = L % T->M;

  assert((L>T->length)?0:1) ;
  return T->data[i][j];
}


static inline int free_coot(COOTemporary *T) {
  int L = T->length;
  int i = L/T->N;




  

  for (;i>=0; i--) {
    free(T->data[i]);
    T->data[i] = NULL;
  }
  T-> length = 0;
  free(T->data);
  return 1;
}


int append_coot(COOTemporary *T, COOE val) {
  int L = (T->length);
  int i = L/T->M;
  int j = L % T->M;
  int add=0;
  if (DEBUG && !T->data[i])  printf(" append_coot L =%d i=%i j=%d ARR=%d \n",L, i, j, (int)T->data[i]);
  if (!T->data[i]) {
    T->data[i] = (COOE*) malloc(T->N*sizeof(COOE));
    assert(T->data[i] );
    add = 1;
  }
  T->data[i][j] = val;
  T->length ++;
  return add;
}




COO matmul_coo(COO C,COO A,COO B) {

  int i, j, t, l,row, col;
    COOTemporary T = { NULL, 0, C.M, C.N};
    initialize_coot(&T);
    COO CT = { NULL, 0, C.M, C.N }; 

    l  = 0; // C and T runner 
    i = 0; // A runner

    while (i<A.length) {
      
      int iii = collectrow(A.data, A.length,i,A.data[i].m); // from i to iii there is the A[row] vector 
      int ii=i;
      row = A.data[i].m;
      if (DEBUG) printf("i = %d Row %d ii=%d iii=%d \n",i,row,ii,iii);
      // filling entire rows from C till we have the first row of A
      for (; C.data[l].m<row ; l++)  {
	int res = append_coot(&T, C.data[l]);
	if (DEBUG) printf("\t l=%d CR append row res=%d at %d C (%d,%d,%d) \n",
			  l,res,row,C.data[l].m,C.data[l].m,C.data[l].value);
      }

      j = 0;  // B runner
      while (j<B.length) {
	COOE temp = { row, B.data[j].n, e_a}; // temporary to hold the product 
	int jjj = collectcol(B.data, B.length,j,B.data[j].n); // from j to jjj there is the B[col] vector 
	int jj=j;
	col = B.data[j].n;
	if (DEBUG) printf("%d\t j=%d Col %d jj=%d jjj=%d \n",l,j,col,jj,jjj);
	// filling the column of C
	if (DEBUG) printf("D%d C.data[%d] %d %d %d  \n",DEBUG,l,C.data[l].m,C.data[l].n,C.data[l].value);
	for (int k=l;
	     C.data[l].n<col &&
	       C.data[l].m==row ;
	     l++,k++)  {
	  int res = 	  append_coot(&T, C.data[l]);
	  if (DEBUG) printf("\t\t l=%d CC append  r=%d at=%d C (%d,%d,%d) \n",l, res,col,C.data[l].m,C.data[l].n,C.data[l].value);
	}

	// filling the element of C
	if (C.data[l].m == row && C.data[l].n == col) {
	  temp = C.data[l];
	  l++;
	  if (DEBUG) printf("\t\t temp <- C.data[l] (%d,%d,%d) \n", temp.m, temp.n,temp.value);
	}
	
	// a_row * b_col is like a merge
	while (ii<iii && jj<jjj) {
	  if (A.data[ii].n == B.data[jj].m)  {

	    temp.value = add(temp.value, mul(A.data[ii].value,B.data[jj].value));
	    if (DEBUG) printf("\t\t Merge  (%d,%d,%d) \n", temp.m, temp.n,temp.value);
	    ii ++;  jj ++;
	  }
	  else { 
	    if (A.data[ii].n < B.data[jj].m)   ii++;
	    else                               jj++;
	  }
	}
	//Done because if either is empty nothing to do e_a*w = e_a

	// if temp is not e_a (identity for +)  
	if (temp.value != e_a) { int res = append_coot(&T, temp);
	  if (DEBUG) printf("\t\t append CT %d temp (%d,%d,%d) \n",res,temp.m,temp.n,temp.value);
	}
	
	j = jjj;  // next column 
	if (DEBUG) printf("\t end jjj %d \n",jjj);
      }
	if (DEBUG) printf("end iii %d \n",iii);
      i =iii; // next row
    }
    
    // we copy the temporary result as a sparse and contiguous
    // matrix and deallocate the temporary file.
    if (DEBUG) printf("Compressing %d \n", T.length);
    CT.length = T.length;
    CT.data = (COOE*) malloc(T.length*sizeof(COOE)); 
    for (t=0; t<T.length;t++)	{
      CT.data[t] = index_coot(&T,t);
    }
    if (DEBUG) printf("Compressed  \n");
	
    free_coot(&T);
    if (DEBUG) printf("free TEMP \n");
    return CT;
    
}


void print_coo(COO B) {

  int rows =B.data[0].m;
  printf("L=%d M=%d N=%d S=%d \n",B.length,B.M, B.N, sizeof(COOE));
  for (int ktemp=0; ktemp<B.length; ktemp++) {
    if (B.data[ktemp].m!= rows) {
      printf("\n");
      rows = B.data[ktemp].m;
    }
    printf("(%d,%d,%d)", B.data[ktemp].m,B.data[ktemp].n,B.data[ktemp].value  );
  }
  printf("\n");

}

void print_coo_c(COO B) {

  int cols =B.data[0].n;
  printf("L=%d M=%d N=%d S=%d \n",B.length,B.M, B.N, sizeof(COOE));
  for (int ktemp=0; ktemp<B.length; ktemp++) {
    if (B.data[ktemp].n!= cols) {
      printf("\n");
      cols = B.data[ktemp].n;
    }
    printf("(%d,%d,%d)", B.data[ktemp].m,B.data[ktemp].n,B.data[ktemp].value  );
  }
  printf("\n");

}





void matmul_f(
	      Mat *C, int cm, int cn,
	      Mat *A, int am, int an,
	      Mat *B, int bm, int bn) {

  if (DEBUG) {
    printf("MAT C <%d, %d > ", cm,cn);
    printf("MAT A <%d, %d > ", am,an);
    printf("MAT B <%d, %d > \n", bm,bn);
    if (cm <20 && cn <20) {
      for (int i=0; i<am; i++) {
	for (int j=0; j< bn; j++) 
	  printf("%f ", C[i*cn + j]);
	    printf("\n");
      }
    }
  }
  for (int i=0; i<am; i++) 
    for (int j=0; j< bn; j++) 
      for (int k=0;k< an; k++) 
	
	C[i*cn + j] = add(C[i*cn + j],
			  mul(A[i*an + k],B[k*bn + j])
			  );
  

  if (DEBUG) {
    printf("MAT RC <%d, %d > \n", cm,cn);
    if (cm <20 && cn <20) {
      
      for (int i=0; i<am; i++) {
	for (int j=0; j< bn; j++) 
	  printf("%f ", C[i*cn + j]);
	printf("\n");
      }
    
    }
  }
  

  
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







int main() {
  int k, D;
  COO M,MT;
  COO temp;
  
  printf("K"); scanf("%d", &k);
  printf("D"); scanf("%d", &D);
  
  M = buildrandom_coo_list(k, D);

  printf("M %d %d %d \n", M.length, M.M, M.N); 
  START_CLOCK;
  MT = M;

  MT.data = (COOE*) malloc(M.length*sizeof(COOE));
  
  assert(MT.data);
  
  
  for (int i=0;i<M.length;i++) MT.data[i] = M.data[i];

  columnsort(&MT); // we really transpose the data so that the order is

  printf("MT %d %d %d \n", MT.length, MT.M, MT.N); 
  
  if (DEBUG) print_coo(M);

  temp = matmul_coo(M,M,MT);
  END_CLOCK;
  if (DEBUG) print_coo(temp);
  
  free(MT.data);
  free(M.data);
  free(temp.data);

  return 0;
}


