

#define BM_  8 
#define BN_  8 




// sparse matrix is a composition of coordinate (m,n,Bm,Bn, value)
struct coo_type_block {
  int m;  // < M
  int n;  // < N
  Mat value[BM_*BN_];
} ;
typedef struct coo_type_block COOB;


#define EMPTY_BLOCK  { \
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,		       \
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,		       \
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,		       \
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,		       \
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,		       \
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,		       \
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,		       \
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0		       \
  }							       \


    
static inline void add_b(COOB C,  COOB A , COOB B){
  
  Mat *c = C.value; Mat *a = A.value; Mat *b = B.value;

  for (int i=0; i<BM_; i++) { 
    c[i*BN_+0] = add(a[i*BN_+0],b[i*BN_+0]);
    c[i*BN_+1] = add(a[i*BN_+1],b[i*BN_+1]);
    c[i*BN_+2] = add(a[i*BN_+2],b[i*BN_+2]);
    c[i*BN_+3] = add(a[i*BN_+3],b[i*BN_+3]);
    c[i*BN_+4] = add(a[i*BN_+4],b[i*BN_+4]);
    c[i*BN_+5] = add(a[i*BN_+5],b[i*BN_+5]);
    c[i*BN_+6] = add(a[i*BN_+6],b[i*BN_+6]);
    c[i*BN_+7] = add(a[i*BN_+7],b[i*BN_+7]);
  }
  
  
}
static inline void mul_b(COOB C,  COOB A , COOB B){
  Mat *c = C.value; Mat *a = A.value; Mat *b = B.value;

  Mat a0,a1,a2,a3,a4,a5,a6,a7;
  Mat t0,t1,t2,t3,t4,t5,t6,t7; 
  for (int i=0; i<BM_; i++){
    a0=a[i*BN_+0];a1=a[i*BN_+1];a2=a[i*BN_+2];a3=a[i*BN_+3];
    a4=a[i*BN_+4];a5=a[i*BN_+5];a6=a[i*BN_+6];a7=a[i*BN_+7]; 
    for (int j=0; j<BN_; j++) {
      t0=mul(a0,b[0*BN_+j]);
      t1=mul(a1,b[1*BN_+j]);
      t2=mul(a2,b[2*BN_+j]);
      t3=mul(a3,b[3*BN_+j]);
      t4=mul(a4,b[4*BN_+j]);
      t5=mul(a5,b[5*BN_+j]);
      t6=mul(a6,b[6*BN_+j]);
      t7=mul(a7,b[7*BN_+j]); 
      c[i*BN_+j] = t0+ t1 +t2 +t3+t4+t5+t6+t7;
    }
  }

}


  
/* 
   Block based sparse matrix .
 */
struct coo_matrix_block {
  COOB *data;
  long unsigned int length; // number of COOB
  long unsigned int ops;    // number of operations
  int M;      // Dimensions Row
  int N;      // Dimensions Column
};

typedef struct coo_matrix_block COOMB;





static inline COOMB initialize_COOMB(COOB *e, long unsigned int L, int M, int N) {
  COOMB C = { e, L, 0, M, N}; 
  return C;
}

/*
  In C we must allocate space in advance. We create chunks of memory
  that can be used and then solidified into a single block
 */
struct cootemp_matrix_b {
  COOB **data;   // data[(L/(3*M)][(L/3 % N)*3]
  long unsigned int length; // number of COOE 
  int M;      // Dimensions Row
  int N;      // Dimensions Column
};
typedef struct cootemp_matrix_b COOBTemporary;

static int initialize_coot_b(COOBTemporary *T) { 
  T->data = (COOB**) calloc(T->M,sizeof(COOB*));
  assert(T->data);
  for (int i=0;i<T->M;i++) {
    T->data[i] = 0;
  }
  T->data[0] = (COOB*) malloc(T->M*sizeof(COOB));
  assert(T->data[0]);
  return 1;
}
static inline int free_coot_b(COOBTemporary *T) {
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

static inline int append_coot_b(COOBTemporary *T, COOB val) {
  long unsigned int L = (T->length);
  int i = L/T->M;
  int j = L % T->M;
  int add=0;
  if (DEBUG && !T->data[i])  printf(" append_coot L =%lu i=%i j=%d  \n",L, i, j);
  if (!T->data[i]) {
    T->data[i] = (COOB*) malloc(T->M*sizeof(COOB));
    assert(T->data[i] );
    add = 1;
  }
  T->data[i][j] = val;
  T->length ++;
  return add;
}
