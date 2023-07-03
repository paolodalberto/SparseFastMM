
#ifdef EIGHT
#define BM_  8 
#define BN_  8 
#else
#define BM_  2 
#define BN_  2 
#endif


// sparse matrix is a composition of coordinate (m,n,Bm,Bn, value)
struct coo_type_block {
  int m;  // < M
  int n;  // < N
  Mat value[BM_*BN_];
} ;
typedef struct coo_type_block COOB;

#ifdef EIGHT
#define EMPTY_BLOCK  {					       \
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,		       \
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,		       \
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,		       \
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,		       \
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,		       \
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,		       \
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,		       \
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0		       \
      }							       
#else
#define EMPTY_BLOCK  {				\
    0.0, 0.0,						       \
      0.0, 0.0						       \
      }	
#endif

static COOB EA = {0, 0, EMPTY_BLOCK};

static inline int e_ab(COOB *B) {
  Mat *b = B->value;
  for (int i=0; i<BM_; i++) 
    for (int j=0; j<BN_; j++) 
      if (b[i*BN_+j] != EA.value[i*BN_+j])  return 0; 
  
  return 1;

}

static inline void add_b(COOB *C,  COOB *A , COOB *B){
  
  Mat *c = C->value; Mat *a = A->value; Mat *b = B->value;
  
  for (int i=0; i<BM_; i++) 
    for (int j=0; j<BN_; j++) 
      c[i*BN_+j] = add(a[i*BN_+j],b[i*BN_+j]);
}
static inline void mul_b(COOB *C,  COOB *A , COOB *B){
  Mat *c = C->value; Mat *a = A->value; Mat *b = B->value;
  Mat t ;
  for (int i=0; i<BM_; i++) 
    for (int j=0; j<BN_; j++){ 
      t = 0 ;

      for (int k=0; k<BM_; k++) 
	t= add(t,mul(a[i*BN_+k],b[k*BN_+j]));
      
      c[i*BN_+j] = t;
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
